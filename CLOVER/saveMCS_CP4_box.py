import numpy as np
from scipy.ndimage.measurements import label
import xarray as xr
import os
import ipdb
import glob
from scipy.interpolate import griddata
import pandas as pd
import ipdb
import itertools
from collections import OrderedDict
from utils import constants as cnst


def olr_to_bt(olr):
    sigma = 5.670373e-8
    return ((olr/sigma)**0.25)-273.15

def griddata_lin(data, x, y, new_x, new_y):

    """
    :param x: current x variables (1 or 2d, definitely 2d if irregular!)
    :param y: current y variables (1 or 2d, definitely 2d if irregular!)
    :param new_x: target x vars
    :param new_y: target y vars
    :return:  triangulisation lookup table, point weights, 2d shape - inputs for interpolation func
    """

    if x.ndim == 1:
        grid_xs, grid_ys = np.meshgrid(x, y)
    else:
        grid_xs = x
        grid_ys = y

    if new_x.ndim == 1:
        new_xs, new_ys = np.meshgrid(new_x, new_y)
    else:
        new_xs = new_x
        new_ys = new_y

    points = np.array((grid_xs.flatten(), grid_ys.flatten())).T
    inter = np.array((np.ravel(new_xs), np.ravel(new_ys))).T
    shape = new_xs.shape

    # Interpolate using delaunay triangularization
    data = griddata(points, data.flatten(), inter, method='linear')
    data = data.reshape((shape[0], shape[1]))

    return data


def file_save(cp_dir, out_dir, ancils_dir, vars, datestring, sbox, tthresh):

    keys = vars.keys()

    if 'lw_out_PBLtop' not in keys:
        print('please provide ORL first in dictionary')
        return
    box = [-17, 14, 3.5, 21]

    #load seamask
    landsea_path = glob.glob(ancils_dir+os.sep+'landseamask*.nc')[0]
    landsea = xr.open_dataset(landsea_path, decode_times=False)
    ls = landsea['lsm']

    ls.rlon.values = ls.rlon.values-360
    ls_arr = ls.sel(rlon=slice(box[0], box[1]), rlat=slice(box[2], box[3]))
    pos = np.where(ls_arr[0, 0, :, :] == 0)
    lons, lats = np.meshgrid(ls_arr.rlon.values, ls_arr.rlat.values)
    goodinds = 0

    #create empty dataset
    ds = xr.Dataset()
    # create empty

    #loop through every var
    for v in keys:

        h = (vars[v])[1]
        pl = (vars[v])[0]
        derived = False
        if (v == 'shear') | (v == 'u_mid') | (v == 'u_srfc') | (v == 'theta'):
            if v == 'theta':
                derived = v
                v = 't_pl'

            else:
                derived = v
                v = 'u_pl'

        # try:
        #     filepath = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*'+datestring+'*.nc')[0]
        # except IndexError:
        #     print('No file found, return')
        #     return

        filepath = cp_dir+os.sep+str(v)+os.sep+'*'+str(d['time.year'].values)+str(d['time.month'].values).zfill(2)+'*.nc'

        arr = xr.open_mfdataset(filepath, autoclose=True)

        dar = arr[v].sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
        pdt = pd.to_datetime(datestring.values)
        pdt = pdt.replace(hour=h)
        dar = dar.sel(time=pdt, method='nearest')

        if int(dar['time.hour'])!=h:
            print('Wrong hour')
            return

        if 'pressure' in dar.coords:
            dar.values[dar.values==0] = np.nan # potential missing value maskout
            if len(pl) > 1:
                shear = dar.sel(pressure=pl[0]).values - dar.sel(pressure=pl[1]).values
                dar = dar.sum(dim='pressure').squeeze()
                dar.values = shear
                if derived:
                    v = derived

            else:
                dar = dar.sel(pressure=pl[0]).squeeze()
            if derived:
                v = derived

        # regrid to common grid (unstagger wind, bring to landsea mask grid)
        regrid = griddata_lin(dar.values, dar.longitude, dar.latitude, ls_arr.rlon, ls_arr.rlat)
        da = xr.DataArray(regrid,
                          coords={'time': dar.time, 'latitude': ls_arr.rlat.values,
                                  'longitude': ls_arr.rlon.values, },
                          dims=['latitude', 'longitude'])


        da.attrs = dar.attrs
        da.values[pos[0], pos[1]] = np.nan  # mask sea

        if v == 'lw_out_PBLtop':

            da.values = olr_to_bt(da.values)
            da.values[da.values >= tthresh] = 0  # T threshold maskout
            da.values[np.isnan(da.values)] = 0 # set ocean nans to 0

            try:
                date = da.time.values[0]
            except IndexError:
                date = da.time.values

            labels, numL = label(da.values)

            u, inv = np.unique(labels, return_inverse=True)
            n = np.bincount(inv)

            goodinds = u[n >= 8]  # defines minimum MCS size e.g. 350 km2 = 39 pix at 3x3km res (258 pix at 4.4km is 5000km2) 52 pix is 1000km2 for cp4
            # 8 for 25km
            if not sum(goodinds) > 0:
                return

        if (v == 'lsRain') | (v == 'totRain'):
            da.values = da.values*3600  # rain to mm/h
            da.attrs['units'] = 'mm h-1'

        ds[v] = da

        print('Saved ', v)

    for gi in goodinds:
        if (gi == 0):  # index 0 is always background, ignore!
            continue

        inds = np.where(labels.flatten() == gi)
        #mask = np.where(labels!=gi)

        dbox = ds.copy(deep=True)

#         if np.nanmin(dbox['lw_out_PBLtop'].values[inds]) > -70:
#             continue

        pos = np.argmin((dbox['lw_out_PBLtop'].values.flatten()[inds]))

        tmin = np.nanmin((dbox['lw_out_PBLtop'].values.flatten()[inds]))
        print(tmin)

        # for v in dbox.data_vars:
        #     (dbox[v].values)[mask] = np.nan

        midlat = ((lats.flatten()[inds]))[pos] #latmin + (latmax-latmin)/2
        midlon = ((lons.flatten()[inds]))[pos] #lonmin + (lonmax-lonmin)/2

        # save only storms in defined storm box
        if (midlon<sbox[0]) | (midlon>sbox[1]) | (midlat<sbox[2]) | (midlat>sbox[3]):
                        continue

        point = dbox.sel(latitude=midlat, longitude=midlon, method='nearest')
        plat = point['latitude'].values
        plon = point['longitude'].values

        xpos = np.where(dbox['longitude'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(dbox['latitude'].values == plat)
        ypos = int(ypos[0])

        distx=103 # ca 200km i.e. 45 *4.4km
        disty = 103 # ca 80km i.e 18*4.4

        try:
            ds_box = dbox.isel(latitude=slice(ypos-disty,ypos+disty+1), longitude=slice(xpos-distx, xpos+distx+1))
        except IndexError:
            continue

        if (len(ds_box.latitude) != disty*2+1) | (len(ds_box.longitude) != distx*2+1):
            continue


        savefile = out_dir + os.sep + pd.Timestamp(date).strftime('%Y-%m-%d_%H:%M:%S') + '_' + str(gi) + '.nc'
        try:
            os.remove(savefile)
        except OSError:
            pass
        ds_box.attrs['tmin'] = tmin
        ds_box.to_netcdf(path=savefile, mode='w')
        print('Saved ' + savefile)


        print('Saved MCS no.'+str(gi)+ ' as netcdf.')



### Inputs CP25

data_path = cnst.network_data + 'data/CP4/CLOVER/CP4hist'  # CP4 data directory
ancils_path = cnst.network_data + 'data/CP4/ANCILS' # directory with seamask file inside
out_path = cnst.network_data + 'data/CP4/CLOVER'  # out directory to save MCS files
sbox = [-17, 14, 3.5, 21]  # W- E , S - N geographical storm box
datestring = '19990401'  # set this to date of file

years = np.array(np.arange(1998,2007), dtype=str)
months = np.array([ '03', '04', '05', '09', '10', '11'])
days = np.array(np.arange(1,32), dtype=str)

tthresh = -50 # chosen temperature threshold, e.g. -50, -60, -70

vars = OrderedDict()   # dictionary which contains info on pressure level and hour extraction for wanted variables
vars['lw_out_PBLtop'] = ([], 18)
vars['lsRain'] =  ([], 18)   # pressure levels, hour # totRain for CP25, lsRain for CP
vars['shear'] = ([650, 925], 12) # should use 925 later
vars['u_mid'] = ([650], 12)
vars['u_srfc'] = ([925], 12)
vars['q_pl'] = ([925], 12)  # 925, 650 available
vars['theta'] = ([650,925], 12)

datelist = []
# for y,m,d in itertools.product(years, months, days):
#     datelist.append(y+m+str(d).zfill(2))
# files = glob.glob(data_path+os.sep+'lw_out_PBLtop'+os.sep+'*.nc')

dummy = xr.open_mfdataset(data_path+os.sep+'lw_out_PBLtop'+os.sep+'*.nc', autoclose=True)

time = dummy['lw_out_PBLtop'][dummy['time.hour']==18].time

for d in time:
    file_save(data_path, out_path, ancils_path, vars, d, sbox, tthresh)
