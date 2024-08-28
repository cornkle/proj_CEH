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
import datetime
from collections import OrderedDict
from utils import constants as cnst, u_darrays as uda, u_interpolate as u_int
import matplotlib.pyplot as plt
import sys
from JASMIN import MetUM_variables as mu

def load_file(ffile,var):
    try:
        ds=xr.open_dataset(ffile)[var]
    except:
        #Deal with funny landseamask file
        ds=xr.open_dataset(ffile,decode_times=False)[var]
        ds=ds.rename({"rlat":"latitude","rlon":"longitude"})
    ds=ds.assign_coords({"longitude":(ds.longitude-360)})

    try:
        ds=ds.isel(pressure=slice(None,None,-1))
    except:
        pass
    return ds


def olr_to_bt(olr):
    #Application of Stefan-Boltzmann law
    sigma = 5.670373e-8
    tf = (olr/sigma)**0.25
    #Convert from bb to empirical BT (degC) - Yang and Slingo, 2001
    a = 1.228
    b = -1.106e-3
    Tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return Tb - 273.15


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


def filtering(dar, v, outv, pl):

    if (v == 'q_pl') | (v == 't_pl') | (v == 'lw_out_PBLtop'):
        dar.values = np.array(dar.values / 100).astype(float)

    if (v == 'SM'):
        dar = dar.sel(depth=0.05)
        dar = dar.where(dar< 500, other=np.nan)

    if 'pressure' in dar.coords:

        try:
            if v == 't_pl':
                dar.values[dar.values == -273.15] = np.nan
            else:
                dar.values[dar.values == 0] = np.nan  # potential missing value maskout
        except ValueError:
            ipdb.set_trace()
            print('Pressure value error!')
            return
        if (len(pl) > 1) & (outv == 'shear'):

            shear = dar.sel(pressure=650).values - dar.sel(pressure=925).values
            dar = dar.sum(dim='pressure').squeeze()
            dar.values = shear

        elif (len(pl) == 1) & (outv != 'shear'):
            dar = dar.sel(pressure=pl[0]).squeeze()

    return dar



def file_save(cp_dir, out_dir, vars, datestring, box, tthresh, pos):

    keys = vars.keys()

    if 'lw_out_PBLtop' not in keys:
        print('please provide ORL first in dictionary')
        return


    goodinds = 0

    #create empty dataset
    ds = xr.Dataset()
    # create empty

    #loop through every var
    for idx, outv in enumerate(keys):

        print('Variable ', outv)
        vv = outv

        h = (vars[outv])[1]
        pl = (vars[outv])[0]

        inds = (vars[outv][2])[0]
        weights = (vars[outv][2])[1]
        shape = (vars[outv][2])[2]

        v = (vars[outv])[3]
        grid = (vars[outv])[4]

        try:
            filepath = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*_'+datestring+'*.nc')[0]
        except IndexError:
            #ipdb.set_trace()
            print('No file found, return')
            return

        try:
            arr = xr.open_dataset(filepath)
        except OSError:
            print('Cannot open file, continue! ', filepath)
            return
        dar = arr[v].sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))

        if outv == 'lsRain_noon':
            dum = dar[dar['time.hour'] == h].squeeze()
            dar = dar[(dar['time.hour'] <= h) & (dar['time.hour'] >= h-2)].sum('time').squeeze() # early rain accumulation over 3h
            dar = dar.assign_coords(coords={'time' : dum.time})
            del dum
        else:
            dar = dar[dar['time.hour'] == h].squeeze()

        del arr


        dar = filtering(dar, v, outv, pl)

        ddate = dar.time.values.item()
        dt = datetime.datetime(ddate.year, ddate.month, ddate.day,
                               ddate.hour, ddate.minute, ddate.second)

        if outv not in ['lwout_noon', 'lsRain_noon', 'lw_out_PBLtop', 'lsRain']:
            mean_arr = []
            for dd in [7,6,5,4,3,2,1]:
                ndate = dt - pd.Timedelta(str(dd)+'days')
                if ndate.day >30:
                    dday = 30
                else:
                    dday = ndate.day
                ndatestring = str(ndate.year)+str(ndate.month).zfill(2)+str(dday).zfill(2)
                try:
                    nfile = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*_'+ndatestring+'*.nc')[0]
                except:
                    continue
                try:
                    meanarr = xr.open_dataset(nfile)
                    mmeanarr = meanarr[v].sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
                except OSError:
                    print('Couldnt find clim file, continue')
                    continue
                mmeanarr = mmeanarr[mmeanarr['time.hour'] == h].squeeze()

                mmeanarr = filtering(mmeanarr, v, outv, pl)
                mean_arr.append(mmeanarr)

                ndate = dt + pd.Timedelta(str(dd)+'days')
                if ndate.day >30:
                    dday = 30
                else:
                    dday = ndate.day
                ndatestring = str(ndate.year)+str(ndate.month).zfill(2)+str(dday).zfill(2)
                try:
                    nfile = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*_'+ndatestring+'*.nc')[0]
                except:
                    continue
                try:
                    meanarr = xr.open_dataset(nfile)
                    mmeanarr = meanarr[v].sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))
                except OSError:
                    print('Couldnt find clim file, continue')
                    continue
                mmeanarr = mmeanarr[mmeanarr['time.hour'] == h].squeeze()

                mmeanarr = filtering(mmeanarr, v, outv, pl)
                mean_arr.append(mmeanarr)
            print('Cases in clim mean:', len(mean_arr))

            if len(mean_arr) < 10:
                print('Mean smaller 10, return')
                return
            fullmean = xr.concat(mean_arr, dim='time').mean('time')
            dar = dar - fullmean

        # regrid to common grid (unstagger wind, bring to landsea mask grid)
        if grid == 'srfc':
            try:
                #regrid = griddata_lin(dar.values, dar.longitude, dar.latitude, ls_arr.rlon, ls_arr.rlat)

                regrid = u_int.interpolate_data(dar.values, inds, weights, shape)

            except ValueError:
                ipdb.set_trace()


            da = xr.DataArray(regrid,
                              coords={'time': dar.time, 'latitude': pl_dummy.latitude.values,
                                      'longitude': pl_dummy.longitude.values, },
                              dims=['latitude', 'longitude'])
            da.attrs = dar.attrs
        else:

            da = xr.DataArray(dar.values,
                              coords={'time': dar.time, 'latitude': pl_dummy.latitude.values,
                                      'longitude': pl_dummy.longitude.values, },
                               dims=['latitude', 'longitude'])
            da.attrs = dar.attrs


        da.values[pos[0], pos[1]] = np.nan  # mask sea

        if (v == 'lw_out_PBLtop') & (idx == 0):

            da.values[da.values >= tthresh] = 0  # T threshold maskout
            da.values[np.isnan(da.values)] = 0 # set ocean nans to 0

            date = pd.Timestamp(datestring)
            date = date.replace(hour=h)

            labels, numL = label(da.values)

            u, inv = np.unique(labels, return_inverse=True)
            n = np.bincount(inv)

            ## labels and filtering not needed if object location table is used for identification!
            goodinds = u[n >= 258]  # 258 == 5000km2 area threshold used here. defines minimum MCS size e.g. 350 km2 = 39 pix at 3x3km res (258 pix at 4.4km is 5000km2) 52 pix is 1000km2 for cp4
            if not sum(goodinds) > 0:
                print('No goodinds!')
                return

        if (v == 'lsRain') | (v == 'totRain'):
            da.values = da.values*3600  # rain to mm/h
            da.attrs['units'] = 'mm h-1'

        ds[outv] = da

        print('Saved ', outv, h)

    for gi in goodinds:
        if (gi == 0):  # index 0 is always background, ignore!
            continue
        mask = np.where(labels!=gi)
        dbox = ds.copy(deep=True)

        for vout in dbox.data_vars:
            if vout not in ['lsRain','lw_out_PBLtop']:
               continue
            (dbox[vout].values)[mask] = np.nan

        # Optional rain and cloud filter, keep off for now!
        # filt = dbox.where((dbox['lsRain_noon'] < 0.005) & (dbox['lwout_noon'] > -30))
        # for raw_var in ['lsRain', 'lw_out_PBLtop']:
        #     filt[raw_var] = dbox[raw_var]

        # location of minimum cloud top temperature, could come from table
        tmin = filt.where(dbox['lw_out_PBLtop'] == filt['lw_out_PBLtop'].min(), drop=True)

        # use this point location when minlon, minlat are saved in table as intermediate step - lon,lat will not correspond exactly to use "==" method
        # point = dbox.sel(latitude=tmin.latitude, longitude=tmin.longitude, method='nearest')
        plat = tmin['latitude'].values
        plon = tmin['longitude'].values
        try:
            if plat[0] < MINLAT:
                continue
        except:
            continue
        xpos = np.where(filt['longitude'].values == plon)
        ypos = np.where(filt['latitude'].values == plat)
        try:
            xpos = int(xpos[0])
            ypos = int(ypos[0])
        except TypeError:
            continue

        distx = 57  # 57 = 250 km at 4.4 res, 500km across
        disty = 57
        try:
            filt = filt.isel(latitude=slice(ypos - disty, ypos + disty + 1),
                             longitude=slice(xpos - distx, xpos + distx + 1))
        except IndexError:
            continue

        if (len(filt.latitude) != disty * 2 + 1) | (len(filt.longitude) != distx * 2 + 1):
            print(filt)
            continue
        filt = filt.assign_coords(
            {'longitude': np.arange(distx * -1, distx + 1), 'latitude': np.arange(distx * -1, distx + 1)})

        try:
            if np.nansum(filt['lsRain'])<=0.1:
                return
        except KeyError:
            if np.nansum(filt['totRain'])<=0.1:
                return
        filt.drop_vars(['lsRain_noon', 'lwout_noon'])
        filt.attrs = {'minlat' : plat, 'minlon' : plon}


        #if np.nansum(ds_box['lwout_noon'])<=0:
        #   ipdb.set_trace()
        #   print('lwout is lt 0', np.nansum(ds_box['lwout_noon']))
        #   return

        savefile = out_dir + os.sep + date.strftime('%Y-%m-%d_%H:%M:%S') + '_' + str(gi) + '_lonXlat_'+str(np.round(plon,1))+'_'+str(np.round(plat,1))+'.nc'
        try:
            os.remove(savefile)
        except OSError:
            pass

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        filt.to_netcdf(path=savefile, mode='w', encoding=encoding)
        print('Saved ' + savefile)
        print('Saved MCS no.'+str(gi)+ ' as netcdf.')


### Inputs:
##### Provide hour when you run script!!! as sys.argv[1]
fdir = str(sys.argv[2])
if fdir == 'hist':
    ftag = 'historical'
else:
    ftag = 'future'

main = '/home/users/cornkle/linked_CP4/'
main_lmcs = '/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/'
data_path = main + '/'+fdir
ancils_path = '/home/users/cornkle/impala/shared/CP4A/ncfiles/4km/ANCILS/'
out_path = main_lmcs + 'CP4_box_anom_JASMIN/CP4_'+ftag+'_5000km2_-50_box_anom_v3'
box = [-18, 25, 5, 25]  # W- E , S - N geographical coordinates box
MINLAT = 9

years = np.array(np.arange(1998,2007), dtype=str)
months = (['07', '08', '09'])
days = np.array(np.arange(1,32), dtype=str)

tthresh = -50 # chosen temperature threshold, e.g. -50, -60, -70
h= int(sys.argv[1]) # hour to extract

plglob = glob.glob(data_path + '/q_pl/*.nc')
pl_dummy = load_file(plglob[0], mu.create_CP4_filename('q_pl'))

srfcglob = glob.glob(data_path + '/lw_out_PBLtop/*.nc')
srfc_dummy =load_file(srfcglob[0], mu.create_CP4_filename('lw_out_PBLtop'))

## getting interpolation weights for later use
pl_dummy = pl_dummy.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
srfc_dummy = srfc_dummy.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
# load seamask
landsea_path = glob.glob(ancils_path + os.sep + 'landseamask*.nc')[0]
landsea = xr.open_dataset(landsea_path, decode_times=False)
ls = landsea['lsm']

ls = ls.assign_coords(rlon = ls.rlon.values - 360)
ls_arr = ls.sel(rlon=slice(box[0], box[1]), rlat=slice(box[2], box[3]))


pos = np.where(ls_arr[0, 0, :, :] == 0)
lons, lats = np.meshgrid(pl_dummy.longitude.values, pl_dummy.latitude.values)#np.meshgrid(ls_arr.rlon.values, ls_arr.rlat.values)
inds, weights, shape = u_int.interpolation_weights(srfc_dummy.longitude, srfc_dummy.latitude, pl_dummy.longitude, pl_dummy.latitude)

vars = OrderedDict()   # dictionary which contains info on pressure level and hour extraction for wanted variables
vars['lw_out_PBLtop'] = ([], h, (inds,weights,shape), 'lw_out_PBLtop', 'srfc')  ### Input in BRIGHTNESS TEMPERATURES!! (degC)
vars['lsRain'] =  ([], h, (inds,weights,shape), 'lsRain', 'srfc')   # pressure levels, hour
vars['shear'] = ([650, 925], 12, (0, 0, 0), 'u_pl', '') # (plinds, plweights, plshape) should use 925 later
vars['u_mid'] = ([650], 12, (0, 0, 0), 'u_pl', '')
vars['u_srfc'] = ([925], 12, (0, 0, 0), 'u_pl', '')
vars['v_mid'] = ([650], 12, (0, 0, 0), 'v_pl', '')
vars['v_srfc'] = ([925], 12, (0, 0, 0), 'v_pl', '')
vars['q_mid'] = ([650], 12, (0, 0, 0), 'q_pl', '')  # INPUT IN T * 100!!
vars['t_mid'] = ([650], 12, (0, 0, 0), 't_pl', '')   # INPUT IN T * 100!!
vars['t_srfc'] = ([925], 12, (0, 0, 0), 't_pl', '')
vars['q_srfc'] = ([925], 12, (0, 0, 0), 'q_pl', '')
vars['geoH_srfc'] = ([925], 12, (inds,weights,shape), 'geoH_pl', 'srfc')
vars['tcwv'] = ([], 12, (inds,weights,shape), 'tcwv', 'srfc')
vars['sh'] = ([], 12, (inds,weights,shape), 'sh', 'srfc')
vars['lh'] = ([], 12, (inds,weights,shape), 'lh', 'srfc')
vars['t2'] = ([], 12, (inds,weights,shape), 't2', 'srfc')
vars['q2'] = ([], 12, (inds,weights,shape), 'q2', 'srfc')
vars['lsRain_noon'] =  ([], 12, (inds,weights,shape), 'lsRain', 'srfc')
vars['lwout_noon'] =  ([], 12, (inds,weights,shape), 'lw_out_PBLtop', 'srfc')
vars['SM'] = ([], 12, (inds, weights, shape), 'SM', 'srfc')

datelist = []
for y,m,d in itertools.product(years, months, days):
    datelist.append(y+m+str(d).zfill(2))

for d in datelist:

    testfiles = glob.glob(out_path + os.sep + d[0:4] + '-' + d[4:6] + '-' + d[6:8] + '_' + str(h) + '*.nc')

    if len(testfiles) > 0:
        print(testfiles[0], ' already exists, continue!')
        continue

    print('Doing ', d)
    file_save(data_path, out_path, vars, d, box, tthresh, pos)

