import numpy as np
from scipy.ndimage.measurements import label
import xarray as xr
import os
import ipdb
import glob
from scipy.interpolate import griddata
import pandas as pd
import pdb
import itertools


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

def cut_kernel(array, xpos, ypos, dist_from_point):
    """
     This function cuts out a kernel from an existing array and allows the kernel to exceed the edges of the input
     array. The cut-out area is shifted accordingly within the kernel window with NaNs filled in
    :param array: 2darray
    :param xpos: middle x point of kernel
    :param ypos: middle y point of kernel
    :param dist_from_point: distance to kernel edge to each side
    :return: 2d array of the chosen kernel size.
    """

    if array.ndim != 2:
        raise IndexError('Cut kernel only allows 2D arrays.')

    kernel = np.zeros((dist_from_point*2+1, dist_from_point*2+1)) * np.nan

    if xpos - dist_from_point >= 0:
        xmin = 0
        xmindist = dist_from_point
    else:
        xmin = (xpos - dist_from_point) * -1
        xmindist = dist_from_point + (xpos - dist_from_point)

    if ypos - dist_from_point >= 0:
        ymin = 0
        ymindist = dist_from_point
    else:
        ymin = (ypos - dist_from_point) * -1
        ymindist = dist_from_point + (ypos - dist_from_point)

    if xpos + dist_from_point < array.shape[1]:
        xmax = kernel.shape[1]
        xmaxdist = dist_from_point + 1
    else:
        xmax = dist_from_point - (xpos - array.shape[1])
        xmaxdist = dist_from_point - (xpos + dist_from_point - array.shape[1])

    if ypos + dist_from_point < array.shape[0]:
        ymax = kernel.shape[0]
        ymaxdist = dist_from_point + 1
    else:
        ymax = dist_from_point - (ypos - array.shape[0])
        ymaxdist = dist_from_point - (ypos + dist_from_point - array.shape[0])

    cutk = array[ypos - ymindist: ypos + ymaxdist, xpos - xmindist: xpos + xmaxdist]


    kernel[ymin: ymax, xmin:xmax] = cutk

    return kernel


def cut_box(xpos, ypos, arr, dist=None):
    """

    :param xpos: x coordinate in domain for kernel centre point
    :param ypos: y coordinate in domain for kernel centre point
    :param arr: numpy array (2d)
    :param dist: distance from kernel centre point to kernel edge (total width = 2*dist+1)
    :return: the kernel of dimensions (2*dist+1, 2*dist+1)
    """

    if dist == None:
        'Distance missing. Please provide distance from kernel centre to edge (number of pixels).'
        return

    kernel = cut_kernel(arr,xpos, ypos,dist)
    if kernel.shape != (dist * 2 + 1, dist * 2 + 1):
        print("Please check kernel dimensions, there is something wrong")
        pdb.set_trace()

    return kernel


def file_save(cp_dir, out_dir, ancils_dir, vars, datestring, box, tthresh):

    keys = vars.keys()

    if 'lw_out_PBLtop' not in keys:
        print('please provide ORL first in dictionary')
        return

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
        if (v == 'shear') | (v == 'u_mid') | (v == 'u_srfc'):
            derived = v
            v = 'u_pl'

        try:
            filepath = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*'+datestring+'*.nc')[0]
        except IndexError:
            print('No file found, return')
            return
        arr = xr.open_dataset(filepath)

        dar = arr[v].sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))

        dar = dar[dar['time.hour']==h].squeeze()

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

            goodinds = u[n >= 258]  # defines minimum MCS size e.g. 258 pix at 4.4km is 5000km2, 52 pix is 1000km2
            if not sum(goodinds) > 0:
                return

        if v == 'lsRain':
            da.values = da.values*3600  # rain to mm/h
            da.attrs['units'] = 'mm h-1'

        ds[v] = da

        print('Written ', v)

    for gi in goodinds:
        if (gi == 0):  # index 0 is always background, ignore!
            continue
        inds = np.where(labels == gi)
        #mask = np.where(labels!=gi)

        dbox = ds.copy(deep=True)

        if np.nanmin(dbox['lw_out_PBLtop'].values[inds]) > -70:
            continue

        # for v in dbox.data_vars:
        #     (dbox[v].values)[mask] = np.nan

        # cut a box for every single blob from msg - get min max lat lon of the blob
        latmax, latmin = np.nanmax(lats[inds]), np.nanmin(lats[inds])
        lonmax, lonmin = np.nanmax(lons[inds]), np.nanmin(lons[inds])

        midlat = latmin + (latmax-latmin)/2
        midlon = lonmin + (lonmax-lonmin)/2

        point = dbox.sel(latitude=midlat, longitude=midlon, method='nearest')
        plat = point['latitude'].values
        plon = point['longitude'].values

        xpos = np.where(dbox['longitude'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(dbox['latitude'].values == plat)
        ypos = int(ypos[0])

        distx=45 # ca 200km i.e. 45 *4.4km
        disty = 18 # ca 80km i.e 18*4.4

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

        ds_box.to_netcdf(path=savefile, mode='w')
        print('Saved ' + savefile)



        print('Saved MCS no.'+str(gi)+ ' as netcdf.')



### Inputs:

data_path = '/users/global/cornkle/shared/data/CP4/CLOVER/CP4hist'  # CP4 data directory
ancils_path = '/users/global/cornkle/shared/data/CP4/ANCILS' # directory with seamask file inside
out_path = '/users/global/cornkle/shared/data/CP4/CLOVER/test'  # out directory to save MCS files
box = [-13, 13, 4.5, 10]  # W- E , S - N geographical coordinates box
datestring = '19990401'  # set this to date of file

years = np.array(np.arange(1998,2007), dtype=str)
months = np.array(['10'])
days = np.array(np.arange(1,32), dtype=str)

tthresh = -50 # chosen temperature threshold, e.g. -50, -60, -70

vars = {}   # dictionary which contains info on pressure level and hour extraction for wanted variables
vars['lw_out_PBLtop'] = ([], 18)
vars['lsRain'] =  ([], 18)   # pressure levels, hour
vars['shear'] = ([650, 925], 12) # should use 925 later
vars['u_mid'] = ([650], 12)
vars['u_srfc'] = ([925], 12)
vars['q_pl'] = ([925], 12)  # 925, 650 available
datelist = []
for y,m,d in itertools.product(years, months, days):
    datelist.append(y+m+str(d).zfill(2))

for d in datelist:
    file_save(data_path, out_path, ancils_path, vars, d, box, tthresh)
