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
from utils import constants as cnst, u_met, u_darrays as uda, u_interpolate as u_int
import matplotlib.pyplot as plt

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

def cut_kernel_3d(array, xpos, ypos, dist_from_point):
    """
     This function cuts out a kernel from an existing array and allows the kernel to exceed the edges of the input
     array. The cut-out area is shifted accordingly within the kernel window with NaNs filled in
    :param array: 2darray
    :param xpos: middle x point of kernel
    :param ypos: middle y point of kernel
    :param dist_from_point: distance to kernel edge to each side
    :return: 2d array of the chosen kernel size.
    """

    if array.ndim != 3:
        raise IndexError('Cut kernel3d only allows 3D arrays.')

    kernel = np.zeros((array.shape[0], dist_from_point*2+1, dist_from_point*2+1)) * np.nan

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

    if xpos + dist_from_point < array.shape[2]:
        xmax = kernel.shape[2]
        xmaxdist = dist_from_point + 1
    else:
        xmax = dist_from_point - (xpos - array.shape[2])
        xmaxdist = dist_from_point - (xpos + dist_from_point - array.shape[2])

    if ypos + dist_from_point < array.shape[1]:
        ymax = kernel.shape[1]
        ymaxdist = dist_from_point + 1
    else:
        ymax = dist_from_point - (ypos - array.shape[1])
        ymaxdist = dist_from_point - (ypos + dist_from_point - array.shape[1])

    cutk = array[:, ypos - ymindist: ypos + ymaxdist, xpos - xmindist: xpos + xmaxdist]


    kernel[:, ymin: ymax, xmin:xmax] = cutk

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
    if arr.ndim == 3:
        kernel = cut_kernel_3d(arr, xpos, ypos, dist)
        if kernel.shape != (kernel.size[0], dist * 2 + 1, dist * 2 + 1):
            print("Please check kernel dimensions, there is something wrong")
            ipdb.set_trace()
    else:
        kernel = cut_kernel(arr,xpos, ypos,dist)
        if kernel.shape != (dist * 2 + 1, dist * 2 + 1):
            print("Please check kernel dimensions, there is something wrong")
            ipdb.set_trace()



    return kernel


def file_save(cp_dir, out_dir, ancils_dir, vars, datestring, box, tthresh, pos, lons, lats):

    keys = vars.keys()

    if 'lw_out_PBLtop' not in keys:
        print('please provide ORL first in dictionary')
        return


    goodinds = 0

    #create empty dataset
    ds = xr.Dataset()
    # create empty

    #loop through every var
    for v in keys:

        print('Variable ', v)
        vv = v

        h = (vars[v])[1]
        pl = (vars[v])[0]

        inds = (vars[v][2])[0]
        weights = (vars[v][2])[1]
        shape = (vars[v][2])[2]



        derived = False
        if (v == 'shear') | (v == 'u_mid') | (v == 'u_srfc'):
            derived = v
            v = 'u_pl'

        if (v == 'q_srfc') | (v == 'q_mid') :
            derived = v
            v = 'q_pl'

        if (v == 't_srfc') | (v == 't_mid') |  (v == 'theta'):
            derived = v
            v = 't_pl'


        try:
            filepath = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*'+datestring+'*.nc')[0]
        except IndexError:
            print('No file found, return')
            return

        try:
            arr = xr.open_dataset(filepath)
        except OSError:
            print('Cannot open file, continue! ', filepath)
            return

        dar = arr[v].sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))

        dar = dar[dar['time.hour'] == h].squeeze()


        del arr

        if (v == 'q_pl') | (v=='t_pl') | (v=='lw_out_PBLtop'):

            dar.values = np.array(dar.values/100).astype(float)


        if int(dar['time.hour'])!=h:
            ipdb.set_trace()
            print('Wrong hour')
            return

        if 'pressure' in dar.coords:
            try:
                dar.values[dar.values==0] = np.nan # potential missing value maskout
            except ValueError:
                print('Pressure value error!')
                return
            if (len(pl) > 1) & (vv == 'shear'):

                shear = dar.sel(pressure=650).values - dar.sel(pressure=925).values
                dar = dar.sum(dim='pressure').squeeze()
                dar.values = shear
                if derived:
                    v = derived


            elif (len(pl) == 1)  & (vv != 'shear'):
                dar = dar.sel(pressure=pl[0]).squeeze()


            if derived:
                v = derived

        # regrid to common grid (unstagger wind, bring to landsea mask grid)
        if v in ['lw_out_PBLtop', 'tcwv', 'lsRain']:
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

            #ipdb.set_trace()

        da.values[pos[0], pos[1]] = np.nan  # mask sea

        if v == 'lw_out_PBLtop':


            da.values[da.values >= tthresh] = 0  # T threshold maskout
            da.values[np.isnan(da.values)] = 0 # set ocean nans to 0


            # try:
            #     date = da.time.values[0]
            # except IndexError:
            #     date = da.time.values

            date = pd.Timestamp(datestring)
            date = date.replace(hour=h)

            labels, numL = label(da.values)

            u, inv = np.unique(labels, return_inverse=True)
            n = np.bincount(inv)

            goodinds = u[n >= 258]  # defines minimum MCS size e.g. 350 km2 = 39 pix at 3x3km res (258 pix at 4.4km is 5000km2) 52 pix is 1000km2 for cp4
            if not sum(goodinds) > 0:
                print('No goodinds!')
                return

        if (v == 'lsRain') | (v == 'totRain'):
            da.values = da.values*3600  # rain to mm/h
            da.attrs['units'] = 'mm h-1'

        ds[v] = da

        print('Saved ', v, h)


    for gi in goodinds:
        if (gi == 0):  # index 0 is always background, ignore!
            continue
        inds = np.where(labels == gi)
        # cut a box for every single blob from msg - get min max lat lon of the blob
        latmax, latmin = np.nanmax(lats[inds]), np.nanmin(lats[inds])
        lonmax, lonmin = np.nanmax(lons[inds]), np.nanmin(lons[inds])
        mask = np.where(labels!=gi)

        dbox = ds.copy(deep=True)

        tgrad = dbox['t_srfc'].sel(longitude=slice(lonmin,lonmax)).mean('longitude')

        tmin = np.nanargmin(tgrad.values)
        tmax = np.nanargmax(tgrad.values)
        tgrad = tgrad.isel(latitude=slice(tmin,tmax))
        #ipdb.set_trace()
        lingress = uda.linear_trend_lingress(tgrad)

        dbox.attrs['Tgrad'] = lingress['slope'].values

        tgrad2 = dbox['t_srfc'].sel(longitude=slice(lonmin, lonmax), latitude=slice(10,20)).mean(['longitude', 'latitude'])- \
        dbox['t_srfc'].sel(longitude=slice(lonmin, lonmax), latitude=slice(5,7)).mean(['longitude', 'latitude'])
        dbox.attrs['Tgradbox'] = tgrad2.values


        for v in dbox.data_vars:
            (dbox[v].values)[mask] = np.nan

        ds_box = dbox.sel(latitude=slice(latmin,latmax), longitude=slice(lonmin, lonmax))
        try:
            if np.nansum(ds_box['lsRain'])==0:
                return
        except KeyError:
            if np.nansum(ds_box['totRain'])==0:
                return

        savefile = out_dir + os.sep + date.strftime('%Y-%m-%d_%H:%M:%S') + '_' + str(gi) + '.nc'
        try:
            os.remove(savefile)
        except OSError:
            pass

        ds_box.to_netcdf(path=savefile, mode='w')
        print('Saved ' + savefile)


        print('Saved MCS no.'+str(gi)+ ' as netcdf.')


### Inputs:

data_path = '/media/ck/Elements/Africa/WestAfrica/CP4/CP4fut' #cnst.network_data + 'data/CP4/CLOVER/CP4hist'  # CP4 data directory
ancils_path = cnst.network_data + 'data/CP4/ANCILS' # directory with seamatotRainsk file inside
out_path = cnst.network_data + 'data/CP4/CLOVER/CP4_16-19UTC_future_5000km2_-40C_TCWV'  # out directory to save MCS files
box = [-12, 15, 5, 25]  # W- E , S - N geographical coordinates box
#datestring = '19990301'  # set this to date of file

years = np.array(np.arange(1998,2007), dtype=str)
months = np.array([ '03', '04', '05','06','07','08','09', '10', '11'])#([ '03', '04', '05', '06', '09', '10', '11'])
days = np.array(np.arange(1,32), dtype=str)

tthresh = -40 # chosen temperature threshold, e.g. -50, -60, -70
h=18

plglob = glob.glob(data_path + '/q_pl/*.nc')
pl_dummy = xr.open_dataset(plglob[0])

srfcglob = glob.glob(data_path + '/lw_out_PBLtop/*.nc')
srfc_dummy = xr.open_dataset(srfcglob[0])

pl_dummy = pl_dummy.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
srfc_dummy = srfc_dummy.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
# load seamask
landsea_path = glob.glob(ancils_path + os.sep + 'landseamask*.nc')[0]
landsea = xr.open_dataset(landsea_path, decode_times=False)
ls = landsea['lsm']

ls.rlon.values = ls.rlon.values - 360
ls_arr = ls.sel(rlon=slice(box[0], box[1]), rlat=slice(box[2], box[3]))

pos = np.where(ls_arr[0, 0, :, :] == 0)
lons, lats = np.meshgrid(pl_dummy.longitude.values, pl_dummy.latitude.values)#np.meshgrid(ls_arr.rlon.values, ls_arr.rlat.values)

#plinds, plweights, plshape = u_int.interpolation_weights(pl_dummy.longitude, pl_dummy.latitude, ls_arr.rlon, ls_arr.rlat)
inds, weights, shape = u_int.interpolation_weights(srfc_dummy.longitude, srfc_dummy.latitude, pl_dummy.longitude, pl_dummy.latitude)
#regrid = griddata_lin(dar.values, dar.longitude, dar.latitude, ls_arr.rlon, ls_arr.rlat)

#ipdb.set_trace()

vars = OrderedDict()   # dictionary which contains info on pressure level and hour extraction for wanted variables
vars['lw_out_PBLtop'] = ([], h, (inds,weights,shape))  ### Input in BRIGHTNESS TEMPERATURES!! (degC)
vars['lsRain'] =  ([], h, (inds,weights,shape))   # pressure levels, hour
vars['shear'] = ([650, 925], 12, (0, 0, 0)) # (plinds, plweights, plshape) should use 925 later
vars['u_mid'] = ([650], 12, (0, 0, 0))
vars['u_srfc'] = ([925], 12, (0, 0, 0))
vars['q_mid'] = ([650], 12, (0, 0, 0))  # INPUT IN T * 100!!
vars['t_mid'] = ([650], 12, (0, 0, 0))   # INPUT IN T * 100!!
vars['t_srfc'] = ([850], 12, (0, 0, 0))
vars['q_srfc'] = ([925], 12, (0, 0, 0))
vars['tcwv'] = ([], 12, (inds,weights,shape))


datelist = []
for y,m,d in itertools.product(years, months, days):
    datelist.append(y+m+str(d).zfill(2))

for d in datelist:

    testfiles = glob.glob(out_path + os.sep + d[0:4] + '-' + d[4:6] + '-' + d[6:8] + '_' + str(h) + '*.nc')

    if len(testfiles) > 0:
        print(testfiles[0], ' already exists, continue!')
        continue


    file_save(data_path, out_path, ancils_path, vars, d, box, tthresh, pos, lons, lats)

# for d in datelist[0:10]:
#
#     if (d['time.year']<1998) | (d['time.month']<3) | (d['time.month']>11):
#         continue
#     file_save(data_path, out_path, ancils_path, vars, d, box, tthresh)
