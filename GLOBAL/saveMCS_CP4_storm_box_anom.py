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

        ddate = dar.time.values.item()
        dt = datetime.datetime(ddate.year, ddate.month, ddate.day,
                               ddate.hour, ddate.minute, ddate.second)

        mean_arr = []
        for dd in [7,6,5,4,3,2,1]:
            ndate = dt - pd.Timedelta(str(dd)+'days')
            ndatestring = str(ndate.year)+str(ndate.month)+str(ndate.day).zfill(2)
            nfile = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*_'+ndatestring+'*.nc')[0]
            try:
                meanarr = xr.open_dataset(nfile)
                mmeanarr = meanarr[v].sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
            except OSError:
                print('Couldnt find clim file, continue')
                continue
            mean_arr.append(mmeanarr)

            ndate = dt + pd.Timedelta(str(dd)+'days')
            ndatestring = str(ndate.year)+str(ndate.month)+str(ndate.day).zfill(2)
            nfile = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*_'+ndatestring+'*.nc')[0]
            try:
                meanarr = xr.open_dataset(nfile)
                mmeanarr = meanarr[v].sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))
            except OSError:
                print('Couldnt find clim file, continue')
                continue
            mean_arr.append(mmeanarr)
        print('Cases in clim mean:', len(mean_arr))
        fullmean = xr.concat(mean_arr, dim='time').mean('time')

        dar = dar - fullmean

        if outv == 'lsRain_noon':
            dum = dar[dar['time.hour'] == h].squeeze()
            dar = dar[(dar['time.hour'] <= h) & (dar['time.hour'] >= h-2)].sum('time').squeeze() # early rain accumulation over 3h
            dar = dar.assign_coords(coords={'time' : dum.time})
            del dum
        else:
            dar = dar[dar['time.hour'] == h].squeeze()

        del arr

        if (v == 'q_pl') | (v=='t_pl') | (v=='lw_out_PBLtop'):
            dar.values = np.array(dar.values/100).astype(float)


        #if int(dar['time.hour'])!=h:
        #    ipdb.set_trace()
        #    print('Wrong hour')
        #    return

        if 'pressure' in dar.coords:
            try:
                if v == 't_pl':
                    dar.values[dar.values == -273.15] = np.nan
                else:
                    dar.values[dar.values==0] = np.nan # potential missing value maskout
            except ValueError:
                print('Pressure value error!')
                return
            if (len(pl) > 1) & (outv == 'shear'):

                shear = dar.sel(pressure=650).values - dar.sel(pressure=925).values
                dar = dar.sum(dim='pressure').squeeze()
                dar.values = shear

            elif (len(pl) == 1)  & (outv != 'shear'):
                dar = dar.sel(pressure=pl[0]).squeeze()


        # regrid to common grid (unstagger wind, bring to landsea mask grid)
        if grid == 'srfc':
            try:
                #regrid = griddata_lin(dar.values, dar.longitude, dar.latitude, ls_arr.rlon, ls_arr.rlat)

                regrid = u_int.interpolate_data(dar.values, inds, weights, shape)

            except ValueError:
                ipdb.set_trace()

            # if v == 'sh':  # potential calculation of SH climatology around SH date. Not finished
            #
            #     ddate = dar.time.values.item()
            #     dt = datetime.datetime(ddate.year, ddate.month, ddate.day,
            #                            ddate.hour, ddate.minute, ddate.second)
            #
            #     window1 = dt - pd.Timedelta('7days')
            #     window2 = dt + pd.Timedelta('7days')
            #
            #     doi = pd.date_range(window1, window2)
            #     c1ds = []
            #     for doi_id in doi:
            #
            #         cdate = str(doi_id.year) + str(doi_id.month).zfill(2) + str(doi_id.day).zfill(2) + '0030'
            #
            #         try:
            #             c1 = xr.open_dataarray(glob.glob(data_path + '/sh/*_'+cdate+'-*.nc')[0])
            #             c1ds.append(c1)
            #         except:
            #             pass
            #     if len(c1ds) < 10:
            #         print('Not enough dates to calc sh clim')
            #         return
            #
            #     ipdb.set_trace()
            #     shmean = xr.concat(c1ds, dim='time').mean('time')


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


            # try:
            #     date = da.time.values[0]
            # except IndexError:
            #     date = da.time.values

            date = pd.Timestamp(datestring)
            date = date.replace(hour=h)

            labels, numL = label(da.values)

            u, inv = np.unique(labels, return_inverse=True)
            n = np.bincount(inv)

            goodinds = u[n >= 517]  # 517 == 10000km2. defines minimum MCS size e.g. 350 km2 = 39 pix at 3x3km res (258 pix at 4.4km is 5000km2) 52 pix is 1000km2 for cp4
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
        inds = np.where(labels == gi)
        # cut a box for every single blob from msg - get min max lat lon of the blob
        # latmax, latmin = np.nanmax(lats[inds]), np.nanmin(lats[inds])
        # lonmax, lonmin = np.nanmax(lons[inds]), np.nanmin(lons[inds])
        mask = np.where(labels!=gi)

        dbox = ds.copy(deep=True)

        for vout in dbox.data_vars:
            if vout not in ['lsRain','lw_out_PBLtop']:
               continue
            (dbox[vout].values)[mask] = np.nan

        filt = dbox.where((dbox['lsRain_noon'] < 0.005) & (dbox['lwout_noon'] > -30))

        tmin = filt.where(dbox['lw_out_PBLtop'] == filt['lw_out_PBLtop'].min(), drop=True)

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
        ds.attrs = {'minlat' : plat, 'minlon' : plon}
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
if fdir == 'CP4hist':
    ftag = 'historical'
else:
    ftag = 'future'

main = cnst.other_drive + 'CP4'
main_lmcs = cnst.lmcs_drive + 'CP_models'
data_path = main + '/CP4_WestAfrica/'+fdir
ancils_path = main + '/CP4_WestAfrica/ANCILS'
out_path = main_lmcs + '/MCS_files/MODELS/CP4_box/CP4_allHours_'+ftag+'_5000km2_-50_WAf_box'
box = [-18, 25, 5, 25]  # W- E , S - N geographical coordinates box
MINLAT = 9

years = np.array(np.arange(1998,2007), dtype=str)
months = ([ '06', '07', '08', '09'])
days = np.array(np.arange(1,32), dtype=str)

tthresh = -50 # chosen temperature threshold, e.g. -50, -60, -70
h= int(sys.argv[1]) # hour to extract

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

ls = ls.assign_coords(rlon = ls.rlon.values - 360)
ls_arr = ls.sel(rlon=slice(box[0], box[1]), rlat=slice(box[2], box[3]))
###########
topo_path = glob.glob(ancils_path + os.sep + 'orog_combined*.nc')[0]
topo = xr.open_dataset(topo_path, decode_times=False)
top = topo['ht']

top = top.assign_coords(rlon = top.rlon.values - 360)
top_arr = top.sel(rlon=slice(box[0], box[1]), rlat=slice(box[2], box[3]))
##########

pos = np.where(ls_arr[0, 0, :, :] == 0)
lons, lats = np.meshgrid(pl_dummy.longitude.values, pl_dummy.latitude.values)#np.meshgrid(ls_arr.rlon.values, ls_arr.rlat.values)

#plinds, plweights, plshape = u_int.interpolation_weights(pl_dummy.longitude, pl_dummy.latitude, ls_arr.rlon, ls_arr.rlat)
inds, weights, shape = u_int.interpolation_weights(srfc_dummy.longitude, srfc_dummy.latitude, pl_dummy.longitude, pl_dummy.latitude)
#regrid = griddata_lin(dar.values, dar.longitude, dar.latitude, ls_arr.rlon, ls_arr.rlat)


vars = OrderedDict()   # dictionary which contains info on pressure level and hour extraction for wanted variables
vars['lw_out_PBLtop'] = ([], h, (inds,weights,shape), 'lw_out_PBLtop', 'srfc')  ### Input in BRIGHTNESS TEMPERATURES!! (degC)
vars['lsRain'] =  ([], h, (inds,weights,shape), 'lsRain', 'srfc')   # pressure levels, hour
vars['shear'] = ([650, 925], 12, (0, 0, 0), 'u_pl', '') # (plinds, plweights, plshape) should use 925 later
vars['u_mid'] = ([650], 12, (0, 0, 0), 'u_pl', '')
vars['u_srfc'] = ([925], 12, (0, 0, 0), 'u_pl', '')
vars['q_mid'] = ([650], 12, (0, 0, 0), 'q_pl', '')  # INPUT IN T * 100!!
vars['t_mid'] = ([650], 12, (0, 0, 0), 't_pl', '')   # INPUT IN T * 100!!
vars['t_low'] = ([850], 12, (0, 0, 0), 't_pl', '')
vars['t_srfc'] = ([925], 12, (0, 0, 0), 't_pl', '')
vars['q_srfc'] = ([925], 12, (0, 0, 0), 'q_pl', '')
vars['tcwv'] = ([], 12, (inds,weights,shape), 'tcwv', 'srfc')
vars['sh'] = ([], 12, (inds,weights,shape), 'sh', 'srfc')
vars['lh'] = ([], 12, (inds,weights,shape), 'lh', 'srfc')
vars['t2'] = ([], 12, (inds,weights,shape), 't2', 'srfc')
vars['q2'] = ([], 12, (inds,weights,shape), 'q2', 'srfc')
vars['lsRain_noon'] =  ([], 12, (inds,weights,shape), 'lsRain', 'srfc')
vars['lwout_noon'] =  ([], 12, (inds,weights,shape), 'lw_out_PBLtop', 'srfc')

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

