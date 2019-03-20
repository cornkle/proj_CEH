# -*- coding: utf-8 -*-


import numpy as np
from wavelet import util
from eod import msg, mfg
import xarray as xr
import os
from scipy import ndimage
from utils import u_arrays as ua
import multiprocessing
from utils import u_grid, u_interpolate as u_int
import datetime as dt
from scipy.ndimage.measurements import label
from utils import constants as cnst
import pickle as pkl
import ipdb



def run(datastring):

    #msg_folder = cnst.network_data + 'data/OBS/meteosat_WA30'
    #ext_drive = '/media/ck/Seagate/DIR/'#
    #local_data = ext_drive + 'mymachine/'
    #network_data = ext_drive + 'cornkle/'

    #msg_folder = network_data + 'data/OBS/meteosat_WA30'

    for yy in range(2004,2005):   # (2004,2016)

        for mm in [9]:
            #
            # yy = 1999
            # mm = 9
            # datastring = 'mfg'

            pool = multiprocessing.Pool(processes=3)
            if datastring == 'mfg':
                msg_folder = cnst.network_data + '/data/OBS/MFG/'
                m = mfg.ReadMfg(msg_folder, y1=yy, y2=yy, months=[mm])

            files  = m.fpath

            gridll = pkl.load( open (cnst.network_data + 'data/OBS/saves/VERA_msg_latlon_18W12E_1N17N.p', 'rb'))

            mdic = m.read_data(files[0], llbox=[-25, 20, 2, 25])  #[-14, 2.5, 4, 11.5]

            # make salem grid
            grid = u_grid.make(gridll['lon'].values, gridll['lat'].values, 5000)
            inds, weights, shape = u_int.interpolation_weights_grid(mdic['lon'].values, mdic['lat'].values, grid)
            gridd = (inds,weights,shape, grid)

            files_str = []

            for f in files:
                if datastring == 'msg':
                    files_str.append(f[0:-4])
                else:
                    files_str.append(f)

            files_str = np.unique(files_str)

            passit = []
            for f in files_str:
                passit.append((gridd,m, f,datastring))

            res = pool.map(file_loop, passit)

            # res=[]
            # for l in passit:
            #
            #     res.append(file_loop(l))

            pool.close()

            res = [x for x in res if x is not None]

            ds = xr.concat(res, 'time')
            path =  cnst.network_data + 'MCSfiles/VERA_blobs/' # '/prj/vera/cores/' # cnst.network_data + 'MCSfiles/VERA_blobs/'
            savefile = path + 'cores_powerTest_MFG_-40_700km2_-50points_dominant_'+str(yy) + '_'+str(mm).zfill(2)+'.nc'#'blobMap_-40-700km2_-50-points_dominant_'+str(yy) + '_'+str(mm).zfill(2)+'.nc'

            try:
                os.remove(savefile)
            except OSError:
                pass
            #da.name = 'blob'
            #enc = {'blob': {'complevel': 5, 'zlib': True}}

            comp = dict(zlib=True, complevel=5)
            enc = {var: comp for var in ds.data_vars}

            ds.to_netcdf(path=savefile, mode='w', encoding=enc, format='NETCDF4')
            print('Saved ' + savefile)


def filter_img(inarr):

        outt = inarr.copy()
        print('outmin', np.nanmin(outt), np.nanmax(outt))

        t_thresh_size = -40
        t_thresh_cut = -50

        outt[outt >= t_thresh_size] = 0
        outt[np.isnan(outt)] = 0

        labels, numL = label(outt)

        u, inv = np.unique(labels, return_inverse=True)
        n = np.bincount(inv)

        pix_nb = 28

        badinds = u[(
                    n < pix_nb)]  # all blobs with more than 1000 pixels = 25,000km2 (meteosat regridded 5km), 200pix = 5000km2, 8pix = 200km2
        # scale 30km, radius 15km ca. 700km2 circular area equals 28 pix

        for bi in badinds:
            inds = np.where(labels == bi)
            outt[inds] = 0

        outt[outt >= t_thresh_cut] = 150

        grad = np.gradient(outt)
        outt[outt == 150] = np.nan

        nogood = np.isnan(outt)  # filters edge maxima later, no maxima in -40 edge area by definition!

        # tdiff = np.nanmax(outt) - np.nanmin(outt)  # define background temperature for image
        # if tdiff > 28:  # temp difference of 28 degrees
        #     xmin = 15
        # else:
        #     xmin = 10

        xmin = 10
        outt[nogood] = t_thresh_cut - xmin
        nok = np.where(abs(grad[0]) > 80)
        d = 2
        i = nok[0]
        j = nok[1]
        # edge smoothing for wavelet application
        for ii, jj in zip(i, j):
            kern = outt[ii - d:ii + d + 1, jj - d:jj + d + 1]
            outt[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

        return outt, nogood, t_thresh_size, t_thresh_cut, pix_nb


def find_scales_dominant(wav, outt, core_min=None, no_good=None, dataset=None):
    outt[no_good] = np.nan

    arr = np.array(wav['scales'], dtype=str)

    wll = wav['t']

    power_img = np.sum(wll, axis=0)
    power_img[no_good] = 0

    if dataset == 'mfg':
        smaller = -20
        thresh_p = np.sum((wav['scales'] + smaller) ** .5)
        power_img[(power_img < np.percentile(power_img[power_img > 0.001], 25)) | (power_img < (thresh_p))] = 0
    if dataset == 'msg':
        smaller = -12
        thresh_p = np.sum((wav['scales'] + smaller) ** .5)
        power_img[(power_img < np.percentile(power_img[power_img > 0.001], 25)) | (power_img < (thresh_p))] = 0

    labels, numL = label(power_img)
    u, inv = np.unique(labels, return_inverse=True)

    for inds in u:
        if inds == 0:
            continue

        arr = power_img.copy()
        arr[np.where(labels != inds)] = 0

        power_img.flat[np.argmax(arr)] = -999

    return power_img



def file_loop(passit):

    gridd = passit[0]
    inds_inter = gridd[0]
    weights_inter = gridd[1]
    shape_inter = gridd[2]
    grid_inter = gridd[3]

    m = passit[1]
    file = passit[2]
    tag = passit[3]


    if tag == 'mfg':
            print(' month')
    else:
        strr = file.split(os.sep)[-1]

        if ((strr[-2::]) != '00') & ((strr[-2::]) != '30'):
            print('Skip minute')
            return

    # if not ((np.int(strr[8:10]) >= 20)): #& (np.int(strr[8:10]) <= 19) ): #((np.int(strr[8:10]) > 3)): #not ((np.int(strr[8:10]) >= 16) & (np.int(strr[8:10]) <= 19) ): #& (np.int(strr[8:10]) < 18): #(np.int(strr[4:6]) != 6) & #(np.int(strr[8:10]) != 3) , (np.int(strr[8:10]) > 3)
    #     print('Skip hour')
    #     return

    lon, lat = grid_inter.ll_coordinates

    if tag != 'mfg':
        file = file +'.gra'

    print('Doing file: ' + file)
    try:
        mdic = m.read_data(file, llbox=[-25, 20, 2, 25])
    except FileNotFoundError:
        print('File not found')
        return

    if not mdic:
        print('File missing')
        return

    bloblist = []
    tirlist = []

    for id, timeslice in enumerate(mdic['t']):

        #test = timeslice.copy(deep=True)
        try:
            outt = u_int.interpolate_data(timeslice.values, inds_inter, weights_inter, shape_inter)
            print('Interpolated', id)
        except ValueError:
            print('Interpolation value error!!')
            ipdb.set_trace()
            continue


        outt, nogood, t_thresh_size, t_thresh_cut, pix_nb = filter_img(outt.values)

        wav = util.waveletT(outt, dataset='METEOSAT5K_vera')
        core_min=-50
        outt[nogood] = np.nan

        power_msg = find_scales_dominant(wav, outt, no_good=nogood, core_min=-50, dataset='mfg')

        hour = timeslice['time.hour']
        minute = timeslice['time.minute']
        day = timeslice['time.day']
        month = timeslice['time.month']
        year = timeslice['time.year']

        date = dt.datetime(year, month, day, hour, minute)

        isnan = np.isnan(outt)
        outt[isnan]=0
        new_savet = (outt*100).astype(np.int16)

        bloblist.append(xr.DataArray(power_msg.astype(np.int16), coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon'])) #[np.newaxis, :])
        tirlist.append(xr.DataArray(new_savet, coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon']))

    ds = xr.Dataset()
    ds['blobs'] = xr.concat(bloblist, 'time')
    ds['tir'] = xr.concat(tirlist, 'time')


    ds.attrs['radii']=(np.ceil(wav['scales'] / 2. / 5.)).astype(np.uint8)
    ds.attrs['scales_rounded'] = np.round(wav['scales']).astype(np.uint8)
    ds.attrs['scales_original'] = wav['scales']
    ds.attrs['cutout_T'] = t_thresh_size
    ds.attrs['core_minT'] = core_min
    ds.attrs['cutout_minPixelNb'] = pix_nb

    ds = ds.sel(lat=slice(2,17), lon=slice(-18,13))     #[-14, 2.5, 4, 11.5] cutout to remove dodgy boundaries

    print('Did ', file)

    return (ds)
