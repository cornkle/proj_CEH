# -*- coding: utf-8 -*-


import numpy as np
from wavelet import util
from eod import msg, mfg
import xarray as xr
import os
from wav_mixed import powerBlob_utils
import multiprocessing
from utils import u_grid, u_interpolate as u_int
import datetime as dt
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


        outt, nogood, t_thresh_size, t_thresh_cut, pix_nb = powerBlob_utils.filter_img(outt.values, 5)

        wav = util.waveletT(outt, dataset='METEOSAT5K_vera')
        core_min=-50
        outt[nogood] = np.nan

        power_msg = powerBlob_utils.find_scales_dominant(wav, nogood, dataset='mfg')

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
