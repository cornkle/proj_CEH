# -*- coding: utf-8 -*-


import numpy as np
import xarray as xr
import ipdb
import os
from utils import u_grid, u_interpolate as u_int
import glob
import CCores.cores as cores
import datetime
import multiprocessing
import pickle as pkl
from utils import constants as cnst
from eod import msg_bigDomain as msg


filepath = {

  #  'MSG_JJAS' : [cnst.network_data +'/data/OBS/MSG_WA30/', [6,7,8,9], (2019,2021)]  # 2004, 2018
    'MSG_JJAS' : ['/media/ck/Elements/Africa/grads_MSG_testfiles', [6], (2019,2021)]  # 2004, 2018

}

def run_all():

    for k in ['MFG_JJAS', 'MFG_MAMON', 'MSG_JJAS', 'MSG_MAMON']:
        run(k)



def run(dataset, CLOBBER=False):


    for yy in range((filepath[dataset])[2][0],((filepath[dataset])[2][1])+1):   # (2004,2016)

        for mm in (filepath[dataset])[1]:


            tag = dataset[0:3].upper()

            path =  '/media/ck/Elements/Africa/WestAfrica/cores_testdomain/' #'/prj/vera/cores_bigDomain/' #cnst.network_data + '/MCSfiles/MSG_cores/' #'/prj/vera/cores_bigDomain/' #
            savefile = path + 'coresPower_'+tag.upper()+'_-40_9-130km_-50points_dominant_'+str(yy) + '_'+str(mm).zfill(2)+'.nc'#'blobMap_-40-700km2_-50-points_dominant_'+str(yy) + '_'+str(mm).zfill(2)+'.nc'

            if not CLOBBER:
                if os.path.isfile(savefile):
                    print('File exists, continue!')
                    continue

            pool = multiprocessing.Pool(processes=6)
            print('Reading '+filepath[dataset][0])
            meteosat_folder = (filepath[dataset])[0]


            if tag == 'MSG':
                m = msg.ReadMsg(meteosat_folder, y1=yy, y2=yy, months=[mm])

            files  = m.fpath

            if len(files) == 0:
                continue


            mdic = m.read_data(files[0], llbox=[-25, 30, 2, 25])  #[-14, 2.5, 4, 11.5]

            # make salem grid
            grid = u_grid.make(np.arange(-18,29.9), np.arange(4,21), 3000)
            glon, glat = grid.ll_coordinates
            inds, weights, shape = u_int.interpolation_weights_grid(mdic['lon'].values, mdic['lat'].values, grid)
            gridd = (inds,weights,shape, grid)

            files_str = []

            for f in files:
                if tag == 'MSG':
                    files_str.append(f[0:-4])


            files_str = np.unique(files_str)

            passit = []
            for f in files_str:
                passit.append((gridd, m, f,tag))


            res = pool.map(_loop, passit)

            #res=[]
            #for l in passit:

            #    res.append(_loop(l))

            pool.close()

            res = [x for x in res if x is not None]

            try:
                ds = xr.concat(res, 'time')
            except ValueError:
                return

            try:
                os.remove(savefile)
            except OSError:
                pass
            # da.name = 'blob'
            # enc = {'blob': {'complevel': 5, 'zlib': True}}

            comp = dict(zlib=True, complevel=5)
            enc = {var: comp for var in ds.data_vars}

            ds.to_netcdf(path=savefile, mode='w', encoding=enc, format='NETCDF4')
            print('Saved ' + savefile)



def _loop(passit):

    gridd = passit[0]
    inds, weights, shape, grid = gridd

    m = passit[1]
    file = passit[2]
    tag = passit[3]


    if tag == 'MFG':
            print(' month')
    if tag == 'MSG':
        strr = file.split(os.sep)[-1]

        if ((strr[-2::]) != '00') & ((strr[-2::]) != '30'):
            print('Skip minute')
            return

    file = file + '.gra'
    print('Doing file: ' + file)
    try:
        mdic = m.read_data(file, llbox=[-25, 30, 2, 25])
    except FileNotFoundError:
        print('File not found')
        return

    if not mdic:
        print('File missing')
        return


    dat = mdic['t']

    hour = dat['time.hour']
    minute = dat['time.minute']
    day = dat['time.day']
    month = dat['time.month']
    year = dat['time.year']

    date = [datetime.datetime(int(year), int(month), int(day), int(hour), int(minute))]

    data = dat.squeeze().values
    #ipdb.set_trace()
    try:
        data = u_int.interpolate_data(data, inds, weights, shape)
    except IndexError:
        print('Interpolation problem, continue')
        return
    lon, lat = grid.ll_coordinates


    wObj = cores.dataset('METEOSAT3K_veraLS')
    wObj.read_img(data, lon, lat, edge_smoothing=False)
    wObj.applyWavelet()

    small_scale = wObj.scaleWeighting(wtype='nflicsv2')
    ss = wObj.to_dataarray(date=date, names=['small_scale', 'tir'])
    dominant = wObj.scaleWeighting(wtype='dominant')
    dom = wObj.to_dataarray(date=date, names=['dom', 'tir'])

    newds = xr.merge([ss,dom])

    #comp = dict(zlib=True, complevel=5)
    #enc = {var: comp for var in newds.data_vars}
    #
    #outfile = filepath.replace('108', 'wavelet')
    # # outfile = outfile.replace('real_time_data', 'real_time_wavelet')
    # #
    #path = os.path.dirname(outfile)
    #if not os.path.exists(path):
    #    os.makedirs(path)
    #
    #if os.path.isfile(outfile):
    #   print('File exists, continue')
    #   return
    # #
    #newds.to_netcdf(path=outfile, mode='w', encoding=enc, format='NETCDF4')
    #print('Saved ' + outfile)

    return newds






