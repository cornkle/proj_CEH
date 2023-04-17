import ipdb
import numpy as np
import xarray as xr
from utils import constants as cnst
from utils import u_parallelise as u_par
import matplotlib.pyplot as plt

import glob


def file_loop(chunk_in):

    dummy = xr.open_dataset(chunk_in[0])
    ds_dummy = xr.zeros_like(dummy)
    ds_mean = xr.zeros_like(dummy)
    ds_valid = xr.zeros_like(dummy)

    for ff in chunk_in:
        print('Doing', ff)
        ds = xr.open_dataset(ff)
        for dv in ds.data_vars:
            if dv not in ['lwout_noon', 'lsRain_noon', 'lw_out_PBLtop', 'lsRain']:
                new_mean = ds[dv] - ds[dv].mean()
                ds_dummy[dv] = ds_dummy[dv] + new_mean.where(np.isfinite(new_mean), other=0)
                # if dv in ['tcwv', 'shear']:
                #     f = plt.figure()
                #     plt.pcolormesh(ds_dummy[dv])
                #     plt.colorbar()
                #     plt.title('variable')
            else:
                ds_dummy[dv]= ds_dummy[dv] + (ds[dv].where(np.isfinite(ds[dv]), other=0))
            ds_mean[dv]= ds_mean[dv] + (ds[dv].where(np.isfinite(ds[dv]), other=0))
            ds_valid[dv].values = ds_valid[dv].values + np.isfinite(ds[dv]).astype(int)

            # if dv in ['tcwv', 'shear']:
            #     f = plt.figure()
            #     plt.pcolormesh(ds_valid[dv])
            #     plt.colorbar()
            #     plt.title('valid')


    anom = ds_dummy
    mean = ds_mean
    valid = ds_valid

    # f = plt.figure()
    # plt.pcolormesh(mean['shear']/valid['shear'])
    # plt.colorbar()
    # plt.title('full mean')

    return anom, mean, valid



nb_in_chunk = 2

TIME = 'hist'
if TIME == 'hist':
    tstr = 'historical'
else:
    tstr = 'future'

out = cnst.lmcs_drive + 'CP_models/MCS_files/MODELS/CP4_box/mean_loop/'
for m in range(6, 10):
    ds_files = glob.glob(
        cnst.lmcs_drive + 'CP_models/MCS_files/MODELS/CP4_box/CP4_allHours_'+tstr+'_5000km2_-50_WAf_box/*-' + str(
            m).zfill(2) + '-*_17:*.nc')
    print(len(ds_files))

    chunks = [ds_files[i:i + nb_in_chunk] for i in range(0, len(ds_files), nb_in_chunk)]

    # for cc in chunks:
    #     dic = file_loop(cc)
    #     ipdb.set_trace()

    dic = u_par.run_mixed(5, file_loop, chunks, ['anom', 'mean', 'valid'])

    outdic = {}
    for k in dic.keys():
        outdic[k] = xr.concat(dic[k], dim='chunks')

    anom = outdic['anom'].sum('chunks') / outdic['valid'].sum('chunks')
    std = (outdic['anom'] / outdic['valid']).std('chunks')

    mean = outdic['mean'].sum('chunks') / outdic['valid'].sum('chunks')
    stdm = (outdic['mean'] / outdic['valid']).std('chunks')

    anom.to_netcdf(out + 'CP4'+TIME+'_anom_mean_17h_' + str(m) + '.nc')
    std.to_netcdf(out + 'CP4'+TIME+'_anom_std_17h_' + str(m) + '.nc')

    mean.to_netcdf(out + 'CP4'+TIME+'_mean_17h_' + str(m) + '.nc')
    stdm.to_netcdf(out + 'CP4'+TIME+'_std_17h_' + str(m) + '.nc')