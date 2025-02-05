import numpy as np
import xarray as xr
from utils import u_arrays as ua
import multiprocessing
import ipdb
import pickle as pkl
from utils import constants as cnst
from metpy import calc
from metpy.units import units
import glob
import pandas as pd
import os

TIMETAG = 'fut'
ANOMTAG = 'anom'
MAIN_PATH = '/gws/nopw/j04/lmcs/cklein/CP_models/MCS_files/WAf/CP4_box_JASMIN/mean3h_v2'

def dictionary(dummy):

    dic = {}
    for vn in dummy.data_vars:
        dic[vn+'_0.25deg'] = []  # 30km

    for vn in dummy.data_vars:
        dic[vn+'_1deg'] = []  # 1deg

    for vn in dummy.data_vars:
        dic[vn+'_2deg'] = []  # 2deg

    for vn in dummy.data_vars:
        dic[vn+'_Smean'] = [] # storm mean

    vars = ['hour', 'month', 'year', 'day',
            'lon', 'lat',
            'tmin', 'pmax', 'thetae_srfc', 'cape_proxy', 'cin_proxy', 'rh', 'div',
            'pgt30', 'pgt1', 'isvalid',
            'tgrad', 'SMgrad', 'SHgrad', 'LHgrad', 'tcwv_cl', 'sh_cl', 'lh_cl',
            'area', 'area-70', 'area-80', 'area-60', 'area_1mm', 'area_5mm', 'area_10mm']

    for v in vars:
        dic[v] = []
    return dic




def perSys():

    timetag = TIMETAG
    anom = ANOMTAG
    main_path = MAIN_PATH
    files = glob.glob(main_path+'/'+anom+'_'+timetag+'/*.nc')

    print('Nb files', len(files))
    # for y in range(1998,2007):
    #
    #     yfiles = []
    #     for f in files:
    #         if str(y)+'-' in f:
    #             yfiles.append(f)

    # pool = multiprocessing.Pool(processes=4)
    dummy = xr.open_dataset(files[0])
    #
    mdic = dictionary(dummy) #defaultdict(list)
    # print('NB files', len(files))
    # res = pool.map(file_loop, files)
    # pool.close()

    res=[]
    for f in files[0:50]:
        res.append(file_loop(f))

    keys = mdic.keys()
    for v in res:
        for k in keys:
            try:
                mdic[k].append(v[k])
            except TypeError:
                continue

    df = pd.DataFrame.from_dict(mdic, orient="index")
    df.to_csv(main_path+'/tables/'+timetag+'_table_'+anom+'_JASMIN_3hmeansVersion.csv')

        # pkl.dump(mdic, open(cnst.network_data +'data/LMCS/CP4_study_saves/bulk_CP4/bulk_'+tthresh+'_5000km2_CP4means_hourly_SAHEL_15kmprecip_WA_18W-25E_9-25N_-50C_LMCSfiles_17h_hist_'+str(y)+'.p',
        #                    'wb'))


def file_loop(f):
    print('Doing file: ' + f)
    ds = xr.open_dataset(f)
    out = dictionary(ds)



    outt = ds['lw_out_PBLtop'].values
    outp = ds['lsRain'].values
    outp_noon = ds['lsRain_noon'].values
    outp_noon[outp_noon!=0] = np.nan

    # if np.sum(np.isnan(outp_noon[54:62, 54:62])) > 0.1:
    #     print('Noon rainfall, continue')
    #     return

    outu_srfc = ds['u_srfc'].values
    outu_mid = ds['u_mid'].values
    outshear = ds['shear'].values
    outq = ds['q_srfc'].values
    tmid = ds['t_mid'].values
    tsrfc = ds['t_srfc'].values
    qmid = ds['q_mid'].values
    ####################################################
    pes = units.Quantity(650, 'hPa')
    pes_down = units.Quantity(925, 'hPa')
    tes_up = units.Quantity(tmid+273.15, 'K')
    qes = units.Quantity(outq/1000, 'kg/kg')
    qes_up = units.Quantity(qmid / 1000, 'kg/kg')
    tes_down = units.Quantity(tsrfc+273.15, 'K')

    rh = calc.relative_humidity_from_specific_humidity(pes_down, tes_down, qes)
    dew = calc.dewpoint_from_specific_humidity(pes_down, tes_down, qes)
    theta_e = calc.equivalent_potential_temperature(pes_down, tes_down,dew)
    theta_es = calc.saturation_equivalent_potential_temperature(pes, tes_up)

    dew_up = calc.dewpoint_from_specific_humidity(pes, tes_up, qes_up)
    theta_e_up = calc.equivalent_potential_temperature(pes, tes_up,dew_up)

    cape_proxy = (np.array(theta_e-theta_es))
    cin_proxy = (np.array(theta_e-theta_e_up))

   ######################################################

    u = units.Quantity(ds['u_srfc'].values, 'm/s')
    v = units.Quantity(ds['v_srfc'].values, 'm/s')

    dx = units.Quantity(4400, 'm')
    div = calc.divergence(u, v, dx=dx, dy=dx)
   ########################################################


    out['lon'] = ds.minlon
    out['lat'] = ds.minlat
    out['hour'] = ds['time.hour'].item()
    out['month'] = ds['time.month'].item()
    out['year'] = ds['time.year'].item()
    out['day'] = ds['time.day'].item()


    t_thresh = -50  # -40C ~ 167 W m-2

    mask = (outt<=t_thresh)
    varmask = (outt<=-10) & (np.isfinite(outp_noon))
    if np.sum(mask)<=3:
        return

    rainfield = outp[mask]

    out['area'] = np.sum((outt<=t_thresh))
    out['area-60'] = np.sum((outt <= -60))
    out['area-70'] = np.sum((outt <= -70))
    out['area-80'] = np.sum((outt <= -80))
    out['area_1mm'] = np.sum((outp >=1))
    out['area_5mm'] = np.sum((outp>=5))
    out['area_10mm'] = np.sum((outp>=10))

    maxpos = np.unravel_index(np.nanargmax(outp), outp.shape)
    minpos = np.unravel_index(np.nanargmin(outt), outt.shape)

    out['tmin'] = np.nanmin(outt)
    out['pmax'] = np.nanmean(ua.cut_kernel(outp, maxpos[1], maxpos[0], 1))  # degrade rainfall to 15km
    out['pmax_native'] = np.max(outp[mask])  # rain at 4.4km
    #######
    if 'anom' in f:
        basefile = os.path.basename(f)
        timetag = TIMETAG

        dpath = MAIN_PATH+'/mean_'+timetag+'/*.nc'
        try:
            dcl = xr.open_dataset(dpath + basefile)
            out['tcwv_cl'] = np.nanmean(ua.cut_kernel(dcl['tcwv'].values, minpos[1], minpos[0], 11))
            out['sh_cl'] = np.nanmean(ua.cut_kernel(dcl['sh'].values, minpos[1], minpos[0], 11))
            out['lh_cl'] = np.nanmean(ua.cut_kernel(dcl['lh'].values, minpos[1], minpos[0], 11))
            del dcl
        except:
            return
            #out['tcwv_cl'] = np.nan
   ########
    for boxes in ((11, '_1deg'), (22, '_2deg'), (3, '_0.25deg'), (mask, '_Smean')):

        dist = boxes[0]
        tag = boxes[1]
        yy = minpos[0]
        xx = minpos[1]

        for vn in ds.data_vars:

            if isinstance(dist, int):
                cnt = ds[vn].where(varmask).sel(latitude=slice(dist*-1, dist), longitude=slice(dist*-1, dist)).count(['longitude', 'latitude'])
                if cnt > 0.5 * ((dist*2)**2):
                    data = float(ds[vn].sel(latitude=slice(dist*-1, dist), longitude=slice(dist*-1, dist)).mean(['longitude', 'latitude']).values)  # box of 1deg or 0.25deg
                else:
                    data = np.nan
            else:
                data = float(ds[vn].where(mask).mean(['longitude', 'latitude']).values)  # storm mask
            out[vn + tag] = (float(data))

    out['rh'] = np.nanmean(ua.cut_kernel(rh,minpos[1], minpos[0],11)) # 1deg RH
    out['thetae_srfc'] = np.nanmean(ua.cut_kernel(theta_e, minpos[1], minpos[0], 11))
    out['cape_proxy'] = np.nanmean(ua.cut_kernel(cape_proxy, minpos[1], minpos[0], 11))
    out['cin_proxy'] = np.nanmean(ua.cut_kernel(cin_proxy, minpos[1], minpos[0], 11))
    out['div'] = np.nanmean(ua.cut_kernel(div, minpos[1], minpos[0], 11))

    tgrad_lat = ds.isel(longitude=slice(minpos[1] - 11, minpos[1] + 11)).mean('longitude').squeeze() #1deg strip for gradients
    tgrad = tgrad_lat.polyfit(dim='latitude', deg=1)

    out['tgrad'] = float((tgrad['t2_polyfit_coefficients'])[0])
    try:
        out['SMgrad'] = float((tgrad['SM_polyfit_coefficients'])[0])
    except:
        out['SMgrad'] = np.nan
    out['SHgrad'] = float((tgrad['sh_polyfit_coefficients'])[0])
    out['LHgrad'] = float((tgrad['lh_polyfit_coefficients'])[0])

    out['pgt30'] = np.sum(rainfield>30)
    out['isvalid'] = np.sum(mask)
    out['pgt1'] = np.sum(rainfield>1)

    ds.close()

    return out
