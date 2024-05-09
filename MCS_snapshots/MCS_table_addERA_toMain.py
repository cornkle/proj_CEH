import xarray as xr
import ipdb
import pandas as pd
import glob
import os
from utils import constants as cnst
from GLOBAL import glob_util
import datetime
import numpy as np
import multiprocessing
from utils import u_darrays


MREGIONS = {
'GPlains' : [[-100,-90,32,47], 'nam', -6, (1,7), (5,9), (1,12)], # # 18
'china' : [[105,115,25,40], 'asia', 8 , (1,7), (5,9), (1,12)], # 4
'india' : [[70,90, 5,30], 'asia', 5, (1,7), (5,9), (1,12)], # 7
'WAf' : [[-18,25,4,25], 'spac', 0, (1,7), (5,9), (1,12)], # last is hourly offset to UCT # 12    # [-18,25,4,25]
'australia' : [[120,140,-23, -11], 'asia', 9, (11,3), (11,3), (1,12)], # 3
'SAf' : [[20,35, -35,-15], 'spac', 2, (9,12), (11,3), (11,3)], # 10
'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4, (11,3), (11,3), (1,12)] , # 16
'SA_big' : [[9.5, 52, -36,-5.5], 'spac', 2, (9,12), (11,3), (11,3)],
}
####################### READS 5000km2 files that include TIR/IMERG from 2d fields to add ERA5 here.
def ERA_dictionary(tab):
    dic = {}
    vars = [
        
        'tcwv', 'tgrad2m', 'tgrad2m_resid', 'tgrad925', 'tgrad850', 'smgrad', 'efgrad_acc', 'shgrad', 'shgrad_resid', 'lhgrad', 'shgrad_acc', 'lhgrad_acc',
        'q925', 'q650', 'q850', 'era_precip', 'era_precip_acc', 'sm', 'ef_acc', 'sh_acc', 'lh_acc',
        'u925', 'u650',
        'v925', 'v650',
        'w925', 'w650',
        'rh925', 'rh650',
        't925', 't650', 't850',
        'div925', 'div850', 'div650',
        'pv925', 'pv650','pv850',
        'ushear925_650', 'ushear850_650', 'ushear925_850', 'ushear100m_650',
        'vshear925_650', 'vshear850_650', 'vshear925_850', 'vshear100m_650',
        'shear925_650', 'shear850_650', 'shear925_850', 'shear100m_650',
        'cape', 't2m']

    for v in vars:
        dic[v] = []

    dummy = tab
    for v in dummy.keys():
        dic[v] = []


    return dic




def add_ERA_regional(lt_mcs_hours):

    intab = cnst.lmcs_drive + '/MCS_5000km2_tables/'
    lt_hour_start = lt_mcs_hours[0]
    lt_hour_end = lt_mcs_hours[1]
    mreg = REGION
    for yy in range(2000,2021):
            filename = glob.glob(intab+mreg+'/'+str(yy)+'_MCS_5000km2_*.csv')[0]
            outfile = filename.replace(mreg, 'ERA5_added/'+mreg)
            outfile = outfile.replace('.csv', '_'+mreg+'_ERA5.csv')

            # try:
            #     os.remove(outfile)
            #     print('Removed file')
            # except:
            #     pass

            if os.path.isfile(outfile):
                print('File exists, continue')
                continue

            print('Doing', yy)

            try:
                msg = pd.read_csv(filename)
            except:
                ipdb.set_trace()

            msg = msg.to_dict(orient='list')
            msg['lt_hour'] = []
            msg['lt_datetime'] = []
            msg['utc_date'] = []
            for utc_hour in msg['hour']:
                msg['lt_hour'].append(glob_util.UTC_to_LT_hour(utc_hour,mreg))
            for y, m, d, hh in zip(msg['year'], msg['month'], msg['day'], msg['hour']):
                udatet = datetime.datetime(y, m, d,hh)
                msg['lt_datetime'].append(glob_util.UTC_to_LT_date(udatet, mreg))
                msg['utc_date'].append(udatet.date())

            hour_mask = (np.array(msg['lt_hour'])>= lt_hour_start) & (np.array(msg['lt_hour'])<= lt_hour_end)
            mask = np.where(hour_mask)
            for k in msg.keys():
                msg[k] = np.array(msg[k])[mask]

            ubt = np.unique(msg['utc_date'])
            #msg['utc_date'] = np.array(msg['utc_date'])
            chunks = []
            for ut in ubt:
                daydir = {}

                pos = np.where(msg['utc_date'] == ut)  # [0]).astype(int)
                # ipdb.set_trace()
                for k in msg.keys():
                    daydir[k] = np.array(msg[k])[pos]
                chunks.append(daydir)


            pool = multiprocessing.Pool(processes=5)

            res = pool.map(run_ERA5_regional, chunks)
            pool.close()

            #res = []
            #for f in chunks[0:10]:
           # 
            #    out = run_ERA5_regional(f)
            #    res.append(out)
            #

            mdic = ERA_dictionary(msg)  # defaultdict(list)
            print('Back from multiproc')
            keys = mdic.keys()
            for vv in res:
                if vv == None:
                    continue
                for k in keys:
                                                          
                    try:
                        mdic[k].extend(vv[k])
                    except TypeError:
                        ipdb.set_trace()
                        

            for kk in mdic.keys():
                print(kk, len(mdic[kk]))
            try:
                df = pd.DataFrame.from_dict(mdic)
            except:
                print('Nb problem')
                for k in mdic.keys():
                    print(k, len(mdic[k]))
                    ipdb.set_trace()
            df.to_csv(outfile, index=False)
            print('Saved ', outfile)


def run_ERA5_regional(fi):
    #run per daily chunks
    inpath = cnst.lmcs_drive + '/ERA5/hourly/'
    out = ERA_dictionary(fi)

    dat = fi['lt_datetime'][0]
    mcs_local_time = datetime.datetime(dat.year, dat.month, dat.day)
    era_hour = 12
    etime_local = mcs_local_time.replace(hour=era_hour, minute=0)
    edate = glob_util.LT_to_UTC_date(etime_local, REGION)
    print(REGION, 'MCS_utc_time:', fi['date'][0], 'MCS_lt_time', dat, 'ERA_sampling_utc', edate)
    try:
        era_pl = xr.open_dataset(
            inpath + 'pressure_levels/' + REGION + '/ERA5_' + str(edate.year) + '_' + str(edate.month).zfill(
                2) + '_' + str(edate.day).zfill(2) + '_' + REGION + '_pl.nc')
    except:
        print('ERA5 missing')
        return
    try:
        era_srfc = xr.open_dataset(
            inpath + 'surface/' + REGION + '/ERA5_' + str(edate.year) + '_' + str(edate.month).zfill(
                2) + '_' + str(edate.day).zfill(2) + '_' + REGION + '_srfc.nc')
    except:
        print('ERA5 srfc missing')
        return

    try:
        era_wi100 = xr.open_dataset(
            inpath + 'surface/'+REGION+'/100mWind_ERA5_' + str(edate.year) + '_' + str(edate.month).zfill(
                2) + '_'+ str(edate.day).zfill(2) + '_'+REGION+'_srfc.nc')
    except:
        print('ERA5 missing')
        return


    era_pl = u_darrays.flip_lat(era_pl)
    era_srfc = u_darrays.flip_lat(era_srfc)
    era_wi100 = u_darrays.flip_lat(era_wi100)

    era_pl_day = era_pl.sel(time=edate, method='nearest')
    era_srfc_day = era_srfc.sel(time=edate, method='nearest')
    era_wi100_day = era_wi100.sel(time=edate, method='nearest')
    era_accum = era_srfc.sel(time=slice(edate-datetime.timedelta(hours=4),edate)).mean('time')
    # if isinstance(fi['pf_maxrainrate1'], int):
    #     pr = np.array([fi['pf_maxrainrate1']])
    # else:
    pr = fi['precipitation_max']

    print('Nb per day', pr.size)

    for ids in range(pr.size):

        for fikey in fi.keys():
            try:
                out[fikey].append(fi[fikey][ids])
            except:
                ipdb.set_trace()

        elat = fi['tminlat'][ids]
        elon = fi['tminlon'][ids]
        try:
            era_day = era_pl_day.sel(latitude=slice(elat - 0.35, elat + 0.35), longitude=slice(elon - 0.35, elon + 0.35)).mean(['latitude', 'longitude'])
        except:
            ipdb.set_trace()
        era_day_srfc = era_srfc_day.sel(latitude=slice(elat - 0.35, elat + 0.35), longitude=slice(elon - 0.35, elon + 0.35)).mean(['latitude', 'longitude'])
        era_day_wi100 = era_wi100_day.sel(latitude=slice(elat - 0.35, elat + 0.35), longitude=slice(elon - 0.35, elon + 0.35)).mean(['latitude', 'longitude'])
        era_day_accum = era_accum.sel(latitude=slice(elat - 0.35, elat + 0.35), longitude=slice(elon - 0.35, elon + 0.35)).mean(['latitude', 'longitude'])

        tgrad_lat = era_srfc_day.sel(latitude=slice(elat - 2, elat + 2), longitude=slice(elon - 0.35, elon + 0.35)).mean('longitude').squeeze()
        tgrad = tgrad_lat.polyfit(dim='latitude', deg=1)
        print('tgradlen', era_srfc_day['t2m'].sel(latitude=slice(elat - 2, elat + 2), longitude=slice(elon - 0.35, elon + 0.35)).mean('longitude').squeeze().size)

        tgrad_pl = era_pl_day['t'].sel(latitude=slice(elat - 2, elat + 2), longitude=slice(elon - 0.35, elon + 0.35)).mean('longitude').squeeze().polyfit(dim='latitude', deg=1)
        

        tgrad_lat_accum = era_accum.sel(latitude=slice(elat - 2, elat + 2), longitude=slice(elon - 0.35, elon + 0.35)).mean('longitude').squeeze()
        tgrad_accum = tgrad_lat_accum.polyfit(dim='latitude', deg=1)


        ef_lat = tgrad_lat_accum['mslhf'] / (tgrad_lat_accum['mslhf'] + tgrad_lat_accum['msshf'])
        ef_poly = ef_lat.squeeze().polyfit(dim='latitude', deg=1)
        ef_mean = era_accum['mslhf'].mean() / (era_accum['mslhf'].mean() + era_accum['msshf']).mean()
        e925 = era_day.sel(level=925).mean()
        e650 = era_day.sel(level=650).mean()
        e850 = era_day.sel(level=850).mean()
        srfc = era_day_srfc.mean()
        wi100 = era_day_wi100.mean()
        accum = era_day_accum.mean()

        del era_day
        del era_day_srfc

        try:
            out['q925'].append(float(e925['q']))
        except TypeError:
            return

        out['q650'].append(float(e650['q']))
        out['v925'].append(float(e925['v']))
        out['v650'].append(float(e650['v']))  # this was wrong in v2 dataset previously, saved is v925 for both levels 650 and 925
        out['u925'].append(float(e925['u']))
        out['u650'].append(float(e650['u']))
        out['w925'].append(float(e925['w']))
        out['w650'].append(float(e650['w']))
        out['rh925'].append(float(e925['r']))
        out['rh650'].append(float(e650['r']))
        out['t925'].append(float(e925['t']))
        out['t650'].append(float(e650['t']))
        out['t850'].append(float(e850['t']))
        out['pv925'].append(float(e925['pv']))
        out['pv650'].append(float(e650['pv']))
        out['pv850'].append(float(e850['pv']))
        out['div925'].append(float(e925['d']))
        out['div650'].append(float(e650['d']))
        out['div850'].append(float(e850['d']))
        out['q850'].append(float(e850['q']))
        out['tcwv'].append(float(srfc['tcwv']))
        out['cape'].append(float(srfc['cape']))
        out['t2m'].append(float(srfc['t2m']))
        out['era_precip'].append(float(srfc['mtpr']) * 3600)
        out['era_precip_acc'].append(float(accum['mtpr']) * 3600) # since 8am
        out['sm'].append(float(srfc['swvl1']))
        out['ef_acc'].append(float(ef_mean))
        out['sh_acc'].append(float(accum['msshf'])*-1)
        out['lh_acc'].append(float(accum['mslhf'])*-1)
        out['tgrad2m'].append(float(tgrad['t2m_polyfit_coefficients'][0]))
        out['tgrad925'].append(float(tgrad_pl.sel(level=925)['polyfit_coefficients'][0]))
        out['tgrad850'].append(float(tgrad_pl.sel(level=850)['polyfit_coefficients'][0]))
        
        out['shgrad_acc'].append(float(tgrad_accum['msshf_polyfit_coefficients'][0]) * -1)
        out['lhgrad_acc'].append(float(tgrad_accum['mslhf_polyfit_coefficients'][0]) * -1)  # era fluxes opposite sign convention
        out['efgrad_acc'].append(float(ef_poly['polyfit_coefficients'][0]))

        out['shgrad'].append(float(tgrad['msshf_polyfit_coefficients'][0]) * -1)
        out['lhgrad'].append(float(tgrad['mslhf_polyfit_coefficients'][0]) * -1)  # era fluxes opposite sign convention
        out['smgrad'].append(float(tgrad['swvl1_polyfit_coefficients'][0]))
       
        out['shgrad_resid'].append(float(tgrad_accum['msshf_polyfit_coefficients'][1]))
        out['tgrad2m_resid'].append(float(tgrad['t2m_polyfit_coefficients'][1]))


        u9256 = float(e650['u'] - e925['u'])
        u8506 = float(e650['u'] - e850['u'])
        u8509 = float(e850['u'] - e925['u'])
        u100 = float(e650['u'] - wi100['u100'])

        v9256 = float(e650['v'] - e925['v'])
        v8506 = float(e650['v'] - e850['v'])
        v8509 = float(e850['v'] - e925['v'])
        v100 = float(e650['v'] - wi100['v100'])

        out['ushear925_650'].append(u9256)
        out['ushear850_650'].append(u8506)
        out['ushear925_850'].append(u8509)
        out['ushear100m_650'].append(u100)

        out['vshear925_650'].append(v9256)
        out['vshear850_650'].append(v8506)
        out['vshear925_850'].append(v8509)
        out['vshear100m_650'].append(v100)

        out['shear925_650'].append(np.sqrt(u9256 ** 2 + v9256 ** 2))
        out['shear850_650'].append(np.sqrt(u8506 ** 2 + v8506 ** 2))
        out['shear925_850'].append(np.sqrt(u8509 ** 2 + v8509 ** 2))
        out['shear100m_650'].append(np.sqrt(u100 ** 2 + v100 ** 2))
     

    return out


for regs in MREGIONS:
    print('DOing', regs)
    
    REGION = regs
    print('Main region', REGION)
    MONTHS = (MREGIONS[REGION])[4]
    MHOUR_SLICE = (15,21)
    add_ERA_regional(MHOUR_SLICE)




