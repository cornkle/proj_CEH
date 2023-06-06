#!/usr/bin/env python
# coding: utf-8

# In[4]:


import ipdb
import numpy as np
import xarray as xr
from utils import u_grid, u_interpolate as u_int, constants as cnst, u_arrays, u_darrays, u_met
import pandas as pd
import datetime
import glob
from GLOBAL import glob_util
import pickle as pkl
import os
#import matplotlib.patches as patches


MREGIONS = {'WAf' : [[-18,25,4,25], 'spac', 0, (1,7), (8,12), (1,12), [-8,-7,10,11], 'WAf'], # last is hourly offset to UCT # 12    # [-18,25,4,25]
 'SAf' : [[20,35, -35,-15], 'spac', 2, (9,12), (1,5), (1,12), [21,24.5,-28,-24], 'SAf'], # 10
 'india_N' : [[70,90, 5,30], 'asia', 5, (1,7), (8,12), (1,12), [74,76,24,26], 'india'], # 7
 'india_S' : [[70,90, 5,30], 'asia', 5, (1,7), (8,12), (1,12), [76,79,18,21], 'india'], # 7
 'china_W' : [[105,115,25,40], 'asia', 8 , (1,7), (8,12), (1,12), [105,107,29,31], 'china'], # 4
 'china_E' : [[105,115,25,40], 'asia', 8 , (1,7), (8,12), (1,12), [111,113,25,27], 'china'], # 4
 'australia' : [[120,140,-23, -11], 'asia', 9, (10,12), (1,5), (1,12), [130,134, -21,-18], 'australia'], # 3
 'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4, (9,12), (1,5), (1,12), [-63,-60,-33,-30], 'sub_SA'] , # 16
# 'trop_SA' : [[-75, -50, -20, -5], 'spac', -5, (1,12), (1,12), (1,12)], # 17
 'GPlains_S' : [[-100,-90,32,47], 'nam', -6, (1,7), (8,12), (1,12), [-98,-95,36,39], 'GPlains'], # # 18
 'GPlains_N' : [[-100,-90,32,47], 'nam', -6, (1,7), (8,12), (1,12), [-99,-96,42,44], 'GPlains'] # # 18

}

def prepare_dates(s_region):

    REGION = MREGIONS[s_region][7]

    rdics = {}
    for regs in MREGIONS.keys():
        bregs = MREGIONS[s_region][7]
        for ids, y in enumerate(range(2012,2020)):
            test = pd.read_csv(cnst.lmcs_drive+'save_files/'+bregs+'_initTime__mcs_tracks_extc_'+str(y)+'0101_'+str(y)+'1231.csv')
            if ids == 0:
                test2 = pd.DataFrame(test)
            else :
                test2 = pd.concat([test2, test])
        rdics[regs] = test2

    tab = rdics[s_region]
    if tab[(tab['hour'] >= glob_util.LT_to_UTC_hour(16, REGION)) & (tab['hour'] <= glob_util.LT_to_UTC_hour(18, REGION))].size == 0:
        tab = tab[(tab['hour'] >= glob_util.LT_to_UTC_hour(16, REGION)) | (tab['hour'] <= glob_util.LT_to_UTC_hour(18, REGION))]
    else:
        tab = tab[(tab['hour'] >= glob_util.LT_to_UTC_hour(16, REGION)) & (tab['hour'] <= glob_util.LT_to_UTC_hour(18, REGION))]

    print(tab.meanlon.max(), tab.meanlon.min(), tab.meanlat.max(), tab.meanlat.min())
    box = MREGIONS[s_region][6]
    np.sum((tab.meanlon>box[0]) & (tab.meanlon<box[1]) & (tab.meanlat>box[2]) & (tab.meanlat<box[3]))

    tab = tab.loc[(tab.meanlon>box[0]) & (tab.meanlon<box[1]) & (tab.meanlat>box[2]) & (tab.meanlat<box[3])]

    lists = []
    for row in tab.iterrows():
        fdate = glob_util.UTC_to_LT_date(pd.to_datetime(row[1]['base_time']),  REGION)

        lists.append(datetime.datetime(fdate.year, fdate.month, fdate.day))
    tab['lt_date'] = lists
    print(np.unique(tab['lt_date']).size)
    return tab


def run(shift, s_region):

    dic = {
    'q' : [],
    'qup' : [],
    't' : [],
    'u' : [],
    'v' : [],
    'w' : [],
    'u_orig' : [],
    'v_orig' : [],        
    'd' : [],
    't2' : [],
    'u10' : [],
    'v10' : [],
    'u10_orig' : [],
    'v10_orig' : [],
    'slp' : [],
    'cnt' : [],
    'cape' : [],
    'divMoist' : [],
    'ice' : [],
    'ice_orig' : [],
   # 'sh' : [],
    'rh' : [],
    'ushear' : [],
    #'theta' : [],
    'sh' : [],
    'lh' : [],
    'sm' : [],
    'pr' : [],
    'vshear' : [],
    'shear' : []
    }

    tab = prepare_dates(s_region)
    REGION = MREGIONS[s_region][7]
    era5_files = cnst.lmcs_drive+ 'ERA5/hourly/' #cnst.ERA5 + 'hourly/
    af_pl = sorted(glob.glob(era5_files+'pressure_levels/'+REGION+'/'+'*.nc'))
    af_srfc = sorted(glob.glob(era5_files+'surface/'+REGION+'/'+'*.nc'))
    
    for date in tab['lt_date']:

        single = tab[tab['lt_date']==date]
        shour = 12

        lt_dt = pd.to_datetime(date)
        lt_dt = lt_dt.replace(hour=shour)
        dt = glob_util.LT_to_UTC_date(lt_dt, REGION)
        hour = dt.hour
        
        daystring = str(abs(shift))
        dayd = pd.Timedelta(daystring + 'days')
        if shift < 0:
            dt = dt - dayd
        if shift >=0:
            dt = dt + dayd
            
        # window1 = dt - pd.Timedelta('10days')
        # window2 = dt + pd.Timedelta('10days')

        fdate = '_' + str(dt.year) +'_' + str(dt.month).zfill(2) + '_' + str(dt.day).zfill(2)

        try:
            pl_file = era5_files + 'pressure_levels/'+REGION+'/ERA5' + fdate + '_'+REGION+'_pl.nc'
            srfc_file = era5_files + 'surface/'+REGION+'/ERA5' + fdate + '_'+REGION+'_srfc.nc'
          
            lsta = xr.open_dataset(pl_file)
            srfc = xr.open_dataset(srfc_file)
        except:
            #print('File missing', era5_files + 'pressure_levels/'+REGION+'/ERA5' + fdate + '_'+REGION+'_pl.nc')
            continue
            
        lsta = u_darrays.flip_lat(lsta)
        srfc = u_darrays.flip_lat(srfc)

        try:
            lsta_low = lsta.sel(time=str(dt.year)+'-'+str(dt.month)+'-'+str(dt.day)+'T'+str(hour).zfill(2)+':00:00', level=925)
        except:
            continue
        lsta_up = lsta.sel(time=str(dt.year)+'-'+str(dt.month)+'-'+str(dt.day)+'T'+str(hour).zfill(2)+':00:00', level=650)
        srfc_low = srfc.sel(time=str(dt.year)+'-'+str(dt.month)+'-'+str(dt.day)+'T'+str(hour).zfill(2)+':00:00')

        fpl_pos = af_pl.index(pl_file)
        fsrfc_pos = af_srfc.index(srfc_file)
        
        pl_toread = np.array(af_pl)[fpl_pos-10:fpl_pos+10]
        srfc_toread = np.array(af_srfc)[fsrfc_pos-10:fsrfc_pos+10]
        
        lsta10 = xr.open_mfdataset(pl_toread)
        lsta10 = u_darrays.flip_lat(lsta10)
        lsta10 = lsta10.sel(time=lsta10['time.hour']==hour).load() #, longitude=slice(-18, 0), latitude=slice(10,17)).load()   


        lsta10_low = lsta10.sel(level=925)
        lsta10_up = lsta10.sel(level=650)

        print('Climlen', len(srfc_toread))
        srfc10 = xr.open_mfdataset(srfc_toread)
        srfc10 = u_darrays.flip_lat(srfc10)
        srfc10 = srfc10.sel(time=srfc10['time.hour']==hour).load() # , longitude=slice(-18, 0), latitude=slice(10,17)).load() 

        q = lsta_low['q'].squeeze()
        qup = lsta_up['q'].squeeze()
        u = lsta_up['u'].squeeze()
        v = lsta_up['v'].squeeze()
        w = lsta_low['w'].squeeze()
        t = lsta_low['t'].squeeze()
        d = lsta_low['d'].squeeze()
        rh = lsta_low['r'].squeeze()
        t2 = srfc_low['t2m'].squeeze()
        u100 = srfc_low['u10'].squeeze()
        v100 = srfc_low['v10'].squeeze()
        slp = srfc_low['sp'].squeeze()
        cape = srfc_low['cape'].squeeze()
        divMoist = srfc_low['p84.162'].squeeze()
        ice = srfc_low['tciw'].squeeze()
        sh = srfc_low['msshf'].squeeze()
        lh = srfc_low['mslhf'].squeeze()
        sm = srfc_low['swvl1'].squeeze()
        pr = srfc_low['mtpr'].squeeze()


        
        q_clim = lsta10_low['q'].squeeze().mean('time')
        qup_clim = lsta10_up['q'].squeeze().mean('time')
        u_clim = lsta10_up['u'].squeeze().mean('time')
        v_clim = lsta10_up['v'].squeeze().mean('time')
        w_clim = lsta10_low['w'].squeeze().mean('time')
        t_clim = lsta10_low['t'].squeeze().mean('time')
        rh_clim = lsta10_low['r'].squeeze().mean('time')
        d_clim = lsta10_low['d'].squeeze().mean('time')
        t2_clim = srfc10['t2m'].squeeze().mean('time')
        u100_clim = srfc10['u10'].squeeze().mean('time')
        v100_clim = srfc10['v10'].squeeze().mean('time')
        slp_clim = srfc10['sp'].squeeze().mean('time')
        cape_clim = srfc10['cape'].squeeze().mean('time')
        divMoist_clim = srfc10['p84.162'].squeeze().mean('time')
        ice_clim = srfc10['tciw'].squeeze().mean('time')
        sh_clim = srfc10['msshf'].squeeze().mean('time')
        lh_clim = srfc10['mslhf'].squeeze().mean('time')
        sm_clim = srfc10['swvl1'].squeeze().mean('time')
        pr_clim = srfc10['mtpr'].squeeze().mean('time')

   

        print('Doing '+ 'AMSR_' + str(dt.year) + str(dt.month).zfill(2) + str(
        dt.day).zfill(2) + '.nc')
    
        # theta_low = u_met.theta_e(925,lsta_low['t'].squeeze().values-273.15,lsta_low['q'].squeeze().values)
        # theta_high = u_met.theta_e(650,lsta_up['t'].squeeze().values-273.15,lsta_up['q'].squeeze().values)
        # thetadiff = (theta_low-theta_high).squeeze()
        
        # theta10_low = u_met.theta_e(925,lsta10_low['t'].squeeze().mean('time').values-273.15,lsta10_low['q'].squeeze().mean('time').values)
        # theta10_high = u_met.theta_e(650,lsta10_up['t'].squeeze().mean('time').values-273.15,lsta10_up['q'].squeeze().mean('time').values)
        # thetadiff10 = (theta10_low-theta10_high).squeeze()
        
        cnt = np.zeros_like(q.values)
        cnt[np.isfinite(q.values)] = 1

        dic['q'].append(q.values- q_clim.values)
        dic['qup'].append(qup.values- qup_clim.values)
        dic['v'].append(v.values- v_clim.values)
        dic['w'].append(w.values- w_clim.values)
        dic['u'].append(u.values- u_clim.values)
        dic['rh'].append(rh.values- rh_clim.values)
        dic['v_orig'].append(v.values)#
        dic['u_orig'].append(u.values)#
        dic['t'].append(t.values-t_clim.values)
        dic['d'].append(d.values-d_clim.values)
        dic['t2'].append(t2.values-t2_clim.values)
        dic['u10'].append(u100.values-u100_clim.values)
        dic['v10'].append(v100.values-v100_clim.values)
        dic['u10_orig'].append(u100.values)#
        dic['v10_orig'].append(v100.values)#
        dic['slp'].append(slp.values-slp_clim.values)
        dic['cape'].append(cape.values-cape_clim.values)#s-v100_clim.values)
        dic['divMoist'].append(divMoist.values-divMoist_clim.values)#-slp_clim.values)
        dic['ice'].append(ice.values-ice_clim.values)
        dic['ice_orig'].append(ice.values)
#         ws, wd = u_met.u_v_to_ws_wd(u.values-u100.values, v.values-v100.values)
#         wsclim, wd = u_met.u_v_to_ws_wd(u_clim.values-u100_clim.values, v_clim.values-v100_clim.values)
        dic['ushear'].append((u.values-u100.values)-(u_clim.values-u100_clim.values)) #-wsclim
        dic['cnt'].append(cnt)
        dic['sh'].append(sh.values-sh_clim.values)
        dic['lh'].append(lh.values-lh_clim.values)
        dic['pr'].append(pr.values-pr_clim.values)
        dic['sm'].append(sm.values-sm_clim.values)
        dic['vshear'].append((v.values-v100.values)-(v_clim.values-v100_clim.values))
        vshear = (v.values-v100.values)-(v_clim.values-v100_clim.values)
        ushear = (u.values-u100.values)-(u_clim.values-u100_clim.values)
        dic['shear'].append(np.sqrt(vshear**2+ushear**2))

        lat = lsta_low.latitude.values
        lon = lsta_low.longitude.values
        

    for k in dic.keys():
        #print(k)
        dic[k] = np.nansum(np.stack(dic[k], axis=0), axis=0)

    return dic, lat, lon


def calc(dic):
    
    dics = {}
    for k in dic.keys():
        if k == 'cnt':
            continue
        dics[k] = dic[k] / dic['cnt']
    return dics


for s_region in MREGIONS.keys():
    if os.path.isfile(cnst.DATA+'LMCS/geogComp/'+s_region+'_'+'_ERA5_12LT.p'):
        print('File exists, continue!')

    dic, lat, lon = run(0, s_region)
    dic = calc(dic)
    dic['lat'] = lat
    dic['lon'] = lon

    pkl.dump(dic,open(cnst.DATA+'LMCS/geogComp/'+s_region+'_'+'_ERA5_12LT.p','wb'))


