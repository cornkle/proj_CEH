import numpy as np
from scipy.ndimage import label
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
from JASMIN import MetUM_variables as mu

def load_file(ffile,var):
    try:
        ds=xr.open_dataset(ffile)[var]
    except:
        return
    ds=ds.assign_coords({"longitude":(ds.longitude-360)})

    try:
        ds=ds.isel(pressure=slice(None,None,-1))
    except:
        pass
    return ds


def olr_to_bt(olr):
    #Application of Stefan-Boltzmann law
    sigma = 5.670373e-8
    tf = (olr/sigma)**0.25
    #Convert from bb to empirical BT (degC) - Yang and Slingo, 2001
    a = 1.228
    b = -1.106e-3
    Tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return Tb - 273.15


def filtering(dar, v, outv, pl, filepath, box, h):

    
    if (v == 'lsRain') | (v == 'totRain'):
        dar.values = dar.values*3600  # rain to mm/h
        dar.attrs['units'] = 'mm h-1'

    if (v == 'SM'):
        try:
            dar = dar.sel(depth=0.05)
        except:
            pass
        dar = dar.where(dar < 500, other=np.nan)

    if (v == 'lw_out_PBLtop'):
        dar.values = olr_to_bt(dar.values)

    if 'pressure' in dar.coords:

        try:
            dar.values[dar.values == 0] = np.nan  # potential missing value maskout
        except ValueError:
            ipdb.set_trace()
            print('Pressure value error!')
            return
        if (len(pl) > 1) & (outv == 'shear'):

            shear = dar.sel(pressure=650).values - dar.sel(pressure=925).values
            dar = dar.sum(dim='pressure').squeeze()
            dar.values = shear

    if (outv == 'tcwv'):
        nfilepath = filepath.replace(v, 'colDryMass')
        nfilepath = nfilepath.replace(mu.create_CP4_filename(v), mu.create_CP4_filename('colDryMass'))
        try:
            dar2 = load_file(nfilepath, mu.create_CP4_filename('colDryMass'))
        except:
            print('Cannot find tcw dry mass, return')
           
            return
        dar2 = dar2.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
        dar2 = dar2[dar2['time.hour'] == h].squeeze()
        dar = dar-dar2
        del dar2

    return dar



def file_save(cp_dir, out_dir, vars, datestring, box, tthresh, pos):

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

        orig_v = mu.create_CP4_filename(v)

        try:
            filepath = glob.glob(cp_dir+os.sep+str(v)+os.sep+'*_'+datestring+'*.nc')[0]
        except IndexError:
            
            if (v=='pblH') & ('fut' in cp_dir):
                print('PBLH not available for future, next variable!')
                continue
            else:
                print('No file found, return')
                return

        try:
            arr = load_file(filepath, orig_v)
        except OSError:
            print('Cannot open file, continue! ', filepath)
            return
        dar = arr.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
        if np.sum(pl) != 0:
            dar = dar.sel(pressure=dar.pressure.isin(pl))
            if len(dar.pressure.values) != len(plevs):
                ipdb.set_trace()
                print('Not correct number of pressure levels,return')
                return
        del arr

        if outv == 'lsRain_noon':
            dum = dar[dar['time.hour'] == h].squeeze()
            dar = dar[(dar['time.hour'] <= h) & (dar['time.hour'] >= h-2)].sum('time').squeeze() # early rain accumulation over 3h
            dar = dar.assign_coords(coords={'time' : dum.time})
            del dum
        else:
            dar = dar[dar['time.hour'] == h].squeeze()

        dar = filtering(dar, v, outv, pl,filepath,box,h)

        ddate = dar.time.values.item()
        dt = datetime.datetime(ddate.year, ddate.month, ddate.day,
                               ddate.hour, ddate.minute, ddate.second)
        
        if ATAG == 'anom' :

            if "hist" in filepath:
                ctag = 'historical'
            if "fut" in filepath:
                ctag = 'future'
            
            clim_files = '/gws/nopw/j04/lmcs/cklein/CP_models/MCS_files/climatology/'+ctag
            
            ndatestring = str(ddate.month).zfill(2)+'-'+str(ddate.day).zfill(2)+'_'+str(h)

            nfile = clim_files+os.sep+str(v)+'_'+ndatestring+'_'+ctag+'.nc'

            try:
                meanarr = xr.open_dataset(nfile)[mu.create_CP4_filename(v)]
                mmeanarr = meanarr.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
                if np.sum(pl) != 0:
                    mmeanarr = mmeanarr.sel(pressure=mmeanarr.pressure.isin(pl))
                    if len(dar.pressure.values) != len(plevs):
                        print('Not correct number of pressure levels,return')
                        ipdb.set_trace()
                        return
            except OSError:
                print('Couldnt find clim file, return')
                return
               
            mmeanarr = filtering(mmeanarr, v, outv, pl, filepath, box, h)
            dar = (dar - mmeanarr).squeeze()
            print('Did anom calc!', dar.shape)

        
        if grid == 'srfc':
            try:
                regrid = u_int.interpolate_data(dar.values, inds, weights, shape)

            except ValueError:
                ipdb.set_trace()


            da = xr.DataArray(regrid,
                              coords={'time': dar.time, 'latitude': pl_dummy.latitude.values,
                                      'longitude': pl_dummy.longitude.values, },
                              dims=['latitude', 'longitude'])
            da.values[pos[0], pos[1]] = np.nan  # mask sea

        else:
            da = xr.DataArray(dar.values,
                              coords={'time': dar.time, 'pressure': dar.pressure, 'latitude': pl_dummy.latitude.values,
                                      'longitude': pl_dummy.longitude.values, },
                               dims=['pressure','latitude', 'longitude'])
            da.values[:, pos[0], pos[1]] = np.nan  # mask sea
            
        da.attrs = dar.attrs

        if (v == 'lw_out_PBLtop') & (idx == 0):

            da.values[da.values >= tthresh] = 0  # T threshold maskout
            da.values[np.isnan(da.values)] = 0 # set ocean nans to 0

            date = pd.Timestamp(datestring)
            date = date.replace(hour=h)

            labels, numL = label(da.values)

            u, inv = np.unique(labels, return_inverse=True)
            n = np.bincount(inv)

            ## labels and filtering not needed if object location table is used for identification!
            goodinds = u[n >= 258]  # 258 == 5000km2 area threshold used here. defines minimum MCS size e.g. 350 km2 = 39 pix at 3x3km res (258 pix at 4.4km is 5000km2) 52 pix is 1000km2 for cp4
            if not sum(goodinds) > 0:
                print('No goodinds!')
                return

        ds[outv] = da

        print('Saved ', outv, h)

    for gi in goodinds:
        if (gi == 0):  # index 0 is always background, ignore!
            continue
        mask = np.where(labels!=gi)
        dbox = ds.copy(deep=True)

        for vout in dbox.data_vars:
            if vout not in ['lsRain','lw_out_PBLtop']:
                continue
            (dbox[vout].values)[mask] = np.nan

        try:
            if np.nansum(dbox['lsRain'])<=0.1:
                print('No rain, return')
                return
        except KeyError:
            if np.nansum(dbox['totRain'])<=0.1:
                print('No rain, return')
                return
        
        filt = dbox
        # Optional rain and cloud filter, keep off for now!
        # filt = dbox.where((dbox['lsRain_noon'] < 0.005) & (dbox['lwout_noon'] > -30))
        # for raw_var in ['lsRain', 'lw_out_PBLtop']:
        #     filt[raw_var] = dbox[raw_var]

        # location of minimum cloud top temperature, could come from table
        tmin = filt.where(filt['lw_out_PBLtop'] == filt['lw_out_PBLtop'].min(), drop=True)

        # use this point location when minlon, minlat are saved in table as intermediate step - lon,lat will not correspond exactly to use "==" method
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
            xpos = xpos[0][0]
            ypos = ypos[0][0]
        except TypeError:
            continue

        distxx = 57  # 57 = 250 km at 4.4 res, 500km across
        distyx = 6   # 6 = 25 km at 4.4 res, 50km across
        try:
            filtx = filt.isel(latitude=slice(ypos - distyx, ypos + distyx + 1),
                             longitude=slice(xpos - distxx, xpos + distxx + 1)).mean('latitude').squeeze()
        except IndexError:
            continue
        if (len(filtx.longitude) != distxx * 2 + 1):
            print(filt)
            continue

        #if 10% missing in cross section, don't save
        if (np.sum(np.isnan(filtx['u_cross']))/filtx['u_cross'].size) > 0.1:
            continue

        distyy = 57  # 57 = 250 km at 4.4 res, 500km across
        distxy = 6   # 6 = 25 km at 4.4 res, 50km across
        try:
            filty = filt.isel(latitude=slice(ypos - distyy, ypos + distyy + 1),
                             longitude=slice(xpos - distxy, xpos + distxy + 1)).mean('longitude').squeeze()
        except IndexError:
            continue
        if (len(filty.latitude) != distyy * 2 + 1):
            print(filt)
            continue

        filtx = filtx.assign_coords({'longitude': np.arange(distxx * -1, distxx + 1)})
        filty = filty.assign_coords({'latitude': np.arange(distyy * -1, distyy + 1)})


        #filt = filt.drop_vars(['lsRain_noon', 'lwout_noon'])
        filt.attrs = {'minlat' : plat, 'minlon' : plon}

        #if np.nansum(ds_box['lwout_noon'])<=0:
        #   ipdb.set_trace()
        #   print('lwout is lt 0', np.nansum(ds_box['lwout_noon']))
        #   return
        
        savefile_x = out_dir + os.sep + date.strftime('%Y-%m-%d_%H:%M:%S') + '_XDIR_' + str(gi) + '_lonXlat_'+str(np.round(plon,1))+'_'+str(np.round(plat,1))+'.nc'
        savefile_y = out_dir + os.sep + date.strftime('%Y-%m-%d_%H:%M:%S') + '_YDIR_' + str(gi) + '_lonXlat_'+str(np.round(plon,1))+'_'+str(np.round(plat,1))+'.nc'

        try:
            os.remove(savefile_x)
        except OSError:
            pass
        try:
            os.remove(savefile_y)
        except OSError:
            pass

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in filt.data_vars}
        try:
            filtx.to_netcdf(path=savefile_x, mode='w', encoding=encoding)
        except:
            ipdb.set_trace()
        filty.to_netcdf(path=savefile_y, mode='w', encoding=encoding)
        print('Saved ' + savefile_x)
        print('Saved MCS no.'+str(gi)+ ' as netcdf.')

        del filt
        del dbox
    del ds
    del da


### Inputs:
##### Provide hour when you run script!!! as sys.argv[1]
fdir = str(sys.argv[2])
if fdir == 'hist':
    ftag = 'hist'
if fdir == 'future':
    ftag = 'fut'

ATAG = str(sys.argv[3]) # anom or mean

main = '/home/users/cornkle/linked_CP4/'
main_lmcs = '/gws/nopw/j04/lmcs/cklein/CP_models/MCS_files/'
data_path = main + '/'+fdir
ancils_path = '/home/users/cornkle/impala/shared/CP4A/ncfiles/4km/ANCILS/'
out_path = main_lmcs + 'CP4_box_JASMIN/pl_'+ATAG+'_'+ftag+'/'
box = [-18, 25, 5, 25]  # W- E , S - N geographical coordinates box
MINLAT = 8

years = np.array(np.arange(1998,2007), dtype=str)
months = (['07','08', '09'])
days = np.array(np.arange(1,32), dtype=str)

tthresh = -50 # chosen temperature threshold, e.g. -50, -60, -70
h= int(sys.argv[1]) # hour to extract

plglob = glob.glob(data_path + '/q_pl/*.nc')
pl_dummy = load_file(plglob[0], mu.create_CP4_filename('q_pl'))

srfcglob = glob.glob(data_path + '/lw_out_PBLtop/*.nc')
srfc_dummy =load_file(srfcglob[0], mu.create_CP4_filename('lw_out_PBLtop'))

## getting interpolation weights for later use
pl_dummy = pl_dummy.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
srfc_dummy = srfc_dummy.sel(longitude=slice(box[0],box[1]), latitude=slice(box[2],box[3]))
# load seamask
landsea_path = glob.glob(ancils_path + os.sep + 'landseamask*.nc')[0]
landsea = xr.open_dataset(landsea_path, decode_times=False)
ls = landsea['lsm']

ls = ls.assign_coords(rlon = ls.rlon.values - 360)
ls_arr = ls.sel(rlon=slice(box[0], box[1]), rlat=slice(box[2], box[3]))


pos = np.where(ls_arr[0, 0, :, :] == 0)
lons, lats = np.meshgrid(pl_dummy.longitude.values, pl_dummy.latitude.values)#np.meshgrid(ls_arr.rlon.values, ls_arr.rlat.values)
inds, weights, shape = u_int.interpolation_weights(srfc_dummy.longitude, srfc_dummy.latitude, pl_dummy.longitude, pl_dummy.latitude)

plevs = [950,925,900,850,800,750,700,650,550,500]

vars = OrderedDict()   # dictionary which contains info on pressure level and hour extraction for wanted variables
vars['lw_out_PBLtop'] = ([], h, (inds,weights,shape), 'lw_out_PBLtop', 'srfc')  ### Input in BRIGHTNESS TEMPERATURES!! (degC)
vars['lsRain'] =  ([], h, (inds,weights,shape), 'lsRain', 'srfc')   # pressure levels, hour
vars['lh'] = ([], 12, (inds,weights,shape), 'lh', 'srfc')
vars['sh'] = ([], 12, (inds,weights,shape), 'sh', 'srfc')
vars['lw_net'] = ([], 12, (inds,weights,shape), 'lw_net', 'srfc')
vars['sw_net'] = ([], 12, (inds,weights,shape), 'sw_net', 'srfc')
vars['sw_in'] = ([], 12, (inds,weights,shape), 'sw_in', 'srfc')
vars['t2'] = ([], 12, (inds,weights,shape), 't2', 'srfc')
vars['pblH'] = ([], 12, (inds,weights,shape), 'pblH', 'srfc')
# vars['q2'] = ([], 12, (inds,weights,shape), 'q2', 'srfc')
# vars['shear'] = ([650, 925], 12, (0, 0, 0), 'u_pl', '') # (plinds, plweights, plshape) should use 925 later
# vars['u_mid'] = ([650], 12, (0, 0, 0), 'u_pl', '')
#vars['u_srfc'] = ([925], 12, (0, 0, 0), 'u_pl', '')
# vars['v_mid'] = ([650], 12, (0, 0, 0), 'v_pl', '')
#vars['v_srfc'] = ([925], 12, (0, 0, 0), 'v_pl', '')
# vars['q_mid'] = ([650], 12, (0, 0, 0), 'q_pl', '')  # INPUT IN T * 100!!
# vars['t_mid'] = ([650], 12, (0, 0, 0), 't_pl', '')   # INPUT IN T * 100!!
# vars['t_srfc'] = ([925], 12, (0, 0, 0), 't_pl', '')
# vars['q_srfc'] = ([925], 12, (0, 0, 0), 'q_pl', '')
vars['geoH_srfc'] = ([925], 12, (inds,weights,shape), 'geoH_pl', 'srfc')
vars['geoH_srfc'] = ([750], 12, (inds,weights,shape), 'geoH_pl', 'srfc')
# vars['tcwv'] = ([], 12, (inds,weights,shape), 'colWetMass', 'srfc')
vars['lsRain_noon'] =  ([], 12, (inds,weights,shape), 'lsRain', 'srfc')
#vars['lwout_noon'] =  ([], 12, (inds,weights,shape), 'lw_out_PBLtop', 'srfc')
vars['SM']= ([], 12, (inds, weights, shape), 'SM', 'srfc')
# vars['u_10m_10LT']= ([], 10, (0, 0, 0), 'u10', '')
# vars['v_10m_10LT']= ([], 10, (0, 0, 0), 'v10', '')
# vars['sh_10LT']= ([], 10, (inds, weights, shape), 'sh', 'srfc')
# vars['lh_10LT']= ([], 10, (inds, weights, shape), 'lh', 'srfc')
# vars['SM_10LT']= ([], 10, (inds, weights, shape), 'SM', 'srfc')
###########
vars['v_cross'] = (plevs, 12, (0, 0, 0), 'v_pl', '')
vars['q_cross'] = (plevs, 12, (0, 0, 0), 'q_pl', '')
vars['t_cross'] = (plevs, 12, (0, 0, 0), 't_pl', '')
vars['u_cross'] = (plevs, 12, (0, 0, 0), 'u_pl', '')
vars['omega_cross'] = (plevs, 12, (0, 0, 0), 'omega_pl', '')
vars['v_cross_ST'] = (plevs, h, (0, 0, 0), 'v_pl', '')
vars['q_cross_ST'] = (plevs, h, (0, 0, 0), 'q_pl', '')
vars['t_cross_ST'] = (plevs, h, (0, 0, 0), 't_pl', '')
vars['u_cross_ST'] = (plevs, h, (0, 0, 0), 'u_pl', '')
vars['omega_cross_ST'] = (plevs, h, (0, 0, 0), 'omega_pl', '')

############ ATTENTION - STORM TIME PRESSURE LEVELS ONLY EXIST EVERY 3 HOURS, SO CAN'T SAMPLE FULL DIURNAL CYCLE FOR CROSSSEC

datelist = []
for y,m,d in itertools.product(years, months, days):
    datelist.append(y+m+str(d).zfill(2))

for d in datelist:

    testfiles = glob.glob(out_path + os.sep + d[0:4] + '-' + d[4:6] + '-' + d[6:8] + '_' + str(h) + '*DIR_*.nc')

    if len(testfiles) > 0:
        print(testfiles[0], ' already exists, continue!')
        continue

    print('Doing ', d)
    file_save(data_path, out_path, vars, d, box, tthresh, pos)

