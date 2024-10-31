import numpy as np
import xarray as xr
import pandas as pd
import os
import glob
import datetime
from collections import OrderedDict
import dask.array as da
from utils import u_interpolate as u_int
from JASMIN import MetUM_variables as mu
import sys
import ipdb


def load_file(ffile):
    try:
        ds = xr.open_dataset(ffile)
    except:
        ds = xr.open_dataset(ffile, decode_times=False)
    try:
        ds = ds.rename({"rlat": "latitude", "rlon": "longitude"})
    except:
        pass
    ds = ds.assign_coords({"longitude": (ds.longitude - 360)})
    try:
        ds = ds.isel(pressure=slice(None, None, -1))
    except:
        pass
    return ds.squeeze()


def olr_to_bt(olr):
    sigma = 5.670373e-8
    tf = (olr / sigma) ** 0.25
    a = 1.228
    b = -1.106e-3
    Tb = (-a + np.sqrt(a ** 2 + 4 * b * tf)) / (2 * b)
    return Tb - 273.15


def filtering(dar, v, outv, pl):
    # if (v == 'q_pl') | (v == 't_pl') | (v == 'lw_out_PBLtop'):
    #     dar.values = np.array(dar.values / 100).astype(float)

    if (v == 'SM'):
        dar = dar.sel(depth=0.05)
        dar = dar.where(dar < 500, other=np.nan)

    if (v == 'lw_out_PBLtop'):
        dar.values = olr_to_bt(dar.values)

    if (v == 'lsRain') | (v == 'totRain'):
        dar.values = dar.values * 3600  # rain to mm/h
        dar.attrs['units'] = 'mm h-1'

    if (v == 't2') | (v == 't_925') | (v == 't_650'):
        dar.values = dar.values - 273.15

    if 'pressure' in dar.coords:
        try:
            dar.values[dar.values == 0] = np.nan
        except ValueError:
            print('Pressure value error!')
            return

        if (len(pl) > 1) & (outv == 'shear'):
            shear = dar.sel(pressure=650).values - dar.sel(pressure=925).values
            dar = dar.sum(dim='pressure').squeeze()

            dar.values = shear
        elif (len(pl) == 1) & (outv != 'shear'):
            dar = dar.sel(pressure=pl[0]).squeeze()
    return dar


def initialize_statistics(shape):
    sum_data = da.zeros(shape, chunks=shape)
    sum_sq_data = da.zeros(shape, chunks=shape)
    count_data = da.zeros(shape, chunks=shape)
    return sum_data, sum_sq_data, count_data


def update_statistics(sum_data, sum_sq_data, count_data, new_data):
    valid_data = ~np.isnan(new_data)
    sum_data += da.where(valid_data, new_data, 0)
    sum_sq_data += da.where(valid_data, new_data ** 2, 0)
    count_data += valid_data.astype(int)
    return sum_data, sum_sq_data, count_data


def save_statistics(statistics, output_file):
    ds = xr.Dataset(statistics)
    ds.to_netcdf(output_file)

def process_composites_by_date(date_group, data_path, output_dir, variables, box, pos, pl_dummy, distx, disty, netcdf=True):

    print('Doing date group', date_group)
    sum_data = {}
    sum_sq_data = {}
    count_data = {}
    sum_anomaly_data = {}
    sum_sq_anomaly_data = {}
    count_anomaly_data = {}

    subdomain_shape = (2 * disty + 1, 2 * distx + 1)

    for variable_name in variables.keys():
        sum_data[variable_name], sum_sq_data[variable_name], count_data[variable_name] = initialize_statistics(subdomain_shape)
        sum_anomaly_data[variable_name], sum_sq_anomaly_data[variable_name], count_anomaly_data[variable_name] = initialize_statistics(subdomain_shape)

    date = date_group['date'].iloc[0]
    print('Doing date', date)
    
    year, month, day = date.split('-')
    file_date = f'{year}{month}{day}'

    for variable_name, (pl, h, interp_params, v, grid) in variables.items():
        print('Doing var', variable_name)
        try:
            filepath = glob.glob(os.path.join(data_path, v, f'*_{file_date}*.nc'))[0]
        except IndexError:
            print(f'No file found for {variable_name} on {date}')
            continue

        try:
            arr = load_file(filepath)
        except OSError:
            print(f'Cannot open file {filepath}')
            ipdb.set_trace()

        orig_v = mu.create_CP4_filename(v)
        
        try:
            dar = arr[orig_v].sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3])).load()
        except:
            ipdb.set_trace()

        if variable_name == 'lsRain_noon':
            dum = dar[dar['time.hour'] == h].squeeze()
            dar = dar[(dar['time.hour'] <= h) & (dar['time.hour'] >= h - 2)].sum('time').squeeze()
            dar = dar.assign_coords(coords={'time': dum.time})
            del dum
        else:
            dar = dar[dar['time.hour'] == h].squeeze()

        dar = filtering(dar, v, variable_name, pl)

        ddate = dar.time.values.item()
        dt = datetime.datetime(ddate.year, ddate.month, ddate.day,
                               ddate.hour, ddate.minute, ddate.second)

        ##### START of local anomaly calculation
        
        if variable_name not in ['lwout_noon', 'lsRain_noon',  'lsRain', 'lw_out_PBLtop']:  #
            if mean_arr == None:
                mean_arr = []
                for dd in range(1, 8):
                    for ndate in [dt - pd.Timedelta(days=dd), dt + pd.Timedelta(days=dd)]:
                        ndatestring = ndate.strftime('%Y%m%d')
                        try:
                            nfile = glob.glob(data_path + os.sep + v + os.sep + f'*_{ndatestring}*.nc')[0]
                        except IndexError:
                            continue
                        
                        try:
                            meanarr = load_file(nfile)
                            mmeanarr = meanarr[orig_v].sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))
                        except OSError:
                            print(f'Could not find clim file for {ndatestring}, continue')
                            continue
                        mmeanarr = mmeanarr[mmeanarr['time.hour'] == h].squeeze()
                        mmeanarr = filtering(mmeanarr, v, variable_name, pl)
                        mean_arr.append(mmeanarr)

                if len(mean_arr) < 10:
                    print('Mean smaller than 10, skipping this case')
                    return
                
                fullmean = xr.concat(mean_arr, dim='time').mean('time')
            anomaly = dar - fullmean
            
        else:
            anomaly = dar

        ##### END of local anomaly calculation

        # regrid to common grid (unstagger wind, bring to landsea mask grid)
        inds, weights, shape = interp_params
        
        if grid == 'srfc':
            try:
                regrid = u_int.interpolate_data(dar.values, inds, weights, shape)
            except ValueError:
                print('Error in regridding')
                continue
        else:
            regrid = dar.values

        da = xr.DataArray(regrid,
                          coords={'time': dar.time, 'latitude': pl_dummy.latitude.values,
                                  'longitude': pl_dummy.longitude.values},
                          dims=['latitude', 'longitude'])
        da.attrs = dar.attrs
        da.values[pos[0], pos[1]] = np.nan  # mask sea

        ##################

        ###### anomaly write
        if grid == 'srfc':
            try:
                regrid_anomaly = u_int.interpolate_data(anomaly.values, inds, weights, shape)
            except ValueError:
                print('Error in regridding')
                continue
        else:
            regrid_anomaly = anomaly.values

        anomaly_da = xr.DataArray(regrid_anomaly,
                                  coords={'time': anomaly.time, 'latitude': pl_dummy.latitude.values,
                                          'longitude': pl_dummy.longitude.values},
                                  dims=['latitude', 'longitude'])
        anomaly_da.attrs = anomaly.attrs
        anomaly_da.values[pos[0], pos[1]] = np.nan  # mask sea

        #########################
            
        mean_arr = None
        
        for idx, row in date_group.iterrows():
            lon = row['tminlon']
            lat = row['tminlat']
            
            # Center subdomain on minlon and minlat
            point = da.sel(latitude=lat, longitude=lon, method='nearest')

            xpos = np.where(da['longitude'].values == point['longitude'].values)[0][0]
            ypos = np.where(da['latitude'].values == point['latitude'].values)[0][0]

            try:
                subdomain = da.isel(latitude=slice(ypos - disty, ypos + disty + 1),
                                    longitude=slice(xpos - distx, xpos + distx + 1))
            except IndexError:
                print('IndexError: subdomain out of bounds')
                continue

            if (len(subdomain.latitude) != subdomain_shape[0]) | (len(subdomain.longitude) != subdomain_shape[1]):
                print('Subdomain shape mismatch, skipping this case')
                continue

            sum_data[variable_name], sum_sq_data[variable_name], count_data[variable_name] = update_statistics(
                sum_data[variable_name], sum_sq_data[variable_name], count_data[variable_name], subdomain.values
            )

            ####################
            try:
                subdomain_anomaly = anomaly_da.isel(latitude=slice(ypos- disty, ypos + disty + 1),
                                                    longitude=slice(xpos - distx, xpos + distx + 1))
            except IndexError:
                print('IndexError: anomaly subdomain out of bounds')
                continue

            if (len(subdomain_anomaly.latitude) != subdomain_shape[0]) | (len(subdomain_anomaly.longitude) != subdomain_shape[1]):
                print('Anomaly subdomain shape mismatch, skipping this case')
                continue

            sum_anomaly_data[variable_name], sum_sq_anomaly_data[variable_name], count_anomaly_data[variable_name] = update_statistics(
                sum_anomaly_data[variable_name], sum_sq_anomaly_data[variable_name], count_anomaly_data[variable_name], subdomain_anomaly.values
            )
            
            ##### finish anomaly write

            #### netcdf writeout for testing
            if netcdf==True:
                for dat in zip((subdomain, subdomain_anomaly), ('mean', 'anom')):
                    savefile = output_dir + os.sep + date + '_lonXlat_'+str(np.round(lon,1))+'_'+str(np.round(lat,1))+'_' + dat[1]+'.nc'

                    filt = dat[0].copy(deep=True)
                    filt = filt.assign_coords(
                        {'longitude': np.arange(distx * -1, distx + 1), 'latitude': np.arange(distx * -1, distx + 1)})
                    
                    #filt = filt.to_dataset()
                    filt.to_netcdf(savefile)
                    del filt
            
            ##### finish netcdf writeout 
           
        del arr

    return sum_data, sum_sq_data, count_data, sum_anomaly_data, sum_sq_anomaly_data, count_anomaly_data


def process_composites(location_table, data_path, output_dir, variables, box, pos, pl_dummy, distx, disty):
    sum_data_all = {}
    sum_sq_data_all = {}
    count_data_all = {}
    sum_anomaly_data_all = {}
    sum_sq_anomaly_data_all = {}
    count_anomaly_data_all = {}
    
    for variable_name in variables.keys():
        sum_data_all[variable_name], sum_sq_data_all[variable_name], count_data_all[variable_name] = initialize_statistics((2 * disty + 1, 2 * distx + 1))
        sum_anomaly_data_all[variable_name], sum_sq_anomaly_data_all[variable_name], count_anomaly_data_all[variable_name] = initialize_statistics((2 * disty + 1, 2 * distx + 1))
    
    for date, date_group in location_table.groupby('date'):
        sum_data, sum_sq_data, count_data, sum_anomaly_data, sum_sq_anomaly_data, count_anomaly_data = process_composites_by_date(date_group, data_path, output_dir, variables, box, pos, pl_dummy, distx, disty)
        
        for variable_name in variables.keys():
            sum_data_all[variable_name] += sum_data[variable_name]
            sum_sq_data_all[variable_name] += sum_sq_data[variable_name]
            count_data_all[variable_name] += count_data[variable_name]
            sum_anomaly_data_all[variable_name] += sum_anomaly_data[variable_name]
            sum_sq_anomaly_data_all[variable_name] += sum_sq_anomaly_data[variable_name]
            count_anomaly_data_all[variable_name] += count_anomaly_data[variable_name]

    statistics = {}
    statistics_anomaly = {}
    
    for variable_name in variables.keys():
        statistics[f'{variable_name}_sum'] = (['latitude', 'longitude'], sum_data_all[variable_name].compute())
        statistics[f'{variable_name}_sum_sq'] = (['latitude', 'longitude'], sum_sq_data_all[variable_name].compute())
        statistics[f'{variable_name}_count'] = (['latitude', 'longitude'], count_data_all[variable_name].compute())
        
        statistics_anomaly[f'{variable_name}_sum'] = (['latitude', 'longitude'], sum_anomaly_data_all[variable_name].compute())
        statistics_anomaly[f'{variable_name}_sum_sq'] = (['latitude', 'longitude'], sum_sq_anomaly_data_all[variable_name].compute())
        statistics_anomaly[f'{variable_name}_count'] = (['latitude', 'longitude'], count_anomaly_data_all[variable_name].compute())

    save_statistics(statistics, os.path.join(output_dir, 'intermediate_statistics_10LT.nc'))
    save_statistics(statistics_anomaly, os.path.join(output_dir, 'intermediate_statistics_anomaly_10LT.nc'))


#######################################################################################################
# Main script execution
### Inputs:
##### Provide hour when you run script!!! as sys.argv[1]

box = [-18, 25, 9, 27]  # W- E , S - N geographical coordinates box
distx = 57  # 57 = 250 km at 4.4 res, 500km across
disty = 57

# YEARS = np.array(np.arange(1998,2007), dtype=str)
# days = np.array(np.arange(1,32), dtype=str)

tthresh = -50 # chosen temperature threshold, e.g. -50, -60, -70
h= int(sys.argv[1]) # hour to extract
m= int(sys.argv[2]) # month to extract

fdir = str(sys.argv[3]) # 'hist' or 'future'
if fdir == 'hist':
    ftag = 'historical'
else:
    ftag = 'future'

main = '/home/users/cornkle/linked_CP4/'
main_lmcs = '/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/'
data_path = main + '/'+fdir
ancils_path = '/home/users/cornkle/impala/shared/CP4A/ncfiles/4km/ANCILS/'
out_path = main_lmcs + 'CP4_box_JASMIN/CP4_'+ftag+'_5000km2_-50_box_v3'


location_table_files = glob.glob('/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/CP4_MCS_table_runOnJasmin/*'+str(h).zfill(2)+'h.csv')
# Read and merge all location tables
location_tables = [pd.read_csv(file) for file in location_table_files]
merged_location_table = pd.concat(location_tables, ignore_index=True)
filtered_location_table = merged_location_table[(merged_location_table['hour'] == h) & (merged_location_table['month'] == m)]


plglob = glob.glob(data_path + '/q_pl/*.nc')
pl_dummy = load_file(plglob[0])[mu.create_CP4_filename('q_pl')]

srfcglob = glob.glob(data_path + '/lw_out_PBLtop/*.nc')
srfc_dummy =load_file(srfcglob[0])[mu.create_CP4_filename('lw_out_PBLtop')]

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
lons, lats = np.meshgrid(pl_dummy.longitude.values, pl_dummy.latitude.values)
inds, weights, shape = u_int.interpolation_weights(srfc_dummy.longitude, srfc_dummy.latitude, pl_dummy.longitude, pl_dummy.latitude)

# Define your variables here
variables = OrderedDict({
    'lw_out_PBLtop': ([], 12, (inds, weights, shape), 'lw_out_PBLtop', 'srfc'),
    'lsRain': ([], 12, (inds, weights, shape), 'lsRain', 'srfc'),
    'shear': ([650, 925], 12, (0, 0, 0), 'u_pl', ''),
    # 'u_mid': ([650], 12, (0, 0, 0), 'u_pl', ''),
    # 'u_srfc': ([925], 12, (0, 0, 0), 'u_pl', ''),
    # 'v_mid': ([650], 12, (0, 0, 0), 'v_pl', ''),
    # 'v_srfc': ([925], 12, (0, 0, 0), 'v_pl', ''),
    # 'q_mid': ([650], 12, (0, 0, 0), 'q_pl', ''),
    # 't_mid': ([650], 12, (0, 0, 0), 't_pl', ''),
    # 't_srfc': ([925], 12, (0, 0, 0), 't_pl', ''),
    # 'q_srfc': ([925], 12, (0, 0, 0), 'q_pl', ''),
    # 'geoH_srfc': ([925], 12, (inds, weights, shape), 'geoH_pl', 'srfc'),
    # 'tcwv': ([], 12, (inds, weights, shape), 'tcwv', 'srfc'),
    # 'sh': ([], 12, (inds, weights, shape), 'sh', 'srfc'),
    # 'lh': ([], 12, (inds, weights, shape), 'lh', 'srfc'),
    # 't2': ([], 12, (inds, weights, shape), 't2', 'srfc'),
    # 'q2': ([], 12, (inds, weights, shape), 'q2', 'srfc'),
    # 'lsRain_noon': ([], 12, (inds, weights, shape), 'lsRain', 'srfc'),
    # 'lwout_noon': ([], 12, (inds, weights, shape), 'lw_out_PBLtop', 'srfc'),
    # 'SM': ([], 12, (inds, weights, shape), 'SM', 'srfc'),
    # 'u_10m_10LT': ([], 10, (0, 0, 0), 'u10', ''),
    # 'v_10m_10LT': ([], 10, (0, 0, 0), 'v10', ''),
    # 'sh_10LT': ([], 10, (inds, weights, shape), 'sh', 'srfc'),
    # 'lh_10LT': ([], 10, (inds, weights, shape), 'lh', 'srfc'),
    # 'SM_10LT': ([], 10, (inds, weights, shape), 'SM', 'srfc')
})

process_composites(filtered_location_table[0:10], data_path, out_path, variables, box, pos, pl_dummy, distx, disty)

