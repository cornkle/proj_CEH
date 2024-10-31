import warnings
warnings.filterwarnings("ignore")
import cartopy.geodesic as cgeo
import pandas as pd
import glob
import xarray as xr
import matplotlib.pyplot as plt
import datetime
import ipdb
from scipy.ndimage.measurements import label
import numpy as np
import multiprocessing




stash_dict={"sm":"STASH_m01s08i223",
            "lhfx":"STASH_m01s03i234",
            "shfx":"STASH_m01s03i217",
            "t2":"STASH_m01s03i236",
            "q2":"STASH_m01s03i237",
            "lw_nsfc":"STASH_m01s02i201",
            "sw_nsfc":"STASH_m01s01i201",
            "olr" : "STASH_m01s02i205",
            "u10" : "STASH_m01s03i225_2",
            "v10" : "STASH_m01s03i226_2",
            "prcp" : "STASH_m01s04i203",
           }



def load_filepath(sim,date,hour,file_list):
    date_str="2006%02d%02d" % (date.month,date.day)
    hour=hour-1 #filenames are 1 hour BEHIND data timestamp
    date_str2=date_str
    date2 = date
    date_prev = date
    date_str_prev = date_str

    if (sim == 'sens') & (hour <= 5):
        date = date - datetime.timedelta(days=1)
        date_str = "2006%02d%02d" % (date.month,date.day)
        date_str2 = date_str
   
    ### this accesses day 2 of sensitivity runs
    if hour>23: # this accesses day 2 of sensitivity runs
        hour = hour-24
        if sim == "control":
            print('Control is continuous simulation, hour has to be <=23')
            return
        if (hour > 5):
            date2 = date + datetime.timedelta(days=1)
            date_str2 = "2006%02d%02d" % (date2.month,date2.day)
            date_str_prev = date_str2
        if (hour <= 5):
            date = date - datetime.timedelta(days=1)
            date_str = "2006%02d%02d" % (date.month,date.day)
            date_prev = date_prev + datetime.timedelta(days=1)
            date_str_prev = "2006%02d%02d" % (date_prev.month,date_prev.day)

    if sim=="sens":
        #root_str="/gws/nopw/j04/lmcs/filtered_soil_moisture_runs/u-da520/20060725T0500Z/Sahel/4km/RA3/um/"
        
        root_str=glob.glob("/gws/nopw/j04/lmcs/sensitivity_runs_wg_mcs/{}T0600Z/*/{}T0600Z/Sahel/1p5km/RA3/um/".format(date_str,date_str2))[0]
    elif sim=="control":
        root_str=glob.glob("/gws/nopw/j04/lmcs/u-cy045_control_run/u-cy045/{}T0000Z/Sahel/1p5km/RA3/um/".format(date_str))[0]
            
    return root_str+"{}_{}_T{:02d}.nc".format(file_list,date_str_prev,hour)

def load_file(ffile,var):
    try:
        ds=xr.open_dataset(ffile)[var]
    except:
        #Deal with funny landseamask file
        ds=xr.open_dataset(ffile,decode_times=False)[var]
    try:
        ds=ds.rename({"grid_latitude_t":"latitude","grid_longitude_t":"longitude"})
    except:
        ipdb.set_trace()
    ds=ds.assign_coords({"longitude":(ds.longitude-360)})

    try:
        ds = ds.rename({"T1HR" : "time"})
    except:
        ipdb.set_trace()
    
    try:
        ds=ds.isel(pressure=slice(None,None,-1))
    except:
        pass

    return ds.squeeze()
    

def loop(ff):
    box = [-15,20,9,20]
    #print(box)
    ds = load_file(ff, stash_dict['prcp'])
    da = ds.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))*3600
    da = da.where(da > 1).mean(['latitude','longitude'])
    #dlist.append(da)
    del ds
    return da


def calc_diurn(flist):

    print('Entering multiprocessing')
    pool = multiprocessing.Pool(processes=5)
    dlist = pool.map(loop, flist)
    pool.close()
    
    combined_dataset = xr.concat(dlist, dim='time')
    combined_dataset['hour'] = combined_dataset['time'].dt.hour
    mean_diurnal_cycle = combined_dataset.groupby('hour').mean()
    return mean_diurnal_cycle

def run():
    ############
    ### Day1 
    flist = []
    for month in [7,8]:
        for day in range(1,32):
            for hour in range(1,25):
    
                if (month == 7) & (day <= 25):
                    continue
                if (month == 9) & (day >=2):
                    if day > 2:
                        continue
                    elif hour >= 5:
                        continue
            
                date = "2006-"+str(month).zfill(2)+"-"+str(day).zfill(2)
                
                try:
                    path = load_filepath("sens",pd.Timestamp(date),hour,"surface_vars")
                except:
                    continue
                flist.append(path)
    print(len(flist))
    day1_flist = sorted(flist)
    
    ############
    ### Day2 
    flist = []
    
    for month in [7,8]:
        for day in range(32):
            for hour in range(1,25):
                date = "2006-"+str(month).zfill(2)+"-"+str(day).zfill(2)

                if (month == 7) & (day <= 25):
                    continue
                if (month == 9) & (day >=2):
                    if day > 2:
                        continue
                    elif hour >= 5:
                        continue
                
                try:
                    path = load_filepath("sens",pd.Timestamp(date),hour+24,"surface_vars")
                except:
                    continue
                flist.append(path)
    print(len(flist))
    day2_flist = sorted(flist)
    
    
    
    ############
    ### Control 
    flist = []
    
    for month in [7,8]:
        for day in range(32):
            for hour in range(1,25):
                date = "2006-"+str(month).zfill(2)+"-"+str(day).zfill(2)
    
                if (month == 7) & (day <= 25):
                    continue
                if (month == 9) & (day >=2):
                    if day > 2:
                        continue
                    elif hour >= 5:
                        continue
                
                try:
                    path = load_filepath("control",pd.Timestamp(date),hour,"surface_vars")
                except:
                    continue
                if path is None:
                    continue
                flist.append(path)
    print(len(flist))
    control_flist = sorted(flist)
    
    d1 = calc_diurn(day1_flist)
    d2 = calc_diurn(day2_flist)
    cont = calc_diurn(control_flist)

    return (d1,d2,cont)

