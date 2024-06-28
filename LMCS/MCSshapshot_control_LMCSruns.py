from ben_notebook_functions import *
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
from MCS_snapshots import MCS_table_main_model_LMCSruns as MCS_table_main
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
            "prcp" : "STASH_m01s04i203"
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


def olr_to_bt(olr):
    #Application of Stefan-Boltzmann law
    sigma = 5.670373e-8
    tf = (olr/sigma)**0.25
    #Convert from bb to empirical BT (degC) - Yang and Slingo, 2001
    a = 1.228
    b = -1.106e-3
    Tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return Tb - 273.15


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
        ds = ds.rename({"T1HR_MN_rad_diag" : "time"})
    except:
        try:
            ds = ds.rename({"T1HR" : "time"})
        except:
            ipdb.set_trace()
    
    try:
        ds=ds.isel(pressure=slice(None,None,-1))
    except:
        pass

    return ds.squeeze()


def dictionary():

    dic = {}
    vars = ['date', 'month', 'hour', 'minute', 'year', 'day', 'area', '70area', '65area','60area', 'tmin',
            'minlon', 'minlat', 'maxlon', 'maxlat', 'clon', 'clat', 'tminlon', 'tminlat',
            'tmin', 'tmean', 'tp1','stormID', 'pmax', 'pmean', 'p99']#, 'cloudMask', 'tir']


    for v in vars:
        dic[v] = []
    return dic



def mcs_define(array, thresh, min_area=None, max_area=None, minmax_area=None):
    """

    :param array: 2d input array
    :param thresh: cloud threshold
    :param min_area: minimum area of the cloud
    :param max_area: maximum area of the cloud
    :param minmax_area: tuple indicating only clouds bigger than tuple[0] and smaller than tuple[1]
    :return: 2d array with labelled blobs
    """
    array[array >= thresh] = 0  # T threshold maskout
    array[np.isnan(array)] = 0  # set ocean nans to 0

    labels, numL = label(array)

    u, inv = np.unique(labels, return_inverse=True)
    n = np.bincount(inv)

    goodinds = u[u!=0]

    if min_area != None:
        goodinds = u[(n>=min_area) & (u!=0)]
        badinds = u[n<min_area]

        # for b in badinds:
        #     pos = np.where(labels==b)
        #     labels[pos]=0

    if max_area != None:
        goodinds = u[(n<=max_area)  & (u!=0)]
        badinds = u[n>max_area]

    if minmax_area != None:
        goodinds = u[(n <= minmax_area[1]) & (u != 0) & (n>=minmax_area[0])]
        badinds = u[(n > minmax_area[1]) | (n < minmax_area[0])]

    if (min_area is not None) | (max_area is not None) | (minmax_area is not None):
        for b in badinds:
            pos = np.where(labels==b)
            labels[pos]=0

    return labels, goodinds


def process_tir_image(ff, data_res=1.5, t_thresh=-50, min_mcs_size=25):
    """
    This function cuts out MCSs. By default, an MCS is defined as contiguous brightness temperature area at <=-50 degC over >= 5000km2.
    :param ctt: brightness temperature image (in degC)
    :param data_res: spatial resolution of input image (approximately, in km - this defines how many pixel are needed to define an MCS)
    :param t_thresh: temperature threshold (in degC) for considered contiguous cloud / MCS check.
    :param min_mcs_size: minimum size of contiguous cloud to be considered an MCS (in km2)
    :return: dictionary with dates and MCS characteristics from cloud top information only.
    """
    print('Doing', ff)
    
    box = [-17,23,6,20]
    orig_var = stash_dict['olr']
    orig_var2 = stash_dict['prcp']
    
    try:
        da = load_file(ff, orig_var).sel(latitude=slice(box[2],box[3]), longitude=slice(box[0],box[1]))
        da2 = load_file(ff, orig_var2).sel(latitude=slice(box[2],box[3]), longitude=slice(box[0],box[1]))
    except:
        print('File not found', ff)
        return
    ctt = olr_to_bt(da.squeeze()) # convert olr to bt
    
    min_pix_nb = min_mcs_size / data_res**2

    max_pix_nb = 500000 / data_res**2  # this is to capture satellite artefacts that come in large contiguous stripes.

    labels, goodinds = mcs_define(ctt.values, t_thresh, minmax_area=[min_pix_nb, max_pix_nb]) # 7.7x7.7km = 64km2 per pix in gridsat? 83 pix is 5000km2
    dic = dictionary()

    for g in goodinds:

        if g==0:
            continue

        pos = np.where(labels==g)
       
        if (np.sum(pos[0]) == 0) | (np.sum(pos[1]) == 0):
            print('COORDINATES ALL 0, PLEASE CHECK')
            continue
        npos = np.where(labels!=g)
        datestr = str(int(ctt['time.year']))+'-'+str(int(ctt['time.month'])).zfill(2)+'-'+str(int(ctt['time.day'])).zfill(2)+'_'+\
                      str(int(ctt['time.hour'])).zfill(2)+':'+str(int(ctt['time.minute'])).zfill(2)

        dic['date'].append(datestr)
        dic['month'].append(int(ctt['time.month']))
        dic['hour'].append(int(ctt['time.hour']))
        dic['year'].append(int(ctt['time.year']))
        dic['day'].append(int(ctt['time.day']))
        dic['minute'].append(int(ctt['time.minute']))

        storm = ctt.copy()
        storm2 = da2.copy()*3600
        
        storm.values[npos] = np.nan
        storm2.values[npos] = np.nan
        
        tmin_pos = np.nanargmin(storm.values)
        tpos_2d = np.unravel_index(tmin_pos, storm.shape)

        latmin = np.nanmin(ctt.latitude.values[pos[0]])
        latmax = np.nanmax(ctt.latitude.values[pos[0]])
        lonmin = np.nanmin(ctt.longitude.values[pos[1]])
        lonmax = np.nanmax(ctt.longitude.values[pos[1]])

        dic['area'].append(np.sum(np.isfinite(storm.values)))
        dic['70area'].append(np.sum(storm.values<=-70))
        dic['60area'].append(np.sum(storm.values<=-60))
        dic['65area'].append(np.sum(storm.values<=-65))
        dic['minlon'].append(np.round(lonmin,3))
        dic['minlat'].append(np.round(latmin,3))
        dic['maxlon'].append(np.round(lonmax,3))
        dic['maxlat'].append(np.round(latmax,3))
        dic['clon'].append(np.round(lonmin + (lonmax - lonmin)/2,3))
        dic['clat'].append(np.round(latmin + (latmax - latmin)/2,3))
        dic['tmin'].append(np.round(np.nanmin(storm),2))
        dic['tminlat'].append(np.round(ctt.latitude[tpos_2d[0]].values,3))
        dic['tminlon'].append(np.round(ctt.longitude[tpos_2d[1]].values,3))
        dic['tmean'].append(np.round(np.nanmean(storm),2))
        dic['tp1'].append(np.round(np.nanpercentile(storm, 1),2))
        dic['stormID'].append(datestr + '_' + str(g))
        dic['pmax'].append(np.round(np.nanmax(storm2.values),2))
        dic['pmean'].append(np.round(np.nanmean(storm2.values),2))
        dic['p99'].append(np.round(np.nanpercentile(storm2.values,99),2))
        # dic['cloudMask'].append(labels==g)
        # dic['tir'].append(storm.values)

        del storm

    # for k in dic.keys():
    #     print(k, len(dic[k]))
    del labels
    del ctt
    del da

    return pd.DataFrame.from_dict(dic)


#def make_table(infiles):

    #all_tab = []
    # for ff in infiles[0:10]:
    #     print('Doing file', ff)
    #     all_tab.append(process_tir_image(da, 1.5))
    # return pd.concat(all_tab)



############
### Control 
flist = []

for month in [7,8]:
    for day in range(32):
        for hour in range(1,25):
            date = "2006-"+str(month).zfill(2)+"-"+str(day).zfill(2)

            if (month == 7) & (day == 25):
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
flist = sorted(flist)

pool = multiprocessing.Pool(processes=5)
all_tab = pool.map(process_tir_image, flist)
pool.close()

# all_tab = []
# for ff in flist[0:480]:
#     all_tab.append(process_tir_image(ff))

tab = pd.concat(all_tab)

# tab.pop('cloudMask')
# tab.pop('tir')

pd.DataFrame(tab).to_csv('SM_runs_control_snapshots.csv')



