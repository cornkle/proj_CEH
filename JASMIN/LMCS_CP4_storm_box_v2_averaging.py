import glob
import ipdb
#import cdo
import os
import xarray as xr
import sys
import multiprocessing


# Function to compute mean from a list of netCDF files
def compute_mean(files):
    datasets = [xr.open_dataset(file) for file in files]
    combined = xr.concat(datasets, dim='cases')  # Assuming time is a common dimension
    #ipdb.set_trace()
    mean = combined.mean(dim='cases')
    return mean


def loop(hh):
    
    FTAG = ['fut','hist'] #'hist', 
    ATAG = ['anom'] #, 'mean'
    MONTHS = (['Jul-Sep'])  #7,8,9, 'Jul-Sep'
    SHAPE = ['pl_', ''] #, 'pl_'
    XY = ['XDIR', 'YDIR']  #'XDIR', 'YDIR', 
    
    main_lmcs = '/gws/nopw/j04/lmcs/cklein/CP_models/MCS_files/WAf/CP4_box_JASMIN/mean3h_v2/'

    print('Doing hour', hh)
    
    for ff in FTAG:
        for aa in ATAG:
            for ss in SHAPE:
                in_path = main_lmcs +ss+aa+'_'+ff+'/'
                for mm in MONTHS:
                    print('Doing month', mm)
                    for xx in XY:
                        if 'DIR' in xx:
                            #ipdb.set_trace()
                            savefile = main_lmcs + ff + '_' + aa + '_' + str(mm).zfill(2) +'_'+str(hh).zfill(2)+'h'+'_'+xx+'.nc'
                        else:
                            savefile = main_lmcs + ff + '_' + aa + '_' + str(mm).zfill(2) +'_'+str(hh).zfill(2)+'h'+'_2d.nc'

                        if os.path.isfile(savefile):
                                print('File exists, continue')
                                continue

                        if (mm == 'Jul-Sep'):
                            files = sorted(glob.glob(in_path + os.sep + '*-' + '*' + '-*_' + str(hh).zfill(2)+':*' + '_'+xx+'*.nc'))
                            #ipdb.set_trace()
                        else:
                            files = sorted(glob.glob(in_path + os.sep + '*-' + str(mm).zfill(2) + '-*_' + str(hh).zfill(2)+':*' + '_'+xx+'*.nc'))
                            
                        #ipdb.set_trace()
                        if len(files) > 0:
                            if '2006' not in files[-1]:
                                print(files[-1], '2006 is not in, continue')
                                #ipdb.set_trace()
                                continue
                            print('File len', len(files))
                            mean = compute_mean(files)
                
                            # try:
                            #     os.remove(savefile)
                            # except OSError:
                            #     pass
                    
                            comp = dict(zlib=True, complevel=5)
                            encoding = {var: comp for var in mean.data_vars}
                            mean.to_netcdf(path=savefile, mode='w', encoding=encoding)
                            print('Saved ' + savefile)

                        else:
                            return

    
        

### Inputs:
HOURS = [19]
print('Doing', HOURS)

pool = multiprocessing.Pool(processes=3)
res = pool.map(loop, HOURS)
pool.close()

# for h in HOURS:
#     loop(h)
                            


