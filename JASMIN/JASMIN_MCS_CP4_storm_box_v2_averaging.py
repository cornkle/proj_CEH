import glob
import ipdb
#import cdo
import os
import xarray as xr
import sys


# Function to compute mean from a list of netCDF files
def compute_mean(files):
    datasets = [xr.open_dataset(file) for file in files]
    combined = xr.concat(datasets, dim='cases')  # Assuming time is a common dimension
    #ipdb.set_trace()
    mean = combined.mean(dim='cases')
    return mean
    

### Inputs:
##### Provide hour when you run script!!! as sys.argv[1]
HOURS = [15,17,18]

FTAG = ['historical', 'future']

ATAG = ['anom', 'mean']

SHAPE = ['', '_DIR']
XY = ['XDIR', 'YDIR']

main_lmcs = '/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/'


MONTHS = ([7,8,9])

for hh in HOURS:
    print('Doing hour', hh)
    for ff in FTAG:
        for aa in ATAG:
            for ss in SHAPE:
                in_path = main_lmcs + 'CP4_box_'+aa+'_JASMIN/CP4_'+ff+'_5000km2_-50_box_'+aa+'_v3'+ss+'/'
                for mm in MONTHS:
                    for xx in XY:
                        
                        savefile = main_lmcs + aa + '_' + str(mm).zfill(2) +'_'+str(hh).zfill(2)+'h'+'_'+xx+'.nc'

                        if os.path.isfile(savefile):
                                print('File exists, continue')
                                continue


                        files = sorted(glob.glob(in_path + os.sep + '*-' + str(mm).zfill(2) + '-*_' + str(hh).zfill(2)+':*' + '_'+xx+'_*.nc'))
                         
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
                            continue
                            


