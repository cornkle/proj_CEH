import glob
import ipdb
import cdo
import os
import xarray as xr

# Function to perform t-test for significance
def t_test(anomaly_files, mean_of_means):
    anomalies = [xr.open_dataset(file).to_array() for file in anomaly_files]
    anomaly_values = np.array([anomaly.values for anomaly in anomalies])
    mean_values = mean_of_means.to_array().values
    
    # Perform t-test along the desired dimension (e.g., time or spatial dimensions)
    t_stat, p_values = stats.ttest_1samp(anomaly_values, mean_values, axis=0)
    
    # Convert p_values to xarray DataArray
    p_values_xr = xr.DataArray(p_values, dims=anomalies[0].dims, coords=anomalies[0].coords)
    return p_values_xr

### Inputs:
##### Provide hour when you run script!!! as sys.argv[1]
HOURS = [15,17,18]

FTAG = ['historical', 'future']

ATAG = ['anom', 'mean']

SHAPE = ['', 'DIR']
XY = ['XDIR', 'YDIR']

main_lmcs = '/home/users/cornkle/lmcs/cklein/CP_models/MCS_files/'

MONTHS = ([7,8,9])

for hh in HOURS:
    for ff in FTAG:
        for aa in ATAG:
            for ss in SHAPE:
                in_path = main_lmcs + 'CP4_box_'+aa+'_JASMIN/CP4_'+ff+'_5000km2_-50_box_'+aa+'_v3'+ss'/'
                for mm in MONTHS:
                    for xx in XY:

                        files = glob.glob(out_path + os.sep + '*-' + str(mm).zfill(2) + '-*_' + str(hh).zfill(2)+':*' + '_'+xx+'_*.nc')
                        if len(files) > 0:
                            cdo.ensmean
                            
                            


