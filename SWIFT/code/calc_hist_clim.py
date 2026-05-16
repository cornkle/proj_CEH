#########################
#Main code for calculating and outputing
#A: Historical climitology informaion
#1. mean of final masked data at each time with 15 days either side and all historical years
#2. Standard deviation of data at each time with 15 days either side and all historical years
#B: Anomoly from the historiccal mean
#
#########################

#########################
#Create SRA 02/10/2019
#########################
import matplotlib as mpl
mpl.use('Agg')
import os
import datetime
import h5py
import sys
import numpy as np
import pandas as pd
import mask_data_functions as fns

if len(sys.argv)>1:
    time=sys.argv[1]
else:   
   time="0800"

nday=15
#fdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj/"
#source_name="itmasked_adj1"
fdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50/"  #after first itteration of removal
source_name="adj_3back1forward_sig2_count50"
outdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_clim/"
plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_clim_plots/"

n_missing=5000
fstart="HDF5_LSASAF_MSG_LST_MSG-Disk_" #start of file names
outstart="HDF5_LSASAF_MSG_LST_MSG-Disk_HistStats_"


#dates: only need for days so times don't matter
start_date=datetime.datetime(2000,5,1,int(time[:2]),int(time[2:]))
end_date=datetime.datetime(2000,10,31,int(time[:2]),int(time[2:]))
dates_all=pd.date_range(start_date,end_date,freq="24h")
use_dates=[str(date).replace("-","").replace("T","").replace(":","")[4:8] for date in dates_all]

years=[str(year) for year in range(2004,2016,1)]

#read in all the data for this time. 
#loop over dates. Statistics are produced for each date
data_all=[]
for date in use_dates:
    print(date)
    #loop over all years
    year_dat=[]
    for year in years:    
        #print("Processing year:",year)        
        testf=os.system("h5ls "+fdir+year+date[:2]+"/"+fstart+year+date+time)
        if (testf!=256):
            fid = h5py.File(fdir+year+date[:2]+"/"+fstart+year+date+time , "r")  #assume files are already copied!!!!!
            dat_now=fid["LST"][...].astype("float")
            dat_now[np.where(dat_now==fid["LST"].attrs["MISSING_VALUE"])]=np.nan
            dat_now[np.where(dat_now==fid["LST"].attrs["MISSING_VALUE_ADJ"])]=np.nan
            if "MISSING_VALUE_HIST" in fid["LST"].attrs.keys():
                dat_now[np.where(dat_now==fid["LST"].attrs["MISSING_VALUE_HIST"])]=np.nan
            if "MISSING_VALUE_COUNT" in fid["LST"].attrs.keys():
                dat_now[np.where(dat_now==fid["LST"].attrs["MISSING_VALUE_COUNT"])]=np.nan
            year_dat.append(dat_now)           
            fid.close()
        else:
            print("Skipping "+fstart+year+date+time)
            year_dat.append(np.zeros((n_missing,n_missing))*np.nan)
    #sort out sizes for missing data
    nrow=np.amin([np.shape(d)[0] for d in year_dat])  #data ny unless all data are missing
    ncol=np.amin([np.shape(d)[1] for d in year_dat])  #data ny unless all data are missing
    year_dat=[dat[:nrow,:ncol] for dat in year_dat]
    data_all.append(np.array(year_dat))

data_all=np.array(data_all) #all the masked data for the year in question

#calculate the mean, standard deviation
#write each file to fout.

#pp = PdfPages(outdir+"MeanStdev_"+outstart+"_"+time+"-"+time_m1+'.pdf')

#for i_mask in len(dates_all)
fns.output_mean_stdev(data_all,plotdir,outdir,use_dates,dates_all,time,time,nday,nrow,ncol,source_name,outstart,years,False)

