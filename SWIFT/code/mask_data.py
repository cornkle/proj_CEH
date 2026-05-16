#########################
#Main code for creating a better cloud mask for SEVIRI LST data
#########################
#Based on the gradient of cosecutive images. In the early afternoon ("heating up time")
#a negative gradient suggests unmasked cloud. Later in the day a below expected gradient
#suggests unmasked cloud.
#
#consider each time seperately (there are 96 times of day, 1 every 15min)
#e.g consider 08:00 and nt times each side. Here nt=3
#
#  07:15 07:30 07:45 08:00 08:15 08:30 08:45
#   t-3   t-2   t-1    t    t+1   t+2   t+3
#
#1. Remove all points where any of the images nt*2+1 are missing (cloud mask start and end times often wrong)
#2. Calculate distibution of gradients at each grid point, combining points from 
#   A. All historical years
#   B. All dates within nday days of the date in question. e.g if nday=15 and we are considering the 15th September
#      this would be 29th Aug - 30th September inclusive giving 31 days
#   C. All gradients between consecutive nt. e.g. for example above
#      (t-2) - (t-3), (t-1) - (t-2), (t) - (t-1),(t+1) - (t), (t+2) - (t+1), (t+3) - (t+2)
#3. Save mean and standard deviation of this distribution which will then be used to mask points in real time images
#4. Remove all points with t-(t-1) less than nsig of the standard deviation

#########################
#Create SRA 02/10/2019
#########################
import matplotlib as mpl
mpl.use('Agg')
import glob,os,datetime,h5py,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mask_data_functions as fns

if len(sys.argv)>1:
    time=sys.argv[1]
else:   
   time="0800"

nday=15

########## data to calculate statistics from ############
#fdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward/"  #after first itteration of removal
fdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50/"  #after first itteration of removal
#fdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count100/"  #after first itteration of removal
#fdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50_sig3_count50/"  #after first itteration of removal
#fdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_stats_sig2_count50_mask_sig2_count50/"



########## settings for statistics output: name of source data, directory for stats, directory for plots ############
#source_name="mask_adj_3back1forward"
#outdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward/"
#plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_plots/"
source_name="mask_adj_3back1forward_sig2_count50"
outdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_sig2_count50/"
plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_sig2_count50_plots/"
#source_name="mask_adj_3back1forward_sig2_count100"
#outdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_sig2_count100/"
#plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_sig2_count100_plots/"
#source_name="mask_adj_3back1forward_sig2_count50_sig3_count50"
#outdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_sig2_count50_sig3_count50/"
#plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_sig2_count50_sig3_count50_plots/"
#source_name="mask_adj_3back1forward_stats_sig2_count50_mask_sig2_count50"
#outdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_stats_sig2_count50_mask_sig2_count50/"
#plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_stats_sig2_count50_mask_sig2_count50_plots/"


n_missing=5000
fstart="HDF5_LSASAF_MSG_LST_MSG-Disk_" #start of file names
outstart="HDF5_LSASAF_MSG_LST_MSG-Disk_HistStats_"


#dates: only need for days so times don't matter
start_date=datetime.datetime(2000,5,1,int(time[:2]),int(time[2:]))
end_date=datetime.datetime(2000,10,31,int(time[:2]),int(time[2:]))
#dates: only need for days so times don't matter
dates_all=pd.date_range(start_date,end_date,freq="24h")
use_dates=[str(date).replace("-","").replace("T","").replace(":","")[4:8] for date in dates_all]
time_m1=str(start_date-pd.Timedelta(minutes=15)).replace("-","").replace("T","").replace(":","")[-6:-2]

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
        testf_m1=os.system("h5ls "+fdir+year+date[:2]+"/"+fstart+year+date+time_m1)
        if (testf!=256) & (testf_m1 !=256):
            fid = h5py.File(fdir+year+date[:2]+"/"+fstart+year+date+time , "r")  #assume files are already copied!!!!!
            fid_m1 = h5py.File(fdir+year+date[:2]+"/"+fstart+year+date+time_m1 , "r")  #assume files are already copied!!!!!
            dat_now=fid["LST"][...].astype("float")
            dat_m1=fid_m1["LST"][...].astype("float")
            dat_now[np.where(dat_now==fid["LST"].attrs["MISSING_VALUE"])]=np.nan
            dat_now[np.where(dat_now==fid["LST"].attrs["MISSING_VALUE_ADJ"])]=np.nan
            dat_m1[np.where(dat_m1==fid_m1["LST"].attrs["MISSING_VALUE"])]=np.nan
            dat_m1[np.where(dat_m1==fid_m1["LST"].attrs["MISSING_VALUE_ADJ"])]=np.nan
            if "MISSING_VALUE_HIST" in fid["LST"].attrs.keys():
                dat_now[np.where(dat_now==fid["LST"].attrs["MISSING_VALUE_HIST"])]=np.nan
                dat_m1[np.where(dat_m1==fid_m1["LST"].attrs["MISSING_VALUE_HIST"])]=np.nan
            if "MISSING_VALUE_COUNT" in fid["LST"].attrs.keys():
                dat_now[np.where(dat_now==fid["LST"].attrs["MISSING_VALUE_COUNT"])]=np.nan
                dat_m1[np.where(dat_m1==fid_m1["LST"].attrs["MISSING_VALUE_COUNT"])]=np.nan
            year_dat.append(dat_now-dat_m1)   
            fid.close()
            fid_m1.close()
        else:
            print("Skipping "+fstart+year+date+time)
            year_dat.append(np.zeros((n_missing,n_missing))*np.nan)
    #sort out sizes for missing data
    nrow=np.amin([np.shape(d)[0] for d in year_dat])  #data ny unless all data are missing
    ncol=np.amin([np.shape(d)[1] for d in year_dat])  #data ny unless all data are missing
    year_dat=[dat[:nrow,:ncol] for dat in year_dat]
    data_all.append(np.array(year_dat))
data_all=np.array(data_all)

    
#calculate the mean, standard deviation, and non_missing count for each day
#write each file to fout.

from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages(outdir+"MeanStdev_"+outstart+"_"+time+"-"+time_m1+'.pdf')

#for i_mask in len(dates_all)
fns.output_mean_stdev(data_all,plotdir,outdir,use_dates,dates_all,time,time_m1,nday,nrow,ncol,source_name,outstart,years)

