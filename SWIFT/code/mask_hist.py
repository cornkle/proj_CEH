#########################
#Main code for masking SEVIRI_LST data based on historicl distribution of points
#historical mean and stdev grids must be make first by running mask_data.py
#
#One job for each historical year. Loop over dates in that year
#1. Read in the origional non-masked grids. Diff the grids to get arrays of (t)-(t-1) and (t+1)-t (read only)
#2. Read in the mean and stdev grids for this day(read only)
#3. Read in the masked_adj grids for this day (read-write)
#4. Modify the masked_adj grids 
#	- mask points with count<50 with MASKED_COUNT=-6000
#	- mask points where value<mean - a*stdev with MASKED_HIST1=-5000
#	- add attributes for these masked values

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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

if len(sys.argv)>1:
    year=int(sys.argv[1])               #one job per year (loops over dates)
else:   
    year=2012

########data to be masked################
#fdir_out="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50/"
#fdir_out="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count100/"
#fdir_out="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_stats_sig2_count50_mask_sig2_count50/"
fdir_out="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50/"
#################the raw data################
fdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_sorted/"

#################what mask to use################
#fdir_hist="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward_sig2_count50/"
fdir_hist="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_stats_adj_3back1forward/"

################plots of the masked data################
plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50_plots/"
#plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count100_plots/"
#plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50_sig3_count50_plots/"
#plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50_plots_count100_plots/"
#plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_stats_sig2_count50_mask_sig2_count50_plots/"

count_thresh=50
#count_thresh=100
missing_count=-6000
missing_hist=-5000
asig=2


fstart="HDF5_LSASAF_MSG_LST_MSG-Disk_" #start of file names
fstart_stats="HDF5_LSASAF_MSG_LST_MSG-Disk_HistStats_"
#logfile to output the number of missing points
logfile=open("/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/logs/Mask_data_count_hist_"+str(year)+".log","w+")
n_missing=5000

#dates: only need for days so times don't matter
start_date=datetime.datetime(year,5,1,8,0)
end_date=datetime.datetime(year,10,31,8,0)
#end_date=datetime.datetime(year,7,31,8,0) #TESTING!
#end_date=datetime.datetime(year,5,1,8,0) #TESTING!

dates_all=pd.date_range(start_date,end_date,freq="24h")
use_dates=[str(date).replace("-","").replace("T","").replace(":","")[:8] for date in dates_all]

#loop over all dates in year to be processed  

for use_date in use_dates:  #20120625, 20120626,20120627  [56:58]

    pp = PdfPages(plotdir+"Masked_LST_"+use_date+'.pdf')
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams.update({'font.size': 8})

    print(use_date)   
    logfile.write(use_date+",")
    #get times of day and files#l to load
    start_time=datetime.datetime(int(use_date[:4]),int(use_date[4:6]),int(use_date[6:]),7,0)
    end_time=datetime.datetime(int(use_date[:4]),int(use_date[4:6]),int(use_date[6:]),19,45)
    #end_time=datetime.datetime(int(use_date[:4]),int(use_date[4:6]),int(use_date[6:]),8,45)  #TESTING!!!

    dates_all=pd.date_range(start_time,end_time,freq="15min")
    toload=[date.replace("-","").replace("T","").replace(":","")[:6]+"/"+\
            fstart+date.replace("-","").replace("T","").replace(":","")[:12] for date in np.datetime_as_string(dates_all)]
    toload_stats=[np.datetime_as_string(np.datetime64(date)).replace("-","").replace("T","").replace(":","")[4:6]+"/"+\
                  fstart_stats+np.datetime_as_string(np.datetime64(date)).replace("-","").replace("T","").replace(":","")[4:8]+"_"+\
                  np.datetime_as_string(np.datetime64(date)).replace("-","").replace("T","").replace(":","")[8:12]+"-"+\
                  np.datetime_as_string(np.datetime64(date-pd.Timedelta(minutes=15))).replace("-","").replace("T","").replace(":","")[8:12]\
                  for date in dates_all[1:]]

    #read in 1 day of data
    data_all=[] #this holds the file handles (these will be modified)
    data_all_adj=[]
    LST_data=[] #this holds the actual LST data
    LST_data_adj=[] #this holds the actual LST data
    data_hist_all=[]

    for i_f,f in enumerate(toload):

        if i_f>0:   #there is no historical stats file for the earliest time in the day as this is t-[t-1] product
            fid_hist= h5py.File(fdir_hist+toload_stats[i_f-1], "r")
            data_hist_all.append(fid_hist)

        testf=os.system("h5ls "+fdir+f)
        if testf!=256:
            fid = h5py.File(fdir+f, "r")  #assume files are already copied!!!!!
            fid_adj = h5py.File(fdir_out+f, "r+")  #assume files are already copied!!!

            if "MISSING_VALUE_ADJ" in fid["LST"].attrs.keys():
                sys.exit("ERROR! this data has already been masked:please check that the specified fdir is correct!")
            else: 
                #keep the fids for access of data
                data_all.append(fid)
                data_all_adj.append(fid_adj)
                #for LST data keep a copy of the data for safety (need to reshape/modify)
                LST_data.append(fid["LST"])
                LST_data_adj.append(fid_adj["LST"])
        else:
            data_all.append(np.nan)  #don't want to mask when we have a full missing image in the image range
            data_all_adj.append(np.nan)  #don't want to mask when we have a full missing image in the image range
            LST_data.append(np.zeros((n_missing,n_missing)))
            LST_data_adj.append(np.zeros((n_missing,n_missing)))

          
    #sort out sizes for missing data
    nx=np.amin([np.shape(d)[0] for d in LST_data])  #data ny unless all data are missing
    ny=np.amin([np.shape(d)[1] for d in LST_data])  #data ny unless all data are missing

    if (nx==n_missing)& (ny==n_missing):
        print("all missing data for ",date)
    else:
        LST_data=np.array([data[:nx,:ny] for data in LST_data]) #create numpy array
        LST_data_adj=np.array([data[:nx,:ny] for data in LST_data_adj]) #create numpy array
        #create difference arrays from LST_data (unmasked version). Do once for the full day
        LST_diff=np.diff(LST_data,axis=0)

        for i_mask,outmask in enumerate(toload):            #loop over times of day and mask 
            if not type(data_all_adj[i_mask])==np.float:    #don't try to mask a missing image
                already_missing=(data_all_adj[i_mask]["LST"][...]==data_all_adj[i_mask]["LST"].attrs["MISSING_VALUE"]) |\
                                (data_all_adj[i_mask]["LST"][...]==data_all_adj[i_mask]["LST"].attrs["MISSING_VALUE_ADJ"])

	            #1. Mask the ooints where count<count_thresh
                if i_mask==0:  #can only look at forward difference for first image of the day
                    for key in data_all_adj[i_mask].keys():
                        data_all_adj[i_mask][key][(data_hist_all[i_mask]["count"][...]<count_thresh)& (already_missing==False)]=missing_count
                        data_all_adj[i_mask][key].attrs.create("MISSING_VALUE_COUNT",missing_count)
                    already_missing[np.where(data_hist_all[i_mask]["count"][...]<count_thresh)]=True

                elif i_mask==(len(toload)-1): #can only look at backward difference for last image of the day 
                    for key in data_all_adj[i_mask].keys():
                        data_all_adj[i_mask][key][(data_hist_all[i_mask-1]["count"][...]<count_thresh)& (already_missing==False)]=missing_count
                        data_all_adj[i_mask][key].attrs.create("MISSING_VALUE_COUNT",missing_count)
                    already_missing[np.where(data_hist_all[i_mask-1]["count"][...]<count_thresh)]=True 

                else: #look forwards and backwards for other images
                    for key in data_all_adj[i_mask].keys():
                        data_all_adj[i_mask][key][(data_hist_all[i_mask]["count"][...]<count_thresh)& (already_missing==False)]=missing_count
                        data_all_adj[i_mask][key][(data_hist_all[i_mask-1]["count"][...]<count_thresh)& (already_missing==False)]=missing_count
                        data_all_adj[i_mask][key].attrs.create("MISSING_VALUE_COUNT",missing_count)
                    #add these count<missing_count values to already_missing
                    already_missing[np.where(data_hist_all[i_mask]["count"][...]<count_thresh)]=True
                    already_missing[np.where(data_hist_all[i_mask-1]["count"][...]<count_thresh)]=True 
 
                #2. Mask based on historical count rate
                if i_mask==0: #can only look at forward difference for first image of the day
                    #removal of upper outliers only
                    tomask=LST_diff[i_mask,...]>data_hist_all[i_mask] ["Mean"][...]+asig*data_hist_all[i_mask]["Stdev"][...]
                    tomask[already_missing==True]=False
                    for key in data_all_adj[i_mask].keys():
                        data_all_adj[i_mask][key][tomask==True]=missing_hist
                        data_all_adj[i_mask][key].attrs.create("MISSING_VALUE_HIST",missing_hist)

                elif i_mask==(len(toload)-1): #can only look at backward difference for last image of the day 
                    #removal of lower outliers only
                    tomask=LST_diff[i_mask-1,...]<data_hist_all[i_mask-1]["Mean"][...]-asig*data_hist_all[i_mask-1]["Stdev"][...]
                    tomask[already_missing==True]=False
                    for key in data_all_adj[i_mask].keys():
                        data_all_adj[i_mask][key][tomask==True]=missing_hist
                        data_all_adj[i_mask][key].attrs.create("MISSING_VALUE_HIST",missing_hist)

                else: #look forwards and backwards for other images
                    tomask_up=LST_diff[i_mask,...]>data_hist_all[i_mask]["Mean"][...]+asig*data_hist_all[i_mask]["Stdev"][...]
                    tomask_down=LST_diff[i_mask-1,...]<data_hist_all[i_mask-1]["Mean"][...]-asig*data_hist_all[i_mask-1]["Stdev"][...]
                    tomask_up[already_missing==True]=False
                    tomask_down[already_missing==True]=False
                    for key in data_all_adj[i_mask].keys():
                        data_all_adj[i_mask][key][tomask_up==True]=missing_hist
                        data_all_adj[i_mask][key][tomask_down==True]=missing_hist
                        data_all_adj[i_mask][key].attrs.create("MISSING_VALUE_HIST",missing_hist)

                #calculate total number of missing values for output in log file
                n_missing_old=np.size(np.where(data_all_adj[i_mask]["LST"][...]==data_all_adj[i_mask]["LST"].attrs["MISSING_VALUE"])[0])
                n_missing_adj=np.size(np.where(data_all_adj[i_mask]["LST"][...]==data_all_adj[i_mask]["LST"].attrs["MISSING_VALUE_ADJ"])[0])
                n_missing_count=np.size(np.where(data_all_adj[i_mask]["LST"][...]==data_all_adj[i_mask]["LST"].attrs["MISSING_VALUE_COUNT"])[0])
                n_missing_hist=np.size(np.where(data_all_adj[i_mask]["LST"][...]==data_all_adj[i_mask]["LST"].attrs["MISSING_VALUE_HIST"])[0])       
                n_tot=float(np.size(data_all_adj[i_mask]["LST"]))

                #simple plot for checking  
                if i_mask in [12,13,14,15]:  #plot images at 10:00, 10:15,10:30 and 10:45 for checking
                    F = plt.figure(figsize=(4, 2))
                    toplt=np.copy(data_all_adj[i_mask]["LST"][...]).astype(float)
                    toplt[toplt<-100]=np.nan
                    toplt=toplt/100
                    plt.imshow(toplt,vmin=10,vmax=50)
                    plt.title(outmask[-4:],pad=2)
                    plt.colorbar(shrink=0.6,extend="both")

                    pp.savefig()


                logfile.write(outmask[-4:]+":"+str(int(round(n_missing_old*100/n_tot))).zfill(2)+\
                                           ":"+str(int(round(n_missing_adj*100/n_tot))).zfill(2)+\
                                           ":"+str(int(round(n_missing_count*100/n_tot))).zfill(2)+\
                                           ":"+str(int(round(n_missing_hist*100/n_tot))).zfill(2)+\
                                            ",") #HHMM:existing missing: additional missing:missing count<count_thresh:hist missing,
            
            else:   #we are not masking as image is missing. No File existed to be opened.              
                logfile.write(outmask[-4:]+":NaN:NaN:NaN:NaN,") #HHMM:existing missing: additional missing

        for i_mask,outmask in enumerate(toload):            #loop over times of day and mask             
            if not type(data_all_adj[i_mask])==np.float:    #don't try to mask a missing image
                data_all_adj[i_mask].close()
                data_all[i_mask].close()
                if i_mask>0:
                    data_hist_all[i_mask-1].close()
    pp.close()
    logfile.write("\n")#new line for each date
logfile.close() 




