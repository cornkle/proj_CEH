
##########################
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
import glob,os,datetime,h5py,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.backends.backend_pdf import PdfPages

if len(sys.argv)>1:
    year=int(sys.argv[1])               #one job per year (loops over dates)
else:   
    year=2012


fdir_out="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50/"
fdir_hist="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_clim/"
plotdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_clim_anom/"

fstart="HDF5_LSASAF_MSG_LST_MSG-Disk_" #start of file names
fstart_stats="HDF5_LSASAF_MSG_LST_MSG-Disk_HistStats_"
#logfile to output the number of missing points
n_missing=5000

#dates: only need for days so times don't matter
start_date=datetime.datetime(year,5,1,8,0)
end_date=datetime.datetime(year,10,31,8,0)
#end_date=datetime.datetime(year,5,1,8,0) #TESTING!
dates_all=pd.date_range(start_date,end_date,freq="24h")
use_dates=[str(date).replace("-","").replace("T","").replace(":","")[:8] for date in dates_all]

#loop over all dates in year to be processed  
for use_date in use_dates:  #20120625, 20120626,20120627  [56:58]

    pp = PdfPages(plotdir+"anom_LST_"+use_date+'.pdf')
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams.update({'font.size': 8})

    print(use_date)   
    #get times of day and files to load
    start_time=datetime.datetime(int(use_date[:4]),int(use_date[4:6]),int(use_date[6:]),7,0)
    end_time=datetime.datetime(int(use_date[:4]),int(use_date[4:6]),int(use_date[6:]),19,45)

    dates_all=pd.date_range(start_time,end_time,freq="15min")
    #print(dates_all)
    toload=[date.replace("-","").replace("T","").replace(":","")[:6]+"/"+\
            fstart+date.replace("-","").replace("T","").replace(":","")[:12] for date in np.datetime_as_string(dates_all)]
    toload_stats=[np.datetime_as_string(np.datetime64(date)).replace("-","").replace("T","").replace(":","")[4:6]+"/"+\
                  fstart_stats+np.datetime_as_string(np.datetime64(date)).replace("-","").replace("T","").replace(":","")[4:8]+"_"+\
                  np.datetime_as_string(np.datetime64(date)).replace("-","").replace("T","").replace(":","")[8:12]\
                  for date in dates_all[0:]]
    #print(toload)
    #print(toload_stats)
    anom_all=[] 
    data_all_adj=[]
    data_clim_all=[]
    for i_f,f in enumerate(toload): #loop over times of day
        #read in the stats data (on for each time averaged over years )
        fid_hist= h5py.File(fdir_hist+toload_stats[i_f], "r")
        data_clim_all.append(fid_hist)
        #read in the lst image data
        testf=os.system("h5ls "+fdir_out+f)
        #if testf!=256:        
        if testf==0:
            fid_adj = h5py.File(fdir_out+f, "r+")  #assume files are already copied!!!
            data_all_adj.append(fid_adj)
        else:
            data_all_adj.append(np.nan)  #don't want to mask when we have a full missing image in the image rang

        
    for i_mask,outmask in enumerate(toload):            #loop over times of day and calculate anomolies 
        #print(i_mask)
        print(data_all_adj[i_mask])
        if not type(data_all_adj[i_mask])==np.float:    #don't try to do anything to a missing image
            #create anomoly dataset
            im_now=np.copy(data_all_adj[i_mask]["LST"][...]).astype(float)
            im_now[np.where(im_now<-1000)]=np.nan
            clim_now=np.copy(data_clim_all[i_mask]["Mean"][...]).astype(float)
            toout_im=im_now-clim_now
            toout_im[np.where(np.isnan(toout_im)==True)]=-9999

            #running sum over the day of anomolies (to calculate average anomoly upto this time)   
            tot=np.copy(toout_im)
            count=np.ones(np.shape(toout_im))
            count[np.where(toout_im==-9999)]=0
            tot[np.where(toout_im<=-9999)]=0.0

            if i_mask==0:
                running_tot=np.copy(tot)
                running_count=np.copy(count)
            else:                
                running_tot=running_tot+np.copy(tot)
                running_count=running_count+np.copy(count)

            #mask for this particular image
            toout_av=np.copy(running_tot)
            toout_av[np.where(running_count==0)]=np.nan #all cloud mask
            #toout_av[np.where(toout_im<=-9999)]=np.nan   #masked in the ith image
            
            #mask from previous (i-1) image
            #if (i_mask>=1): 
            #    if not type(data_all_adj[i_mask-1])==np.float:
            #        toout_av[data_all_adj[i_mask-1]["lsta"][...]<-1000]=np.nan 
            #    else:
            #        toout_av[...]=np.nan
            #            #mask from previous (i-2) image 
            #if (i_mask>=2):
            #    if not type(data_all_adj[i_mask-2])==np.float:
            #        toout_av[data_all_adj[i_mask-2]["lsta"][...]<-1000]=np.nan
            #    else:
            #        toout_av[...]=np.nan
                
            toout_av=toout_av/(running_count)
            toout_av[np.where(np.isnan(toout_av)==True)]=-9999  
           
            if not "lsta" in data_all_adj[i_mask].keys():      
                data_all_adj[i_mask].create_dataset("lsta",data=toout_im.astype("int16"))
                data_all_adj[i_mask]["lsta"].attrs.create("MISSING_VALUE",-9999)
                data_all_adj[i_mask]["lsta"].attrs.create("SCALING_FACTOR",100.0)
                data_all_adj[i_mask]["lsta"].attrs.create("OFFSET",0.0)
            else:
                data_all_adj[i_mask]["lsta"][...]=toout_im.astype("int16")

            if not "lsta_av" in data_all_adj[i_mask].keys():      
                data_all_adj[i_mask].create_dataset("lsta_av",data=toout_av.astype("int16"))
                data_all_adj[i_mask]["lsta_av"].attrs.create("MISSING_VALUE",-9999)
                data_all_adj[i_mask]["lsta_av"].attrs.create("SCALING_FACTOR",100.0)
                data_all_adj[i_mask]["lsta_av"].attrs.create("OFFSET",0.0)
            else:
                data_all_adj[i_mask]["lsta_av"][...]=toout_av.astype("int16")

            if not "lsta_av_count" in data_all_adj[i_mask].keys():      
                data_all_adj[i_mask].create_dataset("lsta_av_count",data=running_count.astype("int16"))
                data_all_adj[i_mask]["lsta_av_count"].attrs.create("MISSING_VALUE",-9999)
                data_all_adj[i_mask]["lsta_av_count"].attrs.create("SCALING_FACTOR",0.0)
                data_all_adj[i_mask]["lsta_av_count"].attrs.create("OFFSET",0.0)
            else:
                data_all_adj[i_mask]["lsta_av_count"][...]=running_count.astype("int16")





          #simple plot for checking  
            if i_mask==40:  
                F = plt.figure(figsize=(4, 4))
                plt.subplot(2,1,1)
                toplt=np.copy(toout_im)
                toplt[np.where(toplt<-1000)]=np.nan
                toplt=toplt/100
                plt.imshow(toplt,vmin=-10,vmax=10,cmap="bwr",interpolation="nearest")
                plt.title("LST anomaly "+ outmask[-4:],pad=2)
                plt.colorbar(shrink=0.6,extend="both")

                plt.subplot(2,1,2)
                toplt=np.copy(toout_av)
                toplt[np.where(toplt<-1000)]=np.nan
                toplt=toplt/100
                plt.imshow(toplt,vmin=-10,vmax=10,cmap="bwr",interpolation="nearest")
                plt.title("LST anomoly av upto "+outmask[-4:],pad=2)
                plt.colorbar(shrink=0.6,extend="both")
                pp.savefig()

    for i_mask,outmask in enumerate(toload):            #loop over times of day and mask             
        if not type(data_all_adj[i_mask])==np.float:    #don't try to mask a missing image
            data_all_adj[i_mask].close()
            data_clim_all[i_mask-1].close()
    pp.close()
    plt.close()

