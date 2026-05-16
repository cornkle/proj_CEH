#########################
#Main code for creating a better cloud mask for SEVIRI LST data
#First pass: remove all points where adjacent in time images have missing data
#In practise this masks further from the cloud edges, identified as a problem area
#from visually inspecticng cloud masks
#########################
#Related scripts
# mask_data_wrap.py creates jobs (when adj=True)
# mask_data 
#########################
#consider each time seperately (there are 96 times of day, 1 every 15min)
#e.g consider 08:00 and nt times each side. Here nt=3
#
#  07:15 07:30 07:45 08:00 08:15 08:30 08:45
#   t-3   t-2   t-1    t    t+1   t+2   t+3
#
#1. Remove all points where any of the images nt*2+1 are missing (cloud mask start and end times often wrong)

#########################
#Create SRA 02/10/2019
#########################

import os
import datetime
import h5py
import sys
import numpy as np
import pandas as pd

if len(sys.argv)>1:
    nt_back=int(sys.argv[1])
    nt_for=int(sys.argv[2]) 
    year=int(sys.argv[3])
else:   
    nt_back=3 #mask any points in grid
    nt_for=1 #mask any points in grid
    year=2012

n_missing=5000
fstart="HDF5_LSASAF_MSG_LST_MSG-Disk_" #start of file names
logfile=open("/scratch/cornkle/Mask_data_adjacent_3back_1forward_"+str(year)+".log","w+")
fdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward/"

#dates: only need for days so times don't matter
start_date=datetime.datetime(year,5,1,8,0)
end_date=datetime.datetime(year,10,31,8,0)
dates_all=pd.date_range(start_date,end_date,freq="24h")
use_dates=[str(date).replace("-","").replace("T","").replace(":","")[:8] for date in dates_all]

#loop over all dates in year to be processed  
  
for use_date in use_dates:  #20120625, 20120626,20120627  [56:58]
    logfile.write(use_date+",")
    print(use_date)
    #get times of day and files to load
    start_time=datetime.datetime(int(use_date[:4]),int(use_date[4:6]),int(use_date[6:]),7,0)
    end_time=datetime.datetime(int(use_date[:4]),int(use_date[4:6]),int(use_date[6:]),19,45)
    dates_all=pd.date_range(start_time,end_time,freq="15min")
    toload=[date.replace("-","").replace("T","").replace(":","")[:6]+"/"+fstart+date.replace("-","").replace("T","").replace(":","")[:12] for date in np.datetime_as_string(dates_all)]

    #read in 1 day of data
    data_all=[] #this holds the file handles (these will be modified)
    LST_data=[] #this holds the actual LST data

    for f in toload:
        testf=os.system("h5ls "+fdir+f)
        #print(testf)
        if testf!=256:
            fid = h5py.File(fdir+f, "r+")  #assume files are already copied!!!!!
            if "MISSING_VALUE_ADJ" in fid["LST"].attrs.keys():
                sys.exit("ERROR! this data has already been masked:please check that the specified fdir is correct!")
            else: 
                data_all.append(fid)
                LST_data.append(fid["LST"])
        else:
            print("Skipping "+f)
            data_all.append(np.nan)  #don't want to mask when we have a full missing image in the image range
            LST_data.append(np.zeros((n_missing,n_missing)))

    #sort out sizes for missing data
    nx=np.amin([np.shape(d)[0] for d in LST_data])  #data ny unless all data are missing
    ny=np.amin([np.shape(d)[1] for d in LST_data])  #data ny unless all data are missing

    if (nx==n_missing)& (ny==n_missing):
        print("all missing data for ",date)
    else:
    #keep a copy of theh data to use for masking: must mask each time using the raw data
        LST_data=np.array([data[:nx,:ny] for data in LST_data]) #create numpy array
        #loop over dates and take 7 grids for each. combine mask and apply 
        #loop over times and create 
        for i_mask,outmask in enumerate(toload):
            if not type(data_all[i_mask])==np.float:
                if i_mask<(nt_back)+1:
                    i_min=0
                    i_max=i_mask+nt_for+1
                elif (len(toload)-i_mask-1)<nt_for:
                    i_min=i_mask-nt_back
                    i_max=len(toload)
                else:  
                    i_min=i_mask-nt_back
                    i_max=i_mask+nt_for+1
                #print("i_maxk, i_min, i_max ",i_mask,i_min,i_max)
                datasub=np.copy(LST_data[i_min:i_max]).astype(float)
                #print("shape datasub",np.shape(datasub))
                datasub[datasub<-100]=np.nan #mask all elements
                any_missing=np.any(np.isnan(datasub),axis=0)
                #add attribute for second missing value
                data_all[i_mask]["LST"].attrs.create("MISSING_VALUE_ADJ",(data_all[i_mask]["LST"].attrs["MISSING_VALUE"]+1000))
                data_all[i_mask]["Q_FLAGS"].attrs.create("MISSING_VALUE_ADJ",(data_all[i_mask]["Q_FLAGS"].attrs["MISSING_VALUE"]+1000))
                data_all[i_mask]["errorbar_LST"].attrs.create("MISSING_VALUE_ADJ",(data_all[i_mask]["errorbar_LST"].attrs["MISSING_VALUE"]+1000))
                #mask missing data values in fdir dataset. 
                data_all[i_mask]["LST"][(any_missing==True) & (data_all[i_mask]["LST"][...]!=data_all[i_mask]["LST"].attrs["MISSING_VALUE"])]=\
                                    data_all[i_mask]["LST"].attrs["MISSING_VALUE"]+1000
                data_all[i_mask]["Q_FLAGS"][(any_missing==True) & (data_all[i_mask]["Q_FLAGS"][...]!=data_all[i_mask]["Q_FLAGS"].attrs["MISSING_VALUE"])]=\
                                    data_all[i_mask]["Q_FLAGS"].attrs["MISSING_VALUE"]+1000   
                data_all[i_mask]["errorbar_LST"][(any_missing==True) & (data_all[i_mask]["errorbar_LST"][...]!=data_all[i_mask]["errorbar_LST"].attrs["MISSING_VALUE"])]=\
                                    data_all[i_mask]["errorbar_LST"].attrs["MISSING_VALUE"] +1000

                n_missing_old=np.size(np.where(data_all[i_mask]["LST"][...]==data_all[i_mask]["LST"].attrs["MISSING_VALUE"])[0])
                n_missing_add=np.size(np.where(data_all[i_mask]["LST"][...]==data_all[i_mask]["LST"].attrs["MISSING_VALUE_ADJ"])[0])
                n_tot=float(np.size(data_all[i_mask]["LST"]))
               # print(n_missing_old,n_missing_add,n_tot)
                data_all[i_mask].close()
                logfile.write(outmask[-4:]+":"+str(int(round(n_missing_old*100/n_tot))).zfill(2)+":"+str(int(round(n_missing_add*100/n_tot))).zfill(2)+",") #HHMM:existing missing: additional missing,
            else:                
                logfile.write(outmask[-4:]+":NaN:NaN,") #HHMM:existing missing: additional missing,

    logfile.write("\n")#new line for each date

logfile.close() 




