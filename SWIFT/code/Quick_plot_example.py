import numpy as np
import matplotlib as mpl
import glob,os,datetime,h5py,sys
import matplotlib.pyplot as plt
import mask_data_functions as fns

file_lat="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/ancillary/hdf5_lsasaf_msg_lat_msg-disk_4bytesprecision"
file_lon="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/ancillary/hdf5_lsasaf_msg_lon_msg-disk_4bytesprecision"
#fid=h5py.File("/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50/200608/HDF5_LSASAF_MSG_LST_MSG-Disk_200608011700")
fid=h5py.File("/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50/200607/HDF5_LSASAF_MSG_LST_MSG-Disk_200607311700")
dat_in=fid["lsta_av"]
dat_count=fid["lsta_av_count"]


toplt=dat_in[...].astype("float")
toplt[toplt==-9999]=np.nan
toplt=toplt/100

toplt_count=dat_count[...].astype("float")
toplt_count[toplt_count<=0]=np.nan

fns.regular_plot(file_lon, file_lat, toplt,toplt_count,"20060731 Example LSTA_AV",'/data/hmf/projects/NFLICS/test_comparison_lsta_av_20060731_testint.png')



F = plt.figure(figsize=(4,2))

file="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_sorted/200506/HDF5_LSASAF_MSG_LST_MSG-Disk_200506071030"
data=h5py.File(file)["LST"]
plt.subplot(3,2,1)
toplt=np.copy(data).astype(float)
toplt[toplt<-50]=np.nan
toplt=toplt/100
plt.imshow(toplt,vmin=10,vmax=50)
plt.title("103006/07/2015 Origional",pad=2)
plt.colorbar(shrink=0.6,extend="both")


file="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj/200506/HDF5_LSASAF_MSG_LST_MSG-Disk_200506071030"
data=h5py.File(file)["LST"]
plt.subplot(3,2,3)
toplt=np.copy(data).astype(float)
toplt[toplt<-50]=np.nan
toplt=toplt/100
plt.imshow(toplt,vmin=10,vmax=50)
plt.title("103006/07/2015 Mask ADJ",pad=2)
plt.colorbar(shrink=0.6,extend="both")


file="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask2/200506/HDF5_LSASAF_MSG_LST_MSG-Disk_200506071030"
data=h5py.File(file)["LST"]
plt.subplot(3,2,5)
toplt=np.copy(data).astype(float)
toplt[toplt<-50]=np.nan
toplt=toplt/100
plt.imshow(toplt,vmin=10,vmax=50)
plt.title("103006/07/2015 Mask Count<50 & 3 Sig",pad=2)
plt.colorbar(shrink=0.6,extend="both")


file="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_2sig/200506/HDF5_LSASAF_MSG_LST_MSG-Disk_200506071030"
data=h5py.File(file)["LST"]
plt.subplot(3,2,2)
toplt=np.copy(data).astype(float)
toplt[toplt<-50]=np.nan
toplt=toplt/100
plt.imshow(toplt,vmin=10,vmax=50)
plt.title("103006/07/2015 Mask 2 Sig",pad=2)
plt.colorbar(shrink=0.6,extend="both")


file="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_1sig/200506/HDF5_LSASAF_MSG_LST_MSG-Disk_200506071030"
data=h5py.File(file)["LST"]
plt.subplot(3,2,4)
toplt=np.copy(data).astype(float)
toplt[toplt<-50]=np.nan
toplt=toplt/100
plt.imshow(toplt,vmin=10,vmax=50)
plt.title("103006/07/2015 Mask 1 Sig",pad=2)
plt.colorbar(shrink=0.6,extend="both")
plt.tight_layout()

file="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_2sig_final/200506/HDF5_LSASAF_MSG_LST_MSG-Disk_200506071030"
data=h5py.File(file)["LST"]
plt.subplot(3,2,6)
toplt=np.copy(data).astype(float)
toplt[toplt<-50]=np.nan
toplt=toplt/100
plt.imshow(toplt,vmin=10,vmax=50)
plt.title("103006/07/2015 Mask 2 Sig FINAL",pad=2)
plt.colorbar(shrink=0.6,extend="both")

plt.savefig("/data/hmf/projects/NFLICS/200506071030_compare_masksd_FINAL.png")

toplt=np.copy(data).astype(float)/100

t1=np.ma.masked_less(toplt,-80)
t2=np.ma.masked_less(toplt,-70)
t3=np.ma.masked_less(toplt,-60)
t4=np.ma.masked_less(toplt,-50)
plt.imshow(t1,cmap='RdPu')
plt.imshow(t2,cmap='spring')
plt.imshow(t3,cmap='autumn')
plt.imshow(t4,cmap='viridis')
plt.title("10:30 06/07/2015 \n pink:origional, magenta:adj, red: N<50, navy: outside of 2sigma",pad=2)
plt.savefig("/data/hmf/projects/NFLICS/200506071030_mask_example.png")
