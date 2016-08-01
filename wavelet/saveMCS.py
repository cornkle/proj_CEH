# -*- coding: utf-8 -*-


import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage.measurements import label 
import datetime as dt
from eod import read_eod as re
import xarray as xr
import os
from wavelet import util
from skimage.feature import match_template
from utils import u_arrays as ua
import pickle as pkl
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import cartopy.crs as ccrs
from scipy import ndimage
from mapping import creategrid as cg

HOD=range(24)   # hours of day
YRANGE=range(2004,2014)

def savetMCS():
   
    trmm_folder= "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder='/users/global/cornkle/data/OBS/meteosat_SA15'
    
    # make a salem grid
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')

    t=re.trmm(trmm_folder, yrange=YRANGE, area=[-10, 10, 10, 20]) 
    m=re.msg(msg_folder)
    
    cnt=0    
    
    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered      
    for _y, _m, _d, _h, _mi in zip(t.dates.y, t.dates.m, t.dates.d, t.dates.h, t.dates.mi):
                      
       tdic=t.getDData(_y, _m, _d, _h, _mi, cut=[9,21])
                                   
       #define the "0 lag" frist
       arr=np.array([15,30,45,60, 0])
       dm = arr - _mi
       ind=(np.abs(dm)).argmin()
       
       # set zero shift time for msg
       date=dt.datetime(_y, _m, _d, _h, _mi)
       print(date)
       
       dt0=dm[ind] 
       ndate = date + dt.timedelta(minutes=int(dt0))                                           
       mdic=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[tdic['lon'].min(), tdic['lon'].max(), tdic['lat'].min(), tdic['lat'].max() ])        
       if not mdic:
           print('Date missing')
           continue
       lon1=mdic['lon']
       lat1=mdic['lat']
       mdic['t'][mdic['t']>-40]=0
       labels, numL = label(mdic['t'])
           
       u , inv = np.unique(labels, return_inverse=True)
       n = np.bincount(inv)
          
       goodinds = u[n>2500]  # all blobs with more than 2500 pixels - size threshold
       print(goodinds)
       if not sum(goodinds) > 0:
            continue
     
       for gi in goodinds:
              if gi == 0:    # index 0 is always background, ignore!
                  continue
              
              inds = np.where(labels == gi)
 
              # cut a box for every single blob from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
              latmax, latmin = mdic['lat'][inds].max() , mdic['lat'][inds].min()
              lonmax, lonmin = mdic['lon'][inds].max() , mdic['lon'][inds].min()
              mmeans=np.percentile(mdic['t'][inds], 90)
              td = t.getDData(_y, _m, _d, _h, _mi, cut=[latmin-0.2, latmax+0.2])
              
              #ensure minimum trmm rainfall in area
              if len(np.where(td['p']>0)[0]) < 100:      #  at least 100 pixel with rainfall 
                  print('Kickout: TRMM min pixel = 100')
                  continue              
              
              dt0=dm[ind] 
              ndate = date + dt.timedelta(minutes=int(dt0) )                
             # print('Date1', ndate)                              
              ml0=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[lonmin-0.3, lonmax+0.3, latmin-0.25, latmax+0.25])        
              if not ml0:
                  continue
              
              dt1=dm[ind]-15
              ndate = date + dt.timedelta(minutes=int(dt1))       
           #   print('Date2', ndate)  
              ml1=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[lonmin-0.3, lonmax+0.3, latmin-0.25, latmax+0.25])               
              if not ml1:
                  continue

              dt2=dm[ind]-30
              ndate = date + dt.timedelta(minutes=int(dt2))       
           #   print('Date2', ndate)  
              ml2=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[lonmin-0.3, lonmax+0.3, latmin-0.25, latmax+0.25])               
              if not ml1:
                  continue

              dt3=dm[ind]-45
              ndate = date + dt.timedelta(minutes=int(dt3))       
           #   print('Date2', ndate)  
              ml3=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[lonmin-0.3, lonmax+0.3, latmin-0.25, latmax+0.25])               
              if not ml1:
                  continue   

              dtx=dm[ind]+45
              ndate = date + dt.timedelta(minutes=int(dtx))       
           #   print('Date2', ndate)  
              mlx=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[lonmin-0.3, lonmax+0.3, latmin-0.25, latmax+0.25])               
              if not ml1:
                  continue                      
              
              #create grid surrounding the blob
              # Transform lon, lats to the mercator projection
              x, y = pyproj.transform(salem.wgs84, proj, ml0['lon'], ml0['lat'])
              # take the min and max
              xmax, xmin = np.max(x), np.min(x)
              ymax, ymin = np.max(y), np.min(y)
              # Count the number of pixels
              dx = 5000
              nx, r = divmod(xmax - xmin, dx)
              ny, r = divmod(ymax - ymin, dx)
              # Here one could add + 1 to be sure that the last pixel is always included
              grid = salem.Grid(nxny=(nx, ny), dxdy=(dx, dx), ll_corner=(xmin, ymin), proj=proj)                 
            
              # interpolate TRM and MSG to salem grid
              xi, yi = grid.ij_coordinates
              lon, lat = grid.ll_coordinates

              # Transform lons, lats to grid
              xm, ym = grid.transform( ml0['lon'].flatten(), ml0['lat'].flatten(), crs=salem.wgs84)
              xt, yt = grid.transform(td['lon'].flatten(), td['lat'].flatten(), crs=salem.wgs84)
        
              # Convert for griddata input 
              mpoints = np.array((ym, xm)).T
              tpoints = np.array((yt, xt)).T
              inter = np.array((np.ravel(yi), np.ravel(xi))).T
        
              # Interpolate using delaunay triangularization 
              dummy = griddata(mpoints, ml0['t'].flatten(), inter, method='linear')
              dummy = dummy.reshape((grid.ny, grid.nx))
              outl = np.full_like(dummy, np.nan)
              xl, yl = grid.transform(lon1[inds], lat1[inds], crs=salem.wgs84, nearest=True, maskout=True)
                            
                            
              print('Lag0')
              outl[yl.compressed(),xl.compressed()] = dummy[yl.compressed(), xl.compressed()]
        
              # Interpolate using delaunay triangularization 
              dummyt = griddata(tpoints, td['p'].flatten(), inter, method='linear')
              outt = dummyt.reshape((grid.ny, grid.nx))      
              if len(np.where(outt>0)[0]) < 100:   #  at least 100 pixel with rainfall 
                  print('Kickout: TRMM wavelet min pixel pcp = 100')
                  continue  
                                        
              tmask = np.isfinite(outt)
              mmask = np.isfinite(outl)
              mask2 = np.isfinite(outl[tmask])
              
              if sum(mask2.flatten()) < 200:#sum(mmask.flatten())*0.3:
                  print('Kickout: TRMM MSG overlap less than 0.3 of cloud area')
                  continue                             
              
              print('Hit:', gi)
                                     
              # lag -1
             
               # Interpolate using delaunay triangularization 
              dummy1 = griddata(mpoints, ml1['t'].flatten(), inter, method='linear')
              dummy1 = dummy1.reshape((grid.ny, grid.nx))
              outl1 = np.full_like(dummy1, np.nan)
              print('Lag1')
              outl1[yl.compressed(),xl.compressed()] = dummy1[yl.compressed(), xl.compressed()]        
              
               # lag -2
             
               # Interpolate using delaunay triangularization 
              dummy2 = griddata(mpoints, ml2['t'].flatten(), inter, method='linear')
              dummy2 = dummy2.reshape((grid.ny, grid.nx))
              outl2 = np.full_like(dummy2, np.nan)
              print('Lag2')
              outl2[yl.compressed(),xl.compressed()] = dummy2[yl.compressed(), xl.compressed()] 
              
              # lag -3
             
               # Interpolate using delaunay triangularization 
              dummy3 = griddata(mpoints, ml3['t'].flatten(), inter, method='linear')
              dummy3 = dummy3.reshape((grid.ny, grid.nx))
              outl3 = np.full_like(dummy3, np.nan)
              print('Lag3')
              outl3[yl.compressed(),xl.compressed()] = dummy3[yl.compressed(), xl.compressed()] 
              
               # lag x
             
               # Interpolate using delaunay triangularization 
              dummyx = griddata(mpoints, mlx['t'].flatten(), inter, method='linear')
              dummyx = dummyx.reshape((grid.ny, grid.nx))
              outlx = np.full_like(dummyx, np.nan)
              print('Lagx')
              outlx[yl.compressed(),xl.compressed()] = dummyx[yl.compressed(), xl.compressed()] 
                                                               
              
              da = xr.Dataset({'p' : (['x', 'y'], outt),  
                                 't_lag0'  : (['x', 'y'], dummy),
                                 'tc_lag0'  : (['x', 'y'], outl),
                                 't_lag1'  : (['x', 'y'], dummy1),
                                 'tc_lag1'  : (['x', 'y'], outl1),
                                 't_lag2'  : (['x', 'y'], dummy2),
                                 'tc_lag2'  : (['x', 'y'], outl2),
                                 't_lag3'  : (['x', 'y'], dummy3),
                                 'tc_lag3'  : (['x', 'y'], outl3),
                                 't_lagx'  : (['x', 'y'], dummyx),
                                 'tc_lagx'  : (['x', 'y'], outlx),
                                 'tmask'   : (['x', 'y'], mmask.astype(int)),
                                 'pmask'   : (['x', 'y'], tmask.astype(int))}, 
                            coords= {'lon' : (['x', 'y'], lon),
                                     'lat' : (['x', 'y'], lat),                                    
                                     'time': date })
              da.attrs['lag0']=dt0
              da.attrs['lag1']=dt1
              da.attrs['meanT']=np.mean(outl[mmask])
              da.attrs['T90perc']=mmeans
              da.attrs['meanT_cut']=np.mean(outl[tmask][mask2])
              da.attrs['area']=sum(mmask.flatten())
              da.attrs['area_cut']=sum(mask2)
              da.close()  
              savefile = '/users/global/cornkle/MCSfiles/'+date.strftime('%Y-%m-%d_%H:%M:%S')+'_'+str(gi)+'.nc'
              try:
                  os.remove(savefile)
              except OSError:
                  pass
              da.to_netcdf(path=savefile, mode='w')
              print('Saved '+savefile)
                       
              cnt=cnt+1
             
    print('Saved '+str(cnt)+' MCSs as netcdf.')    


def readMCS_relateScales():
             
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/')
    
#    arr=np.array([15,   16,   17,   18,   19,   20,   21,   22,   24,
#         25,   27,   28,   30,   32,   34,   36,   38,   40,
#         42,   45,   48,   50,   53,   57,   60,   64,   67,
#         71,   76,   80,   85,   90,   95,  101,  107,  113,
#        120,  127,  135,  143,  151,  160,  170,  180,  190,  202], dtype=str)
    
    arr=np.array([15,   16,   17,   18,   19,   20,   21,   22,   24,
         25,   27,   28,   30,   32,   34,   36,   38,   40,
         42,   45,   48,   50,   53,   57,   60], dtype=str)

    
    wave={}
    for a in arr:    
        wave[a]={}
        for b in arr:
            wave[a][b]={}
            wave[a][b]['corr']=[]
            wave[a][b]['x']=[]
            wave[a][b]['y']=[]
            wave[a][b]['torig']=[]
            wave[a][b]['porig']=[]
            wave[a][b]['twavelet']=[]
            wave[a][b]['pwavelet']=[]
            
                    
 #   return wave   
    
    for f in files:
        print('Doing file: '+f)
        dic = xr.open_dataset(f)
               
        outt1=dic['tc_lag0'].values.copy()
        outt2=dic['tc_lag1'].values.copy()
        
        mmeans=np.percentile(outt1[np.isfinite(outt1)], 30)
        
        maxi=np.nanmin(outt1)
        thresh=maxi+15
        
        
        outp=dic['p'].values.copy()
        mmeans=dic.meanT_cut
        outp[np.isnan(outp)]=-10**-5
        outt1[np.isnan(outt1)]=mmeans
        outt1[outt1>thresh]=mmeans
        outt2[np.isnan(outt2)]=mmeans
        outt2[outt2>thresh]=mmeans
        
        wav1 = util.waveletTP(outt1, outp, 5)
        wav2 = util.waveletTP(outt2, outp, 5)
        
        for pos in range(arr.size):  

            print(arr[pos])                      
        
            tt1=wav1['t'][pos,:,:]
            tt2=wav2['t'][pos,:,:]
            
            tt1[np.where(dic['pmask'].values==0)]=0
            tt2[np.where(dic['pmask'].values==0)]=0
            
            for ppos in range(arr.size):  
                
                pp=wav1['p'][ppos,:, :]
        
                corr1=match_template(tt1, pp[5:-5, 5:-5])  
                y1,x1 = np.unravel_index(np.argmax(corr1), corr1.shape)              
                corr2=match_template(tt2, pp[5:-5, 5:-5])
                y2,x2 = np.unravel_index(np.argmax(corr2), corr2.shape)            

         #   if (np.max(corr1)>np.max(corr2)) & y1!=0 & x1!=0 & y1!=10 & x1!=10: 
         #       wave[arr[pos]].append(np.max(corr1))
         #   else:
         #       if y2!=0 & x2!=0 & y2!=10 & x2!=10:
         #           wave[arr[pos]].append(np.max(corr2))
         #       else:
         #           wave[arr[pos]].append(np.nan)

                if (np.max(corr1)>np.max(corr2)):
                    wave[arr[pos]][arr[ppos]]['corr'].append(np.max(corr1))  # correlation value
                    wave[arr[pos]][arr[ppos]]['x'].append(x1)        # x value of cross correlaion to compute shift 
                    wave[arr[pos]][arr[ppos]]['y'].append(y1)          # y value of cross correlaion to compute shift 
                    wave[arr[pos]][arr[ppos]]['torig'].append(outt1)    # original t map with masking (used for wavelet)
                    wave[arr[pos]][arr[ppos]]['porig'].append(outp)     # original p map with masking (used for wavelet)
                    wave[arr[pos]][arr[ppos]]['twavelet'].append(tt1)   # wavelet transformed t
                    wave[arr[pos]][arr[ppos]]['pwavelet'].append(pp)    # wavelet transformed p
                else:     
                    wave[arr[pos]][arr[ppos]]['corr'].append(np.max(corr2))     # same as above but for -15 minutes MSG if cross correlation is higher
                    wave[arr[pos]][arr[ppos]]['x'].append(x2)
                    wave[arr[pos]][arr[ppos]]['y'].append(y2)
                    wave[arr[pos]][arr[ppos]]['torig'].append(outt2)
                    wave[arr[pos]][arr[ppos]]['porig'].append(outp)
                    wave[arr[pos]][arr[ppos]]['twavelet'].append(tt2)
                    wave[arr[pos]][arr[ppos]]['pwavelet'].append(pp)            
    
    pkl.dump(wave, open('/users/global/cornkle/MCSfiles/save/MCScorr_scales_allyears.p', 'wb'))
    
    print('Saved!')

   
def readMCS_getWavelet():
             
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/')
    
#    arr=np.array([15,   16,   17,   18,   19,   20,   21,   22,   24,
#         25,   27,   28,   30,   32,   34,   36,   38,   40,
#         42,   45,   48,   50,   53,   57,   60,   64,   67,
#         71,   76,   80,   85,   90,   95,  101,  107,  113,
#        120,  127,  135,  143,  151,  160,  170,  180,  190,  202], dtype=str)
    
    arr=np.array([10, 11, 11, 12,13, 13, 14,15,   16,   17,   18,   19,   20,   21,   22,   24,
         25,   27,   28,   30,   32,   34,   36,   38,   40,
         42,   45,   48,   50,   53,   57,   60], dtype=str)
         
    ext=np.array([10,15,20,30,60], dtype=str)      
    rrange=list(range(arr.size))
    scales=np.in1d(arr, ext)
    rpos=np.array(rrange)[scales]
    
    print(rpos)
    
    wave={}
    wave['torig0_min']=[]
    wave['torig1_min']=[]
    wave['torig2_min']=[]
    wave['torig3_min']=[]
    wave['torigx_min']=[]
    wave['porig_max']=[]
    
    wave['torig0']=[]    
    wave['torig1']=[]
    wave['torig2']=[]
    wave['torig3']=[]
    wave['torigx']=[]
    wave['porig']=[]
            
    for a in ext:              
        wave[a]={}           
        
        wave[a]['twavelet0']=[]
        wave[a]['twavelet1']=[]
        wave[a]['twavelet2']=[]
        wave[a]['twavelet3']=[]
        wave[a]['twaveletx']=[]
        wave[a]['pwavelet']=[]           
        
        wave[a]['twavelet0_max']=[]
        wave[a]['twavelet1_max']=[]
        wave[a]['twavelet2_max']=[]
        wave[a]['twavelet3_max']=[]
        wave[a]['twaveletx_max']=[]
        wave[a]['pwavelet_max']=[]
               
 #   return wave   
    
    for f in files:
        print('Doing file: '+f)
        dic = xr.open_dataset(f)
        
     #   if (dic['time.hour'].values<15) or (dic['time.hour'].values>21):
          # print('smaller') 
     #      continue       
        
        outt0=np.array(dic['tc_lag0'].values.copy())
        outt1=np.array(dic['tc_lag1'].values.copy())
        outt2=np.array(dic['tc_lag2'].values.copy())
        outt3=np.array(dic['tc_lag3'].values.copy())
        outtx=np.array(dic['tc_lagx'].values.copy())
        
        outp=np.array(dic['p'].values.copy())
        outp[np.isnan(outp)]=-10**-5
        wave['porig_max'].append(np.percentile(outp[outp>1], 99))
        wave['porig'].append(outp)
        
        looparr=[outt0, outt1, outt2, outt3, outtx]
        strarr=['0', '1', '2', '3', 'x']
    
        for outt, strr in zip(looparr, strarr):
            
            wave['torig'+strr].append(outt)    
                        
            outt[np.isnan(outt)]=-40

            wav = util.waveletTP(outt, outp, 5) 
                
            outt[np.where(dic['pmask'].values==0)]=-40
     
            wave['torig'+strr+'_min'].append(np.percentile(outt0[outt0<-40], 1))        
      
              
            for pos in rpos:  

                print(arr[pos])            
                tt=np.array(wav['t'][pos,:,:])                  
                
                if strr=='0':
                    pp=np.array(wav['p'][pos,:,:])
                    wave[arr[pos]]['pwavelet'].append(pp)  
                    wave[arr[pos]]['pwavelet_max'].append(np.percentile(pp[pp>0], 99))
            
                tt[np.where(dic['pmask'].values==0)]=0                                                                 
                wave[arr[pos]]['twavelet'+strr].append(tt)         
                wave[arr[pos]]['twavelet'+strr+'_max'].append(np.percentile(tt[tt>0], 99))                                                 
    
    pkl.dump(wave, open('/users/global/cornkle/MCSfiles/save/MCS_wavelet_allyears_perc.p', 'wb'))
    
    print('Saved!')   
    

def readMCS_getWavelet_tresh():
             
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/')
    
#    arr=np.array([15,   16,   17,   18,   19,   20,   21,   22,   24,
#         25,   27,   28,   30,   32,   34,   36,   38,   40,
#         42,   45,   48,   50,   53,   57,   60,   64,   67,
#         71,   76,   80,   85,   90,   95,  101,  107,  113,
#        120,  127,  135,  143,  151,  160,  170,  180,  190,  202], dtype=str)
    
    arr=np.array([10, 11, 11, 12,13, 13, 14,15,   16,   17,   18,   19,   20,   21,   22,   24,
         25,   27,   28,   30,   32,   34,   36,   38,   40,
         42,   45,   48,   50,   53,   57,   60], dtype=str)
         
    ext=np.array([10,15,20,30,60], dtype=str)      
    rrange=list(range(arr.size))
    scales=np.in1d(arr, ext)
    rpos=np.array(rrange)[scales]
    
    print(rpos)
    
    wave={}
    wave['torig0_min']=[]
    wave['torig1_min']=[]
    wave['torig2_min']=[]
    wave['torig3_min']=[]
    wave['torigx_min']=[]
    wave['porig_max']=[]
    
    wave['torig0']=[]    
    wave['torig1']=[]
    wave['torig2']=[]
    wave['torig3']=[]
    wave['torigx']=[]
    wave['porig']=[]
            
    for a in ext:              
        wave[a]={}           
        
        wave[a]['twavelet0']=[]
        wave[a]['twavelet1']=[]
        wave[a]['twavelet2']=[]
        wave[a]['twavelet3']=[]
        wave[a]['twaveletx']=[]
        wave[a]['pwavelet']=[]           
        
        wave[a]['twavelet0_max']=[]
        wave[a]['twavelet1_max']=[]
        wave[a]['twavelet2_max']=[]
        wave[a]['twavelet3_max']=[]
        wave[a]['twaveletx_max']=[]
        wave[a]['pwavelet_max']=[]
               
 #   return wave   
    
    for f in files:
        print('Doing file: '+f)
        dic = xr.open_dataset(f)
        
     #   if (dic['time.hour'].values<15) or (dic['time.hour'].values>21):
          # print('smaller') 
     #      continue       
        
        outt0=np.array(dic['tc_lag0'].values.copy())
        outt1=np.array(dic['tc_lag1'].values.copy())
        outt2=np.array(dic['tc_lag2'].values.copy())
        outt3=np.array(dic['tc_lag3'].values.copy())
        outtx=np.array(dic['tc_lagx'].values.copy())
    
        mmeans=np.percentile(outt0[np.isfinite(outt0)], 30)
        
        print(mmeans)
        
        maxi=np.nanmin(outt0)
        thresh=maxi+15        
        
        outp=np.array(dic['p'].values.copy())
     
        outp[np.isnan(outp)]=-10**-5
        
        outt0[np.isnan(outt0)]=mmeans   # -40
        outt0[outt0>thresh]=mmeans
        outt1[np.isnan(outt1)]=mmeans
        outt1[outt1>thresh]=mmeans
        outt2[np.isnan(outt2)]=mmeans
        outt2[outt2>thresh]=mmeans
        outt3[np.isnan(outt3)]=mmeans
        outt3[outt3>thresh]=mmeans
        outtx[np.isnan(outtx)]=mmeans 
        outtx[outtx>thresh]=mmeans
    #    outt2[outt2>thresh]=-40#mmeans
        
        if np.mean(outt0)==mmeans:
            continue
        if np.mean(outt1)==mmeans:
            continue
        if np.mean(outt2)==mmeans:
            continue
        if np.mean(outt3)==mmeans:
            continue
        if np.mean(outtx)==mmeans:
            continue
        
        if not outtx[outtx<mmeans].any():
            continue
        
        print('Wavelet start')
                       
        wav0 = util.waveletTP(outt0, outp, 5) 
        
        wav1 = util.waveletTP(outt1, outp, 5)
        
        wav2 = util.waveletTP(outt2, outp, 5)

        wav3 = util.waveletTP(outt3, outp, 5)
 
        wavx = util.waveletTP(outtx, outp, 5)

        #print(wav1['scales'])
  
        
      #  outt0[np.where(dic['pmask'].values==0)]=mmeans
      #  outt1[np.where(dic['pmask'].values==0)]=mmeans   
      #  outt2[np.where(dic['pmask'].values==0)]=mmeans
     #   outt3[np.where(dic['pmask'].values==0)]=mmeans      
     #   outtx[np.where(dic['pmask'].values==0)]=mmeans      
        
        wave['torig0_min'].append(np.percentile(outt0[outt0<mmeans], 1))        
        wave['torig1_min'].append(np.percentile(outt1[outt1<mmeans], 1))
        wave['torig2_min'].append(np.percentile(outt2[outt2<mmeans], 1))
        wave['torig3_min'].append(np.percentile(outt3[outt3<mmeans], 1))
        wave['torigx_min'].append(np.percentile(outtx[outtx<mmeans], 1))
        wave['porig_max'].append(np.percentile(outp[outp>1], 99))
    
        wave['torig0'].append(outt0)
        wave['torig1'].append(outt1)
        wave['torig2'].append(outt2)
        wave['torig3'].append(outt3)
        wave['torigx'].append(outtx)
        wave['porig'].append(outp)
            
        
        for pos in rpos:  

            print(arr[pos])            
            tt0=np.array(wav0['t'][pos,:,:])                  
            tt1=np.array(wav1['t'][pos,:,:])
            tt2=np.array(wav2['t'][pos,:,:])
            tt3=np.array(wav3['t'][pos,:,:])
            ttx=np.array(wavx['t'][pos,:,:])
            pp=np.array(wav2['p'][pos,:,:])
            
            tt0[np.where(dic['pmask'].values==0)]=0                      
            tt1[np.where(dic['pmask'].values==0)]=0
            tt2[np.where(dic['pmask'].values==0)]=0      
            tt3[np.where(dic['pmask'].values==0)]=0
            ttx[np.where(dic['pmask'].values==0)]=0                                
           
            wave[arr[pos]]['twavelet0'].append(tt0)            
            wave[arr[pos]]['twavelet1'].append(tt1)
            wave[arr[pos]]['twavelet2'].append(tt2)
            wave[arr[pos]]['twavelet3'].append(tt3)
            wave[arr[pos]]['twaveletx'].append(ttx)
            wave[arr[pos]]['pwavelet'].append(pp)     
            
            wave[arr[pos]]['twavelet0_max'].append(np.percentile(tt0[tt0>0], 99))
            wave[arr[pos]]['twavelet1_max'].append(np.percentile(tt1[tt1>0], 99))
            wave[arr[pos]]['twavelet2_max'].append(np.percentile(tt2[tt2>0], 99))
            wave[arr[pos]]['twavelet3_max'].append(np.percentile(tt3[tt3>0], 99))
            wave[arr[pos]]['twaveletx_max'].append(np.percentile(ttx[ttx>0], 99))
            wave[arr[pos]]['pwavelet_max'].append(np.percentile(pp[pp>0], 99))
    
    pkl.dump(wave, open('/users/global/cornkle/MCSfiles/save/MCS_wavelet_allyears_perc_thresh.p', 'wb'))
    
    print('Saved!')  
    
def readMCS_getWavelet_label():
             
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/')
    
    wave={}
    
    strarr=['0', '1', '2', '3', 'x']
       
    wave['porig']=[]
    wave['pw']=[]
    wave['tw0']=[]
   # wave['scales']=[]
    
    for st in strarr:    
        wave['torig'+st]=[]         
        wave['tw'+st+'_max']=[]     # wavelet value at max point
        wave['pw'+st+'_max']=[]    # max p wavelet in radius      
        wave['p'+st+'_max']=[]    # max p  in radius       
        wave['p'+st+'_mean']=[]    # mean p in radius
        wave['t'+st+'_mean']=[]     # t mean in radius   
        wave['t'+st+'_min']=[]     # t min in radius   
        wave['pw'+st+'_mean']=[]    # mean p in radius
        wave['tw'+st+'_mean']=[]     # t mean in radius  
        wave['scales'+st]=[] 
        wave['pnb'+st]=[] 
                      
               
    cntmax=0 
    cntin=0      
       
    for f in files:
        print('Doing file: '+f)
        dic = xr.open_dataset(f)
        
     #   if (dic['time.hour'].values<15) or (dic['time.hour'].values>21):
          # print('smaller') 
     #      continue       
                
        outp=np.array(dic['p'].values.copy())                
        wave['porig'].append(outp)
        
                
        for strr in strarr:
            
            outt=np.array(dic['tc_lag'+strr].values.copy())
            outp[np.isnan(outp)]=-10**-5
            wave['torig'+strr].append(outt)    
                        
            outt[np.isnan(outt)]=150
            outt[outt>-40]=150
            grad=np.gradient(outt)
            outt[outt>-40]=-55
            o2=outt.copy()
            nok = np.where(abs(grad[0]) > 80)
            d=2
            i=nok[0]
            j=nok[1]    
        
            for ii,jj in zip(i,j):    
                kernel=o2[ii-d:ii+d+1, jj-d:jj+d+1]
              #  if not kernel.any():
                 #   continue
         #   else:    
                o2[ii-d:ii+d+1, jj-d:jj+d+1]=ndimage.gaussian_filter(kernel, 3, mode='nearest')
                  
            wav = util.waveletTP_localMax(o2, outp, 5) 
            o2[np.where(dic['pmask'].values==0)]=np.nan   
                        
            if strr=='0':
              outp[np.where(dic['pmask'].values==0)]=np.nan  
              wave['pw'].append(wav['p'])
              wave['tw'+strr].append(wav['t'])  
              #print(wav['scales'])              
              cntmax = cntmax+len(wav['z'])
            
            maxoutt = (outt == ndimage.minimum_filter(outt,5, mode='constant',cval=np.amax(outt)+1))
            maxoutt = maxoutt.astype(int)
            ypks,xpks=np.where((maxoutt==1) & (outt < -55))
           
            z = wav['z'] 
            y = wav['y']  
            x = wav['x']  
             
            for i in range(len(z)):  
                                                
                zz=z[i]
                xx=x[i]
                yy=y[i]                

                                                
                if dic['pmask'][yy,xx]==0:    # if maximum falls in region where no TRMM exists, continue                   
                    continue
                
                sc=wav['scales'][zz]                                   
                                    
                if strr=='0':                       
                    cntin = cntin+1
                                                  
                iscale = (np.ceil(wav['scales'][zz]/2./5.)).astype(int)
                
                tw = wav['t'][zz, :, :].copy()
                pw = wav['p'][zz, :, :].copy()   #copy??
                
                tw[np.isnan(tw)]=0
                pw[np.isnan(pw)]=0
                tw[np.where(dic['pmask'].values==0)]=np.nan            
                pw[np.where(dic['pmask'].values==0)]=np.nan     
                
#                pw[np.isnan(pw)]=1000
#
#                ax=plt.axes(projection=ccrs.PlateCarree())                
#                plt.contour(dic['lon'], dic['lat'], pw, levels=np.arange(500,1001,100), transform=ccrs.PlateCarree())                                     
#                plt.show()    
                           
                twmax=tw[yy,xx]     
                print(twmax)
               
                #Find all indices within the local circle of radius iscale...
                # ... Then average over those indices
                xloc1 = np.arange(xx-iscale,xx+iscale+1)
                yloc1 = np.arange(yy-iscale,yy+iscale+1)
                xloc,yloc = np.meshgrid(xloc1,yloc1)
                distloc = ( (xloc-xx)**2 + (yloc-yy)**2 ) ** .5

                indloc = (distloc <= iscale).nonzero()
                ycirc = indloc[0] - iscale + yy
                xcirc = indloc[1] - iscale + xx   
                
              #  print('pwshape',pw.shape[0], pw.shape[1] )
              #  print('twshape',tw.shape[0], tw.shape[1] )

                noky=np.where(ycirc>=pw.shape[0])   # if the circle is off the edge                               
                if noky[0].size>0:
                    ycirc=np.delete(ycirc,noky)
                    xcirc=np.delete(xcirc,noky)
                    
                nokx=np.where(xcirc>=pw.shape[1])                                  
                if nokx[0].size>0:
                    ycirc=np.delete(ycirc,nokx)
                    xcirc=np.delete(xcirc,nokx)    
                               
                tmean=np.nanmean(dic['tc_lag'+strr].values[ycirc, xcirc])
                pmean=np.nanmean(dic['p'].values[ycirc, xcirc])
                twmean=np.nanmean(tw[ycirc, xcirc])
                pwmean=np.nanmean(pw[ycirc, xcirc])
                pmax=np.nanmax(outp[ycirc, xcirc])
                pwmax=np.nanmax(pw[ycirc, xcirc])
                tmin=np.nanmin(outt[ycirc, xcirc])
                pnb=ycirc.size
              
                wave['tw'+strr+'_max'].append(twmax)            
                wave['pw'+strr+'_max'].append(pwmax)
                wave['tw'+strr+'_mean'].append(twmean)            
                wave['pw'+strr+'_mean'].append(pwmean)
                wave['p'+strr+'_max'].append(pmax)
                wave['p'+strr+'_mean'].append(pmean)
                wave['t'+strr+'_mean'].append(tmean)
                wave['t'+strr+'_min'].append(tmin)
                wave['scales'+strr].append(sc)
                wave['pnb'+strr].append(pnb)
               
                
                #just append all the variables into a dictionary now!                                                 
    
    for k in wave:
        if isinstance(wave[k][0], np.ndarray):
            continue
        print(k)
        wave[k]=np.array(wave[k])
    
    pkl.dump(wave, open('/users/global/cornkle/MCSfiles/save/MCS_wavelet_allyears_label.p', 'wb'))
    
    print('Saved!') 
    print('Found '+str(cntmax)+' maxima in '+str(len(files))+' systems.')
    print(str(cntin)+' maxima coincided with TRMM')
    
#readMCS_getWavelet_label()
     
    
    
def relatePintensities():
    
    wave = pkl.load( open ('/users/global/cornkle/MCSfiles/save/MCScorr_scales.p', 'rb'))
       
    arr=np.array([15,   16,   17,   18,   19,   20,   21,   22,   24,
         25,   27,   28,   30,   32,   34,   36,   38,   40,
         42,   45,   48,   50,   53,   57,   60], dtype=str)
    l=[]
    for k in arr:
        for kk in arr:
            l.append(np.mean(wave[k][kk]['corr']))         
                
    nl=np.array(l)
    nl_resh=nl.reshape(25,25)
    y = np.argmax(nl_resh, axis=0)  #temperature for P, index for getting the right T if you decide on a size for P
    x = np.argmax(nl_resh, axis=1) # precipitation for T, index for getting the right P if you decide on a size for T
    rise=range(25)
     
    PforT=arr[np.array(x)] # p for T
    TforP=arr[np.array(y)] # T for p 
    PforT=PforT.astype(np.str)
    TforP=TforP.astype(np.str)
    
    non = lambda s: s if s<0 else None
    mom = lambda s: max(0,s)
    
    dic={}
    for a in arr:    
        dic[a]={'p': [], 't' : []}

    for tpos, ppos in zip(arr, PforT):  

        print('Tscale: ', tpos)          
        print('Pscale: ', ppos)            
        
        tt=wave[tpos][ppos]['twavelet']
        pp=wave[tpos][ppos]['pwavelet']
        porig=wave[tpos][ppos]['porig']
        xx=wave[tpos][ppos]['x']
        yy=wave[tpos][ppos]['y']
        
        for t, p, po, x, y in zip(tt, pp, porig, xx, yy):
                                    
            ox, oy = x-5, y-5

            shift_p = np.zeros_like(p)
            shift_p[mom(oy):non(oy), mom(ox):non(ox)] = p[mom(-oy):non(-oy), mom(-ox):non(-ox)]    
            
            shift_po = np.zeros_like(po)
            shift_po[mom(oy):non(oy), mom(ox):non(ox)] = po[mom(-oy):non(-oy), mom(-ox):non(-ox)]                  
            ok = np.where(shift_p > 1)
            
            pout=shift_po[ok]
            tout=t[ok]
            
            dic[tpos]['p'].extend(pout)
            dic[tpos]['t'].extend(tout)                       

    pkl.dump(dic, open('/users/global/cornkle/MCSfiles/save/MCS_pintensity.p', 'wb'))
            
    
def plotShit1():
    
    wavelet = pkl.load( open ('/users/global/cornkle/MCSfiles/save/MCS_wavelet_allyears_perc.p', 'rb'))
    #plt.figure(1)
    scale=np.array(['10', '30'])
    
    t0=np.array(wavelet['torig0_min'])*(-1)
    t1=np.array(wavelet['torig1_min'])*(-1)
    t2=np.array(wavelet['torig2_min'])*(-1)
    t3=np.array(wavelet['torig3_min'])*(-1)
    tx=np.array(wavelet['torigx_min'])*(-1)
    p=np.array(wavelet['porig_max'])
    plt100=np.where(p<70)
  #  print(plt100[0:10])
    
    f = plt.figure()
    
   # m,b = np.polyfit(p,t1,1)
    r0=pearsonr(p[plt100],t0[plt100])
    r1=pearsonr(p[plt100],t1[plt100])
    r2=pearsonr(p[plt100],t2[plt100])
    r3=pearsonr(p[plt100],t3[plt100])
    rx=pearsonr(p[plt100],tx[plt100])
        
    title = "Tmin and Pmax (t0)"
    ax = f.add_subplot(361)
    ax.set_title(title,fontsize=12)
    
    ax.text(0.7, 0.8, 'r: '+str(np.round(r0[0], decimals=2)), transform=ax.transAxes)
    ax.scatter(t0[plt100], p[plt100])
    #plt.plot(p, m*p+b, '-')    
    ax.set_xlabel('T')
    ax.set_ylabel('P')
    
    title = "Tmin and Pmax (t-1)"
    ax = f.add_subplot(362)
    ax.set_title(title,fontsize=12)
    ax.text(0.7, 0.8, 'r: '+str(np.round(r1[0], decimals=2)), transform=ax.transAxes)
    ax.scatter(t1[plt100], p[plt100])
    ax.set_xlabel('T')
    ax.set_ylabel('P')    
    
    title = "Tmin and Pmax (t-2)"
    ax = f.add_subplot(363)
    ax.set_title(title,fontsize=12)
    ax.text(0.7, 0.8, 'r: '+str(np.round(r2[0], decimals=2)), transform=ax.transAxes)
    ax.scatter(t2[plt100], p[plt100])
    ax.set_xlabel('T')
    ax.set_ylabel('P')  
    
    title = "Tmin and Pmax (t-3)"
    ax = f.add_subplot(364)
    ax.set_title(title,fontsize=12)
    ax.text(0.7, 0.8, 'r: '+str(np.round(r3[0], decimals=2)), transform=ax.transAxes)
    ax.scatter(t3[plt100], p[plt100])
    ax.set_xlabel('T')
    ax.set_ylabel('P')  
    
    title = "Tmin and Pmax (t+2)"
    ax = f.add_subplot(365)
    ax.set_title(title,fontsize=12)
    ax.text(0.7, 0.8, 'r: '+str(np.round(rx[0], decimals=2)), transform=ax.transAxes)
    ax.scatter(tx[plt100], p[plt100])
    ax.set_xlabel('T')
    ax.set_ylabel('P') 
        
        
    for r in range(scale.size):
        sscale=scale[r]
        t1=np.array(wavelet[sscale]['twavelet0_max'])
        t2=np.array(wavelet[sscale]['twavelet2_max'])
        tx=np.array(wavelet[sscale]['twaveletx_max'])
        p=np.array(wavelet['porig_max'])
        
        r1=pearsonr(p[plt100],t1[plt100])
        r2=pearsonr(p[plt100],t2[plt100])
        rx=pearsonr(p[plt100],tx[plt100])
                                  
        title = "TWmax and Pmax (t0)"+sscale
        ax = f.add_subplot(3,6,r+7)
        ax.set_title(title,fontsize=12)
        ax.text(0.7, 0.8, 'r: '+str(np.round(r1[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(np.sqrt(t1[plt100]),p[plt100])
        ax.set_xlabel('T')
        ax.set_ylabel('P')
        
        title = "TWmax and Pmax (t-2)"+sscale
        ax = f.add_subplot(3,6,r+9)
        ax.set_title(title,fontsize=12)
        ax.text(0.7, 0.8, 'r: '+str(np.round(r2[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(np.sqrt(t2[plt100]), p[plt100])
        ax.set_xlabel('T')
        ax.set_ylabel('P')
        
        title = "TWmax and Pmax (t+2)"+sscale
        ax = f.add_subplot(3,6,r+11)
        ax.set_title(title,fontsize=12)
        ax.text(0.7, 0.8, 'r: '+str(np.round(rx[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(np.sqrt(tx[plt100]), p[plt100])
        ax.set_xlabel('T')
        ax.set_ylabel('P')
        
    for r in range(scale.size):
        
        sscale=scale[r]
        t1=np.array(wavelet[sscale]['twavelet0_max'])
        t2=np.array(wavelet[sscale]['twavelet2_max'])
        tx=np.array(wavelet[sscale]['twaveletx_max'])
        p=np.array(wavelet[sscale]['pwavelet_max'])
        
        r1=pearsonr(np.sqrt(p[plt100]),np.sqrt(t1[plt100]))
        r2=pearsonr(np.sqrt(p[plt100]),np.sqrt(t2[plt100]))
        rx=pearsonr(np.sqrt(p[plt100]),np.sqrt(tx[plt100]))
                                  
        title = "TWmax and PWmax (t0)"+sscale
        ax = f.add_subplot(3,6,r+13)
        ax.set_title(title,fontsize=12)
        ax.text(0.7, 0.8, 'r: '+str(np.round(r1[0], decimals=2)), transform=ax.transAxes)
        ax.set_xlim([-0.5,25])
        ax.set_ylim([-0.5,20])
        ax.scatter(np.sqrt(t1[plt100]), np.sqrt(p[plt100]))
        ax.set_xlabel('T')
        ax.set_ylabel('P')
        
        title = "TWmax and PWmax (t-2)"+sscale
        ax = f.add_subplot(3,6,r+15)
        ax.set_title(title,fontsize=12)
        ax.text(0.7, 0.8, 'r: '+str(np.round(r2[0], decimals=2)), transform=ax.transAxes)
        ax.set_xlim([-0.5,25])
        ax.set_ylim([-0.5,20])
        ax.scatter(np.sqrt(t2[plt100]), np.sqrt(p[plt100]))
        ax.set_xlabel('T')
        ax.set_ylabel('P')
        
        title = "TWmax and Pmax (t+2)"+sscale
        ax = f.add_subplot(3,6,r+17)
        ax.set_title(title,fontsize=12)
        ax.text(0.7, 0.8, 'r: '+str(np.round(rx[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(np.sqrt(tx[plt100]), p[plt100])
        ax.set_xlabel('T')
        ax.set_ylabel('P')
            
        
 #   plt.tight_layout()
    plt.show() 
    
def plotShit_label():
    
    wavelet = pkl.load( open ('/users/global/cornkle/MCSfiles/save/MCS_wavelet_allyears_label.p', 'rb'))
    #plt.figure(1)
    
    strg=['0', '1', '2', '3', 'x']
    f = plt.figure()
    k=361
    cnt=0
    for strr in strg:    
    
        t=np.array(wavelet['t'+strr+'_min'])#*(-1)   
        p=np.array(wavelet['p'+strr+'_max'])
        plt100=np.where(p<70)
        #  print(plt100[0:10])

        # m,b = np.polyfit(p,t1,1)
        r=pearsonr(p[plt100],t[plt100])
            
        title = "Tmin and Pmax (t"+strr+")"
        ax = f.add_subplot(361+cnt)
        ax.set_title(title,fontsize=12)
    
        ax.text(0.7, 0.8, 'r: '+str(np.round(r[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(t[plt100], p[plt100])
        #plt.plot(p, m*p+b, '-')    
        ax.set_xlabel('T')
        ax.set_ylabel('P')
        cnt=cnt+1
    cnt=0     
    for strr in strg:          
        pos=np.where(wavelet['scales'+strr]==18)
        t=np.sqrt(np.array(wavelet['tw'+strr+'_max'])[pos] ) 
        p=np.array(wavelet['p'+strr+'_max'][pos])     
        
       # t=t[np.isfinite(t)]
      #  p=p[np.isfinite(p)]

        # m,b = np.polyfit(p,t1,1)
        r=pearsonr(p,t)
        print(strr)
        print(cnt)
            
        title = "sqrt(TWmax) and Pmax (t"+strr+")"
        ax = f.add_subplot(3,6,7+cnt)
        ax.set_title(title,fontsize=12)
    
        ax.text(0.7, 0.8, 'r: '+str(np.round(r[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(t, p)
        #plt.plot(p, m*p+b, '-')    
        ax.set_xlabel('TWmax 18km')
        ax.set_ylabel('Pmax 18km')
        cnt=cnt+1
    
    cnt=0    
    for strr in strg:    
        pos=np.where(wavelet['scales'+strr]==25)
        t=np.sqrt(np.array(wavelet['tw'+strr+'_max'])[pos] ) 
        p=np.array(wavelet['p'+strr+'_max'][pos])        

        # m,b = np.polyfit(p,t1,1)
        r=pearsonr(p,t)
            
        title = "sqrt(TWmax) and Pmax (t"+strr+")"
        ax = f.add_subplot(3,6,13+cnt)
        ax.set_title(title,fontsize=12)
    
        ax.text(0.7, 0.8, 'r: '+str(np.round(r[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(t, p)
        #plt.plot(p, m*p+b, '-')    
        ax.set_xlabel('TWmax 25km')
        ax.set_ylabel('Pmax 25km')
        
       # tt=np.array(wavelet['t'+strr+'_min'])*(-1)   
     #   pp=np.array(wavelet['p'+strr+'_max'])
        cnt=cnt+1
#        
            
        
 #   plt.tight_layout()
    plt.show() 
#plotShit_label()    
#      