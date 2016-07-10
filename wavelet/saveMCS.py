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
              dummy2 = griddata(mpoints, ml1['t'].flatten(), inter, method='linear')
              dummy2 = dummy2.reshape((grid.ny, grid.nx))
              outl2 = np.full_like(dummy2, np.nan)
              print('Lag1')
              outl2[yl.compressed(),xl.compressed()] = dummy2[yl.compressed(), xl.compressed()]              
                                                               
              
              da = xr.Dataset({'p' : (['x', 'y'], outt),  
                                 't_lag0'  : (['x', 'y'], dummy),
                                 'tc_lag0'  : (['x', 'y'], outl),
                                 't_lag1'  : (['x', 'y'], dummy2),
                                 'tc_lag1'  : (['x', 'y'], outl2),
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
             
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles')
    
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
            wave[a][b]=[]
                    
 #   return wave   
    
    for f in files:
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
                    wave[arr[pos]][arr[ppos]].append(np.max(corr1))
                else:     
                    wave[arr[pos]][arr[ppos]].append(np.max(corr2))
    return wave        
    
    pkl.dump(wave, open('/users/global/cornkle/MCSfiles/save/MCScorr_scales.p', 'wb'))
    
    print('Saved!')
    
    
def relatePintensities():
    
    myDicts = pkl.load( open ('/users/global/cornkle/MCSfiles/save/MCScorr.p', 'rb'))
    arr=np.array([15,   16,   17,   18,   19,   20,   21,   22,   24,
         25,   27,   28,   30,   32,   34,   36,   38,   40,
         42,   45,   48,   50,   53,   57,   60], dtype=str)
    l=[]
    for k in arr:
        for kk in arr:
            l.append(np.mean(myDicts[k][kk]))         
                
    nl=np.array(l)
    nl_resh=nl.reshape(25,25)
    y = np.argmax(nl_resh, axis=0)  #temperature
    x = np.argmax(nl_resh, axis=1) # precipitation
    rise=range(25)
 