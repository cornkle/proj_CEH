# -*- coding: utf-8 -*-
"""
Created on Tue May 31 11:46:10 2016

@author: cornkle
"""

import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage.measurements import label 
import datetime as dt
from eod import read_eod as re
from wavelet import util
from skimage.feature import match_template

HOD=range(24)   # hours of day
YRANGE=range(2004,2014)


def test_msg_tshift():
    
    cdic = {'corr' : [] ,  'dt': []}
    
    trmm_folder= "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder='/users/global/cornkle/data/OBS/meteosat_SA15'
    
    # make a salem grid
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')

    t=re.trmm(trmm_folder, yrange=range(2008, 2009), area=[-10, 10, 10, 20]) 
    m=re.msg(msg_folder)



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
            
       mdic['t'][mdic['t']>-20]=0
       labels, numL = label(mdic['t'])
           
       u , inv = np.unique(labels, return_inverse=True)
       n = np.bincount(inv)
          
       goodinds = u[n>2500]  # all blobs with more than 2500 pixels
       #print(goodinds)
       if not goodinds.any():
            continue
           
       for gi in goodinds:
              if gi == 0:
                  continue
              
              inds = np.where(labels == gi)              
              # cut a box for every single blob from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
              latmax, latmin = mdic['lat'][inds].max() , mdic['lat'][inds].min()
              lonmax, lonmin = mdic['lon'][inds].max() , mdic['lon'][inds].min()
              td = t.getDData(_y, _m, _d, _h, _mi, cut=[latmin-0.2, latmax+0.2])
              
              dt0=dm[ind] 
              ndate = date + dt.timedelta(minutes=int(dt0) )                                           
              ml0=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[lonmin-0.3, lonmax+0.3, latmin-0.25, latmax+0.25])        
       
              dt1=dm[ind]-15
              ndate = date + dt.timedelta(minutes=int(dt1))       
              ml1=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[lonmin-0.3, lonmax+0.3, latmin-0.25, latmax+0.25])
       
              dt2=dm[ind]-30      
              ndate = date + dt.timedelta(minutes=int(dt2)  )     
              ml2=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[lonmin-0.3, lonmax+0.3, latmin-0.25, latmax+0.25])     
              
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
              outm = griddata(mpoints, ml0['t'].flatten(), inter, method='linear')
              outm = outm.reshape((grid.ny, grid.nx))
        
              # Interpolate using delaunay triangularization 
              outt = griddata(tpoints, td['p'].flatten(), inter, method='linear')
              outt = outt.reshape((grid.ny, grid.nx)) 
              
              outt=outt[1:-1, 4:-4]
              outm=outm[1:-1, 4:-4]
                                        
              tmask = np.isfinite(outt)
              mmask = np.isfinite(outm)
              mask2 = np.isfinite(outm[tmask])

              if sum(mask2.flatten()) < sum(mmask.flatten())*0.3:
                  continue         
              
              # zero lag
              outt[np.isnan(outt)]=-10**-5
              outm[np.isnan(outm)]=30     
              dic = util.waveletTP(outm, outt, 5)
              
              tt=dic['t'][20,:,:]
              pp=dic['p'][20,:, :]
              corr = min(match_template(tt, pp[10:-10, 10:-10]).flatten())
              
              cdic['corr'].append(corr)
              cdic['dt'].append(dt0)
              
              
              # lag -1
             
               # Interpolate using delaunay triangularization 
              outm1 = griddata(mpoints, ml1['t'].flatten(), inter, method='linear')
              outm1 = outm1.reshape((grid.ny, grid.nx))
              outm1=outm1[1:-1, 4:-4]
              outm1[np.isnan(outm1)]=30     
              dic = util.waveletTP(outm1, outt, 5)              
              
              tt=dic['t'][20,:,:]
              pp=dic['p'][20,:, :]
              corr = min(match_template(tt, pp[10:-10, 10:-10]).flatten())
              
              cdic['corr'].append(corr)
              cdic['dt'].append(dt1)
              
              # lag -2
              
               # Interpolate using delaunay triangularization 
              outm2 = griddata(mpoints, ml2['t'].flatten(), inter, method='linear')
              outm2 = outm2.reshape((grid.ny, grid.nx))   
              outm2=outm2[1:-1, 4:-4]
              outm2[np.isnan(outm2)]=30                                 
              dic = util.waveletTP(outm2, outt, 5)              
              
              tt=dic['t'][20,:,:]
              pp=dic['p'][20,:, :]
              corr = max(match_template(tt, pp[10:-10, 10:-10]).flatten())
              
              cdic['corr'].append(corr)
              cdic['dt'].append(dt2)
          
    return cdic              
            
           