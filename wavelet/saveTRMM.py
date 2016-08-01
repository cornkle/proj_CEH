# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 16:53:55 2016

@author: cornkle
"""

import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
import datetime as dt
from eod import read_eod as re
import xarray as xr
import os

HOD=range(24)   # hours of day
YRANGE=range(2004,2014)

#BOX=[XL, XR, YL, YU]
def saveTRMM():
   
    trmm_folder= "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"   
    box = [-10, 10, 10, 20]
    # make grid
    #define projection
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')
    #get lower left x,y fr 10W, 4N
    x, y = pyproj.transform(salem.wgs84, proj, [box[0], box[1]], [box[2], box[3]] )  
    dx = 5000    # 5km grid
    nx, r = divmod(x[1] - x[0], dx)
    ny, r = divmod(y[1] - y[0], dx)
    #make salem grid
    grid = salem.Grid(nxny=(nx,ny), dxdy=(5000,5000), ll_corner=(x[0], y[0]), proj=proj)
     
    xi, yi = grid.ij_coordinates
    lon, lat = grid.ll_coordinates

    t=re.trmm(trmm_folder, yrange=YRANGE, area=[box[0], box[1], box[2], box[3]])    
    
    cnt=0    
    
    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered      
    for _y, _m, _d, _h, _mi in zip(t.dates.y, t.dates.m, t.dates.d, t.dates.h, t.dates.mi):
       
       #dummy = np.empty((ny,nx))*-100#np.NAN
       td=t.getDData(_y, _m, _d, _h, _mi, cut=[box[2], box[3]])
                                          
       date=dt.datetime(_y, _m, _d, _h, _mi)
       print(date)
       
                                                 #ensure minimum trmm rainfall in area
       if len(np.where(td['p']>0)[0]) < 100:      #  at least 100 pixel with rainfall 
          print('Kickout: TRMM min pixel = 100')
          continue                                                                             

       # Transform lons, lats to grid              
       xt, yt = grid.transform(td['lon'].flatten(), td['lat'].flatten(), crs=salem.wgs84)
        
       # Convert for griddata input               
       tpoints = np.array((yt, xt)).T
       inter = np.array((np.ravel(yi), np.ravel(xi))).T                    
        
       # Interpolate using delaunay triangularization 
       dummyt = griddata(tpoints, td['p'].flatten(), inter, method='linear')
       outt = dummyt.reshape((grid.ny, grid.nx))        

       
       for nb in range(5):
           boole=np.isnan(outt)
           outt[boole]=-1000
           grad=np.gradient(outt)
           outt[boole]=np.nan
           outt[abs(grad[1])>300]=np.nan
           outt[abs(grad[0])>300]=np.nan  

       attrs = {'units': 'hours since 2000-01-01'}                                                                        
              
       da = xr.Dataset({'p' : (['x', 'y', 'time'], outt),  
                        'time' : ('time', date)       }, 
                            coords= {'lon' : (['x', 'y'], lon),
                                     'lat' : (['x', 'y'], lat),                                    
                                      })
       da.close()  
       savefile = '/users/global/cornkle/TRMMfiles/'+date.strftime('%Y-%m-%d_%H:%M:%S')+'.nc'
       try:
           os.remove(savefile)
       except OSError:
           pass
       da.to_netcdf(path=savefile, mode='w')
       print('Saved '+savefile)
                       
       cnt=cnt+1
             
    print('Saved '+str(cnt)+' TRMM swaths as netcdf.')   