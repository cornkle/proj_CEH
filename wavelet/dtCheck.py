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
import pickle as pkl
from skimage.measure import compare_ssim
from scipy.stats.stats import pearsonr
from scipy import ndimage

HOD=range(24)   # hours of day
YRANGE=range(2004,2014)


def tshift():
    
    cdic = {'corr' : [], 'corrs' : [] ,'corrp' : [] ,  'dt': [], 'date' : [], 'id' : [], 'pmatch' : [], 'cmax' : [],
            'dtmax' : [], 'tdmax' : [], 'mdmax' : [],  'indmax' : [] ,  'idate' : [], 'pears' : [],'shift' : [], 'compare' : []}
    
    trmm_folder= "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder='/users/global/cornkle/data/OBS/meteosat_SA15'
    
    # make a salem grid
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')

    t=re.trmm(trmm_folder, yrange=range(2006, 2011), area=[-10, 10, 10, 20]) 
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
          
       goodinds = u[n>500]  # all blobs with more than 2500 pixels
       print(goodinds)
       if not sum(goodinds) > 0:
            continue
     
       for gi in goodinds:
              if gi == 0:
                  continue
              
              inds = np.where(labels == gi)
                            
              # cut a box for every single blob from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
              latmax, latmin = mdic['lat'][inds].max() , mdic['lat'][inds].min()
              lonmax, lonmin = mdic['lon'][inds].max() , mdic['lon'][inds].min()
              mmeans=np.percentile(mdic['t'][inds], 90)
              td = t.getDData(_y, _m, _d, _h, _mi, cut=[latmin-0.2, latmax+0.2])
              
              #ensure minimum trmm rainfall in area
              if np.size(np.nonzero(td['p']))< 50:
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
              ndate = date + dt.timedelta(minutes=int(dt2)  )     
          #    print('Date3', ndate)  
              ml2=m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute, llbox=[lonmin-0.3, lonmax+0.3, latmin-0.25, latmax+0.25])     
              if not ml2:
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
              outl = np.full_like(dummy, -150)
              xl, yl = grid.transform(lon1[inds], lat1[inds], crs=salem.wgs84, nearest=True, maskout=True)
              
              tpair = (xt+yt)*(xt+yt+1)/2+yt
              bpair = (xl.compressed()+yl.compressed())*(xl.compressed()+yl.compressed()+1)/2+yl.compressed()     
              iinter=np.in1d(tpair, bpair)    
             
              if sum(iinter) < 200:  
                continue
              print('Lag0')
              outl[yl.compressed(),xl.compressed()] = dummy[yl.compressed(), xl.compressed()]
        
              # Interpolate using delaunay triangularization 
              outt = griddata(tpoints, td['p'].flatten(), inter, method='linear')
              outt = outt.reshape((grid.ny, grid.nx)) 
              
              outt=outt[1:-1, 4:-4]
              outl=outl[1:-1, 4:-4]
                                        
              tmask = np.isfinite(outt)
              ttmask = np.isnan(outt)
              mmask = np.isfinite(outl)
              mask2 = np.isfinite(outl[tmask])

              if sum(mask2.flatten()) < sum(mmask.flatten())*0.3:
                  continue   
              
              print('Hit:', gi)
              
              # zero lag
              outt[np.isnan(outt)]=-10**-5
              outl[np.isnan(outl)]=-150     
              
              grad=np.gradient(outl)
              nok = np.where(abs(grad[1]) > 40)
                          
              
              outl[outl<-100]=mmeans
              
              d=2
              i=nok[0]
              j=nok[1]
              for ii,jj in zip(i,j):
    
                    kernel=outl[ii-d:ii+d+1, jj-d:jj+d+1]
                    if not kernel.any():
                        outl[ii,jj]=mmeans
                    else:    
                        outl[ii-d:ii+d+1, jj-d:jj+d+1]=ndimage.gaussian_filter(kernel, 3)

              pos=10
              dic = util.waveletTP(outl, outt, 5)
              
              tt=dic['t'][pos,:,:]
              tt[ttmask]=0
              pp=dic['p'][pos,:, :]
              corr00 = max(match_template(tt, pp[5:-5, 5:-5]).flatten())
              corr0 = compare_ssim(tt, pp)
              corr000=pearsonr(tt.flatten(), pp.flatten())
              
              cdic['corr'].append(corr0)
              cdic['corrs'].append(corr00)
              cdic['corrp'].append(corr000)
              cdic['dt'].append(dt0)
              cdic['date'].append(date)
              cdic['pmatch'].append(sum(mask2.flatten()))
              
              
              # lag -1
             
               # Interpolate using delaunay triangularization 
              dummy = griddata(mpoints, ml1['t'].flatten(), inter, method='linear')
              dummy = dummy.reshape((grid.ny, grid.nx))
              outl = np.full_like(dummy, mmeans)
              print('Lag1')
              outl[yl.compressed(),xl.compressed()] = dummy[yl.compressed(), xl.compressed()]
              outl=outl[1:-1, 4:-4]
              outl[np.isnan(outl)]=mmeans
              
              for ii,jj in zip(i,j):
    
                    kernel=outl[ii-d:ii+d+1, jj-d:jj+d+1]
                    if not kernel.any():
                        outl[ii,jj]=mmeans
                    else:    
                        outl[ii-d:ii+d+1, jj-d:jj+d+1]=ndimage.gaussian_filter(kernel, 3)
                                          
              dic = util.waveletTP(outl, outt, 5)              
              
              tt=dic['t'][pos,:,:]
              tt[ttmask]=0
              pp=dic['p'][pos,:, :]
              corr11 = max(match_template(tt, pp[5:-5, 5:-5]).flatten())
              corr1 = compare_ssim(tt, pp)
              corr111=pearsonr(tt.flatten(), pp.flatten())
              
              cdic['corr'].append(corr1)
              cdic['corrs'].append(corr11)
              cdic['corrp'].append(corr111)
              cdic['dt'].append(dt1)
              cdic['date'].append(date)
              cdic['pmatch'].append(sum(mask2.flatten()))
              
              # lag -2
              
               # Interpolate using delaunay triangularization 
              dummy = griddata(mpoints, ml2['t'].flatten(), inter, method='linear')
              dummy = dummy.reshape((grid.ny, grid.nx))   
              outl = np.full_like(dummy, mmeans)
              print('Lag2')
              outl[yl.compressed(),xl.compressed()] = dummy[yl.compressed(), xl.compressed()]
              outl=outl[1:-1, 4:-4]
              outl[np.isnan(outl)]=mmeans
              
              for ii,jj in zip(i,j):
    
                    kernel=outl[ii-d:ii+d+1, jj-d:jj+d+1]
                    if not kernel.any():
                        outl[ii,jj]=mmeans
                    else:    
                        outl[ii-d:ii+d+1, jj-d:jj+d+1]=ndimage.gaussian_filter(kernel, 3)                             
              dic = util.waveletTP(outl, outt, 5)              
              
              tt=dic['t'][pos,:,:]
              tt[ttmask]=0
              pp=dic['p'][pos,:, :]
              corr22 = max(match_template(tt, pp[5:-5, 5:-5]).flatten())
              corr2 = compare_ssim(tt, pp)
              corr222=pearsonr(tt.flatten(), pp.flatten())              
              
              cdic['corr'].append(corr2)
              cdic['corrs'].append(corr22)
              cdic['corrp'].append(corr222)
              cdic['dt'].append(dt2)
              cdic['date'].append(date)
              cdic['pmatch'].append(sum(mask2.flatten()))                          
              
              compare = np.array([corr0, corr1, corr2])
              dtt = [dt0, dt1, dt2]
              sshift = np.array([corr00, corr11, corr22])
              pears = [corr000, corr111, corr222]
              
              
              dtmax = dtt[np.argmax(sshift)]
              cmax = sshift.max()
              
              cdic['cmax'].append(cmax)
              cdic['dtmax'].append(dtmax)
              cdic['tdmax'].append(date)
              cdic['mdmax'].append(date+dt.timedelta(minutes=int(dtmax)))
              cdic['indmax'].append(gi)
              cdic['idate'].append(date + dt.timedelta(minutes=int(dt0) ))
              cdic['pears'].append(pears[np.argmax(sshift)])
              cdic['shift'].append(sshift[np.argmax(sshift)])
              cdic['compare'].append(compare[np.argmax(sshift)])
              
              cnt=cnt+1
          
    pkl.dump(cdic, open('/users/global/cornkle/timeshift2.p', 'wb'))
    print('Saved '+str(cnt)+' MCSs in timeshift file')             
            
           