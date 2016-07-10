# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 16:26:38 2016

@author: cornkle
"""

import salem
import cleo
import pyproj
import numpy as np
from scipy.interpolate import griddata
from eod import msg
import xarray as xray
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs

def cutout_test_MCS():
    
    # 'nb : dy1,dy2, dx1, dx2, tdy1, td2, tdx1, tdx2'
    
    ind_dic={ '845' : [1,290, 20, 360], '730' : [50, 350, 0, 280], '437' : [50, 400,20, 450], '490' : [0, 450, 20, 450]}    
   
    
    # make a salem grid
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')
    
   # sav_mcs={'845' : {'lat' : [], 'lon' : [], 'dmin' : [], 'date' : [], 'trmm' : [], 'msg' : [] }}
  #  sav_mcs={ }
    HOD=range(24)   # hours of day
    YRANGE=range(2004,2014)

    tpath= "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    tf=msg.extract_TRMMfile(tpath, hod=HOD, yrange=YRANGE, mtresh=1)       
    print(tf)
    tnb=np.array(range(len(tf['fpath'])))
    maxs=[]
    for tp in tnb:
        dmsg, dtrmm=msg.quickreadTrmmMSG(tf, nb=tp)
        td=dtrmm['pcp'].mean().max()
        maxs.append(td)
    sort_max, sort_nbs = zip(*sorted(zip(maxs, tnb), reverse=True))
    
    print(tf['date'].getStr(['845']))
    return
    
        # Transform lons, lats to grid
    print(ind_dic.keys())
    for nb in ind_dic.keys():   
                       
        nbn=int(nb)
        print(nbn)
        mi=tf['tmins'][nbn]
        date=tf['date'].getStr([nbn])
        print(date)
        dmsg, dtrmm=msg.quickreadTrmmMSG(tf, nb=nbn)

        md=dmsg['t']
        td=dtrmm['pcp']
        
        dy1=ind_dic[nb][0]
        dy2=ind_dic[nb][1]
        dx1=ind_dic[nb][2]
        dx2=ind_dic[nb][3]
        print(dy1,dy2,dx1,dx2)
        
        mlat=dmsg['lats'][dy1:dy2, dx1:dx2]
        mlon=dmsg['lons'][dy1:dy2, dx1:dx2]
        tlat=dtrmm['lats']#[tdy1:tdy2, tdx1:tdx2]
        tlon=dtrmm['lons']#[tdy1:tdy2, tdx1:tdx2]                                

        mx, my = pyproj.transform(salem.wgs84, proj, mlon, mlat)
        tx, ty = pyproj.transform(salem.wgs84, proj, tlon, tlat)
        
        ax=mx.flatten().tolist()
        ay=my.flatten().tolist()
        
        # take the min and max
        xmax, xmin = max(ax), min(ax)
        ymax, ymin = max(ay), min(ay) 

        print(xmax, xmin, ymax, ymin)
        
        #Count the number of pixels
        dx = 5000
        nx, r = divmod(xmax - xmin, dx)
        ny, r = divmod(ymax - ymin, dx)
        # Here one could add + 1 to be sure that the last pixel is always included
        grid = salem.Grid(nxny=(nx, ny), dxdy=(dx, dx), ll_corner=(xmin, ymin), proj=proj)     
        cm = cleo.Map(grid)
    #    cm.visualize()
        xi, yi = grid.ij_coordinates
        lon, lat = grid.ll_coordinates
        
        print('Grid coord: '+str(min(lon.flatten()))+str(max(lon.flatten()))+str(min(lat.flatten()))+str(max(lat.flatten())))
        
        # Transform lons, lats to grid
        xm, ym = grid.transform(mlon.flatten(), mlat.flatten(), crs=salem.wgs84)
        xt, yt = grid.transform(tlon.flatten(), tlat.flatten(), crs=salem.wgs84)
        
        # Convert for griddata input 
        mpoints = np.array((ym, xm)).T
        tpoints = np.array((yt, xt)).T
        inter = np.array((np.ravel(yi), np.ravel(xi))).T
        
        # Interpolate using delaunay triangularization 
        outm = griddata(mpoints, md[dy1:dy2, dx1:dx2].flatten(), inter, method='linear')
        outm = outm.reshape((grid.ny, grid.nx))
        
        # Interpolate using delaunay triangularization 
        outt = griddata(tpoints, td.flatten(), inter, method='linear')
        outt = outt.reshape((grid.ny, grid.nx))
        
        dx=10
        outm=outm[dx:-dx, dx:-dx]
        outt=outt[dx:-dx, dx:-dx]
        lon=lon[dx:-dx, dx:-dx]
        lat=lat[dx:-dx, dx:-dx]
        outt[np.isnan(outt)]=-10**-5
        outm[np.isnan(outm)]=30        
            
        da = xray.Dataset({'trmm' : (['x', 'y'], outt), 
                           'msg'  : (['x', 'y'], outm)}, 
                            coords= {'lon' : (['x', 'y'], lon),
                                     'lat' : (['x', 'y'], lat),
                                     'time': date })
        da.attrs['mins_trmm']=mi
        da.close()  
        savefile = '/users/global/cornkle/wtest/'+nb+'_mt_wavelet_test.nc'
        try:
            os.remove(savefile)
        except OSError:
            pass
        da.to_netcdf(path=savefile, mode='w')
        

    print('Saved '+savefile)
   
        
    
def plot_test_MCS(): 
   HOD=range(24)   # hours of day
   YRANGE=range(2004,2014)

   tpath= "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
   tf=msg.extract_TRMMfile(tpath, hod=HOD, yrange=YRANGE, mtresh=1)
   
   tnb=np.array(range(len(tf['fpath'])))
   maxs=[]
   for tp in tnb:
    dmsg, dtrmm=msg.quickreadTrmmMSG(tf, nb=tp)
    td=dtrmm['pcp'].mean().max()
    maxs.append(td)
   sort_max, sort_nbs = zip(*sorted(zip(maxs, tnb), reverse=True))
   top5=sort_nbs[0:5]
   NB=top5[3]
   mi=tf['tmins'][NB]
   dmsg, dtrmm=msg.quickreadTrmmMSG(tf, nb=NB)
   mlat=dmsg['lats']
   mlon=dmsg['lons']
   md=dmsg['t']
   tlat=dtrmm['lats']
   tlon=dtrmm['lons']
   td=dtrmm['pcp']
   
   dy1=0
   dy2=550
   dx1=0
   dx2=500
   tdy1=0
   tdy2=-1

   ax = plt.axes(projection=ccrs.PlateCarree())
   ax.coastlines()
   plt.contourf(tlon[tdy1:tdy2, :], tlat[tdy1:tdy2, :], td[tdy1:tdy2, :],levels=np.arange(0,100, 0.5), transform=ccrs.PlateCarree(), cmap='Pastel1')
   plt.contourf(mlon[dy1:dy2, dx1:dx2], mlat[dy1:dy2, dx1:dx2], md[dy1:dy2, dx1:dx2], levels=np.arange(-80,-30,1), transform=ccrs.PlateCarree())
   cbar=plt.colorbar()
   cbar.set_label('T (degC)', rotation=270, labelpad=+11)
   plt.contourf(tlon, tlat, td,levels=np.arange(4,70, 0.5), transform=ccrs.PlateCarree())

   cbar=plt.colorbar()
   cbar.set_label('PCP (mm h-1)', rotation=270, labelpad=+11)
   plt.savefig('/users/global/cornkle/'+str(NB)+'_msgtd0_p_'+str(mi)+'.pdf')
    
   ax = plt.axes(projection=ccrs.PlateCarree())
   plt.contourf(tlon[tdy1:tdy2, :], tlat[tdy1:tdy2, :], td[tdy1:tdy2, :],levels=np.arange(0,100, 0.5), transform=ccrs.PlateCarree(), cmap='Pastel1')
   plt.contourf(mlon[dy1:dy2, dx1:dx2], mlat[dy1:dy2, dx1:dx2], md[dy1:dy2, dx1:dx2], levels=np.arange(-80,-30,1), transform=ccrs.PlateCarree())
   cbar=plt.colorbar()
   cbar.set_label('T (degC)', rotation=270, labelpad=+11)
   ax.coastlines()
   plt.savefig('/users/global/cornkle/'+str(NB)+'_msgtd0.pdf')
   
   