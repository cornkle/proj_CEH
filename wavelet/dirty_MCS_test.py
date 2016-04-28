# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 16:26:38 2016

@author: cornkle
"""

import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
from eod import msg
import xray

def cutout_test_MCS():
    
    # 'nb : dy1,dy2, dx1, dx2, tdy1, td2, tdx1, tdx2'
    
    ind_dic={ '845' : [0,300, 30, 400, 80, 400, 0, -1], }
    
   
    
    # make a salem grid
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')
    
    sav_mcs={'845' : {'lat' : [], 'lon' : [], 'dmin' : [], 'date' : [], 'trmm' : [], 'msg' : [] }}
  #  sav_mcs={ }
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
    
        # Transform lons, lats to grid
    print(ind_dic.keys())
    for nb in ind_dic.keys():   
        nbn=int(nb)
        print(nbn)
        mi=tf['tmins'][nbn]
        date=tf['date'].get_str([nbn])
        dmsg, dtrmm=msg.quickreadTrmmMSG(tf, nb=nbn)

        md=dmsg['t']
        td=dtrmm['pcp']
        
        dy1=ind_dic[nb][0]
        dy2=ind_dic[nb][1]
        dx1=ind_dic[nb][2]
        dx2=ind_dic[nb][3]
        tdy1=ind_dic[nb][4]
        tdy2=ind_dic[nb][5]
        tdx1=ind_dic[nb][6]
        tdx2=ind_dic[nb][7]
        
        mlat=dmsg['lats'][dy1:dy2, dx1:dx2]
        mlon=dmsg['lons'][dy1:dy2, dx1:dx2]
        tlat=dtrmm['lats'][tdy1:tdy2, tdx1:tdx2]
        tlon=dtrmm['lons'][tdy1:tdy2, tdx1:tdx2]
        
   #     print(min(mlat.flatten()), min(mlon.flatten()))
   #     print(min(tlat.flatten()), min(mlon.flatten()))
                                
        # Transform lon, lats to the mercator projection
        mx, my = pyproj.transform(salem.wgs84, proj, mlon, mlat)
        tx, ty = pyproj.transform(salem.wgs84, proj, tlon, tlat)
        
        ax=mx.flatten().tolist()
        ax.extend(tx.flatten().tolist())
        ay=my.flatten().tolist()
        ay.extend(ty.flatten().tolist())
        
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
        xi, yi = grid.ij_coordinates
        lon, lat = grid.ll_coordinates
        
        print('Grid coord: '+str(min(lon.flatten()))+str(max(lon.flatten()))+str(min(lat.flatten()))+str(max(lat.flatten())))
        
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
        outt = griddata(tpoints, td[tdy1:tdy2, tdx1:tdx2].flatten(), inter, method='linear')
        outt = outt.reshape((grid.ny, grid.nx))

        sav_mcs[nb]['lat']=lat
        sav_mcs[nb]['lon']=lon
        sav_mcs[nb]['dmin']=mi 
        sav_mcs[nb]['date']=date
        sav_mcs[nb]['trmm']=outt
        sav_mcs[nb]['msg']=outm
        
        da = xray.Dataset({'trmm' : (['x', 'y'], outt), 
                           'msg'  : (['x', 'y'], outm)}, 
                            coords= {'lon' : (['x', 'y'], lon),
                                     'lat' : (['x', 'y'], lat),
                                     'time': date })
        da.attrs['mins_trmm']=mi
          
        savefile = '/users/global/cornkle/wtest/'+nb+'_mt_wavelet_test.nc'
        da.to_netcdf(path=savefile, mode='w')

    print('Saved '+savefile)
   
        
    
    