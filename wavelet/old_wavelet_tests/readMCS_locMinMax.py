# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:20:36 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
from utils import u_arrays as ua
import matplotlib.pyplot as plt
from scipy import ndimage
import pickle as pkl
from scipy.stats import pearsonr
    
def plotShit_label_noW():
    
    wavelet = pkl.load( open ('/users/global/cornkle/MCSfiles/save/MCS_allyears_label.p', 'rb'))
    #plt.figure(1)
    
    strg=['0', '1', '2', '3', 'x']
    f = plt.figure()
    k=361
    cnt=0
    for strr in strg:    
    
        t=np.array(wavelet['t'+strr+'_min'])#*(-1)   
        p=np.array(wavelet['p'+strr+'_max'])
        plt100=np.where(t<-40)
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
                  
        t=np.array(wavelet['t'+strr+'_mean'])#*(-1)   
        p=np.array(wavelet['p'+strr+'_max'])
        plt100=np.where(t<-40)
        #  print(plt100[0:10])

        # m,b = np.polyfit(p,t1,1)
        r=pearsonr(p[plt100],t[plt100])
            
        title = "Tmean and Pmax (t"+strr+")"
        ax = f.add_subplot(3,6,7+cnt)
        ax.set_title(title,fontsize=12)
    
        ax.text(0.7, 0.8, 'r: '+str(np.round(r[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(t[plt100], p[plt100])
        #plt.plot(p, m*p+b, '-')    
        ax.set_xlabel('T')
        ax.set_ylabel('P')
        cnt=cnt+1
         
        
            
        
 #   plt.tight_layout()
    plt.show()          