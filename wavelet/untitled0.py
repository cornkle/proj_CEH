# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 15:33:56 2016

@author: cornkle
"""

import numpy as np
import pickle as pkl
from netCDF4 import Dataset as nc

def minitest():
    
    tt=np.array(range(24), dtype=float)
    tt=np.insert(tt, 0, np.nan, axis=0)

    dic={'bla':tt}    
    
    print(tt[0:5])
    
    nf=nc('/users/global/cornkle/mt_wavelet_test_v3.5.nc', 'w', format='NETCDF4_CLASSIC')
    
    nf.createDimension('we', len(tt))
    dat=nf.createVariable('test', 'f4', 'we')
    dat[:]=tt
    
    nf.close()
    
    
    
    
def pkl_load():
   # myDicts = pkl.load( open ('/users/global/cornkle/mt_wavelet_test_v3.5.p', 'rb'))
    myDicts=np.load('/users/global/cornkle/mt_wavelet_test_v3.5.npy')
    return myDicts