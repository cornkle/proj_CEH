#!/bin/env python

import numpy as np

llFile = "/users/global/cmt/msg/tropWA/VERA_MSG_grid.gra"

llShape = (222,1384)
llMDI = np.float32(13.5)
ll = np.fromfile(llFile,dtype=llMDI.dtype)
lon = ll[0:222*1384]
lat = ll[222*1384:]
lat.shape = llShape
lon.shape = llShape

llsavefile = "/users/global/danbel/msg/tropWA/latlon"
np.savez(llsavefile,lon=lon,lat=lat)
