from ecmwfapi import ECMWFDataServer
import numpy as np
import pdb
y = np.arange(1979,2018,1)
stri = ''

for yy in y:

    stri = stri+"{0}0101/{0}0201/{0}0301/{0}0401/{0}0501/{0}0601/{0}0701/{0}0801/{0}0901/{0}1001/{0}1101/{0}1201/".format(yy)

stri = stri[0:-1]

server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": stri,
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "81.162/134.128/137.128/164.128/165.128/166.128/167.128/168.128/207.128",
    "stream": "moda",
    "type": "an",
    "area": "34/-23/-42/55",
    "format": "netcdf",
    "target": "/localscratch/wllf030/cornkle/ERA-I/monthly/monthly_1979-2017_srfc_full.nc"
})




