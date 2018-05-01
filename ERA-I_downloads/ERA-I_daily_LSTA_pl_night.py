from ecmwfapi import ECMWFDataServer
import xarray as xr
import numpy as np

file = "/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCpl_night.nc"
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2006-06-01/to/2010-09-30",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "pl",
    "levelist": "600/700/850/925",
    "param": "60.128/129.128/130.128/131.128/132.128/133.128/135.128/155.128/157.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00",
    "type": "an",
    "area": "22/-12/8/12",
    "format": "netcdf",
    "target": file
})


