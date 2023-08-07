from ecmwfapi import ECMWFDataServer
import xarray as xr
import numpy as np

file = "/localscratch/wllf030/cornkle/ERA-I/daily_2004-2014_pl.nc"
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1983-03-01/to/2014-09-30",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "pl",
    "levelist": "600/650/700/850/925/950",
    "param": "60.128/130.128/131.128/132.128/133.128/135.128/155.128/157.128",
    "step": "0",
    "stream": "oper",
    "time": "12:00:00",
    "type": "an",
    "area": "22/-15/3/15",
    "format": "netcdf",
    "target": file
})


