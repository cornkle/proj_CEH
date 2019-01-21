from ecmwfapi import ECMWFDataServer
import xarray as xr
import numpy as np


def loop(y):
  file = "/localscratch/wllf030/cornkle/ERA-I/daily_"+ str(y) + "_pl_12UTC.nc"
  server = ECMWFDataServer()
  server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": str(y)+"-01-01/to/" + str(y) + "-12-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "pl",
    "levelist": "250/350/450/550/600/650/700/750/800/850/900/925/950",
    "param": "60.128/130.128/131.128/132.128/133.128/135.128/155.128/157.128",
    "step": "0",
    "stream": "oper",
    "time": "12:00:00",
    "type": "an",
    "area": "22/-18/3/15",
    "format": "netcdf",
    "target": file
  })

for y in range(1983,2018):
    loop(y)
