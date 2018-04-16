from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2006-06-01/to/2010-09-30",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "55.162/84.162/134.128/151.128/165.128/166.128/167.128/168.128/235.128",
    "step": "0",
    "stream": "oper",
    "time": "12:00:00",
    "type": "an",
    "area": "22/-12/8/12",
    "format": "netcdf",
    "target": "/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCsrfc.nc"
})
