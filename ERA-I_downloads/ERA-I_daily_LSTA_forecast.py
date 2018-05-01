from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2006-06-01/to/2010-09-30",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "59.128/134.128/143.128/146.128/147.128/151.128/159.128/165.128/166.128/167.128/168.128/182.128/201.128/202.128/228.128/235.128/244.128",
    "step": "6",
    "stream": "oper",
    "time": "12:00:00",
    "type": "fc",
    "area": "22/-12/8/12",
    "format": "netcdf",
    "target": "/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCsrfc_forecast.nc"
})
