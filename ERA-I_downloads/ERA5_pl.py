from ecmwfapi import ECMWFDataServer


def download(year):
    server = ECMWFDataServer()
    file = "/localscratch/wllf030/cornkle/ERA5/ERA5_"+str(year)+"_pl.nc"
    server.retrieve({
        "class": "ea",
        "dataset": "era5",
        "date": str(year)+"-06-01/to/"+str(year)+"-09-30",
        "expver": "1",
        "grid": "0.25/0.25",
        "levtype": "pl",
        "levelist": "600/650/700/850/925/950",
        "param": "60.128/130.128/131.128/132.128/133.128/135.128/155.128/157.128",
        "step": "0",
        "stream": "oper",
        "time": "12:00:00/18:00:00/00:00:00/03:00:00",
        "type": "an",
        "area": "22/-12/8/12",
        "format": "netcdf",
        "target": file
    })



for y in range(2008,2011):
    download(y)
