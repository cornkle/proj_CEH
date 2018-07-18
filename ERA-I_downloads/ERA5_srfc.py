from ecmwfapi import ECMWFDataServer



def download(year):
    server = ECMWFDataServer()
    file = "/localscratch/wllf030/cornkle/ERA5/ERA5_"+str(year)+"_srfc.nc"
    server.retrieve({
        "class": "ea",
        "dataset": "era5",
        "date": str(year)+"-06-01/to/"+str(year)+"-09-30",
        "expver": "1",
        "grid": "0.25/0.25",
        "levtype": "sfc",
        "param": "59.128/63.162/81.162/84.162/129.128/134.128/136.128/137.128/151.128/159.128/165.128/166.128/167.128/168.128/231.128/232.128/235.128/246.228/247.228",
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