from ecmwfapi import ECMWFDataServer
from utils import constants as cnst



def download(year,month):
    server = ECMWFDataServer()
    file = cnst.local_data + "ERA5/ERA5_"+str(year)+"_"+str(month).zfill(2)+"_srfc.nc"
    server.retrieve({
        "class": "ea",
        "dataset": "era5",
        "date": str(year)+"-"+str(month).zfill(2)+"-01/to/"+str(year)+"-"+str(month).zfill(2)+"-31",
        "expver": "1",
        "grid": "0.25/0.25",
        "levtype": "sfc",
        "param":"59.128/60.162/61.162/62.162/63.162/71.162/72.162/78.128/79.128/79.162/80.162/81.162/82.162/83.162/84.162/86.162/89.228/129.128/134.128/137.128/151.128/159.128/164.128/165.128/166.128/167.128/168.128/172.128/186.128/187.128/188.128/231.128/232.128/235.128/246.228/247.228" ,
        "step": "0",
        "stream": "oper",
        "time": "12:00:00/18:00:00/00:00:00/03:00:00",
        "type": "an",
        "area": "22/-18/4/15",
        "format": "netcdf",
        "target": file
    })


for y in range(2000,2006):
    for m in range(1,13):
       download(y,m)
