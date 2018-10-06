from ecmwfapi import ECMWFDataServer


def download(year, month):
    server = ECMWFDataServer()
    file = "/prj/AMMA2050/ERA5/pressure_levels/ERA5_"+str(year)+"_"+str(month).zfill(2)+"_pl.nc"
    server.retrieve({
        "class": "ea",
        "dataset": "era5",
        "date": str(year)+"-"+str(month).zfill(2)+"-01/to/"+str(year)+"-"+str(month).zfill(2)+"-31",
        "expver": "1",
        "grid": "0.25/0.25",
        "levtype": "pl",
        "levelist": "200/250/300/350/400/450/500/550/600/650/700/750/825/850/875/900/925/950/975",
        "param": "60.128/130.128/131.128/132.128/133.128/135.128/155.128/157.128",
        "step": "0",
        "stream": "oper",
        "time": "12:00:00/18:00:00/00:00:00/03:00:00",
        "type": "an",
        "area": "22/-18/4/15",
        "format": "netcdf",
        "target": file
    })


for y in range(2000,2018):
    for m in range(1, 13):
            download(y, m)
