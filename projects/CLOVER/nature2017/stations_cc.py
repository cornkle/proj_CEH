import xarray as xr
import numpy as np
import pandas as pd

mfg=0

def readStation(temp, mfg=mfg):
    root = "/users/global/cornkle/data/OBS/stations/"
    nc = root + "MSG_0415_event_cold_cloud_"+str(temp)+"C_10mm_cellT.nc"
    if mfg == 1: nc = root + "event_cold_cloud_"+str(temp)+"C_10mm_cellT.nc"
    txt = root + "thresholds.txt"

    dat = xr.open_dataset(nc)
    rain = np.squeeze(dat['rain'].values / 10)
    minT = np.squeeze(dat['mint'].values)
    area = np.squeeze(dat['area'].values)

    tab = pd.read_csv(txt, header=None, delim_whitespace=True,
                      names=['id', 'lat', 'lon', 'thresh', 'low_tresh', 'name'])

    lat = tab['lat'].as_matrix()
    thresh = tab['thresh'].as_matrix()

    pos = np.where(lat >= 13)

    rain = rain[pos[0],:]
    minT = minT[pos[0],:]
    area = area[pos[0],:]
    thresh = thresh[pos[0]]

    dic = {'rain' : rain, 'minT' : minT, 'area' : area, 'thresh' : thresh}
    return dic

def readStation_both(temp):
    root = "/users/global/cornkle/data/OBS/stations/"
    nc1 = root + "MSG_0415_event_cold_cloud_"+str(temp)+"C_10mm_cellT.nc"
    nc2 = root + "event_cold_cloud_"+str(temp)+"C_10mm_cellT.nc"
    txt = root + "thresholds.txt"

    dat1 = xr.open_dataset(nc1)


    rain1 = np.squeeze(dat1['rain'].values / 10)
    minT1 = np.squeeze(dat1['mint'].values)
    area1 = np.squeeze(dat1['area'].values)

    dat2 = xr.open_dataset(nc2)
    rain2 = np.squeeze(dat2['rain'].values / 10)
    minT2 = np.squeeze(dat2['mint'].values)
    area2 = np.squeeze(dat2['area'].values)
    year = np.squeeze(dat2['eyear'].values)

    tab = pd.read_csv(txt, header=None, delim_whitespace=True,
                      names=['id', 'lat', 'lon', 'thresh', 'low_tresh', 'name'])

    lat = tab['lat'].as_matrix()
    thresh = tab['thresh'].as_matrix()

    pos = np.where(lat >= 13)

    rain1 = rain1[pos[0],:]
    minT1 = minT1[pos[0],:]
    area1 = area1[pos[0],:]
    rain2 = rain2[pos[0], :]
    minT2 = minT2[pos[0], :]
    area2 = area2[pos[0], :]
    year = year[pos[0], :]
    thresh = thresh[pos[0]]

    rain2[(year==2004) ^ (year==2005)] = np.nan
    minT2[(year == 2004) ^ (year == 2005)] = np.nan
    area2[(year == 2004) ^ (year == 2005)] = np.nan

    rain = np.concatenate((rain1, rain2), axis=1)
    minT = np.concatenate((minT1, minT2), axis=1)
    area = np.concatenate((area1, area2), axis=1)

    dic = {'rain' : rain, 'minT' : minT, 'area' : area, 'thresh' : thresh}
    return dic
