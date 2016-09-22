import xarray as xr
import numpy as np
import pandas as pd


def readStation():
    root = "/users/global/cornkle/data/OBS/stations/"
    nc = root + "MSG_0415_event_cold_cloud_40C_10mm_cellT.nc"
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
