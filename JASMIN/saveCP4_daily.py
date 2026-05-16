# -*- coding: utf-8 -*-


import glob
import xarray as xr


def saveDaily():

    files = glob.glob('/prj/AMMA2050/CP4/historical/25km/precip/*.nc')

    for svar in files:

        out = svar.replace('A1hr_mean', 'Aday_mean')
        out = out.replace('precip', 'precip_day')

        ds = xr.open_dataset(svar)
        ds = ds.groupby('time.day').sum('time')

        ds.to_netcdf(out)