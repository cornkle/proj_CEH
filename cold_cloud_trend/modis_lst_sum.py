import xarray as xr
import pdb
import numpy as np


y1 = np.arange(2006,2010,1)
y2 = np.arange(2012,2016,1)

path = '/users/global/cornkle/data/MODIS/LST_MOD11C3/'

ylist = [y1, y2]

doy = [1,32,60,91,121,152,182,213,244,274,305,335]
smonth = np.arange(1,13)
for yy in ylist:

    #yy = ylist[1]
    for monthy, sm in zip(doy, smonth):

        fill_array_day = np.zeros((len(np.arange(-90, 90, 0.05)), len(np.arange(-180, 180, 0.05)), 12))+np.nan
        fill_array_night = np.zeros((len(np.arange(-90, 90, 0.05)), len(np.arange(-180, 180, 0.05)), 12))+np.nan

        for y in yy:

            month = monthy

            if ((y == 2008) | (y == 2012)) & (month >32):
                month = month+1
            dmonth = "%03d" % month
            file = 'MOD11C3.A'+str(y)+dmonth+'.hdf'
            print(file)
            tfile = path + file

            tdummy = xr.open_dataset(tfile)
            data_d = np.flip(tdummy['LST_Day_CMG'].values, axis=0)
            data_n = np.flip(tdummy['LST_Night_CMG'].values, axis=0)

            fill_array_day[ :, :, sm-1] = data_d
            fill_array_night[ :, :, sm-1]  = data_n

        fill_array_day = np.nanmean(fill_array_day, axis=2)
        fill_array_night = np.nanmean(fill_array_night, axis=2)
        #pdb.set_trace()
        marray_day = xr.DataArray((fill_array_day)[np.newaxis], coords=[('month', np.array([sm])), ('lat', np.arange(-90,90,0.05)+0.05), ('lon', np.arange(-180,180,0.05)+0.05)], dims=['month', 'lat', 'lon'] )
        marray_day.name = 'lst'
        marray_day = marray_day.sel(lon=slice(-15, 0), lat=slice(3.5, 9))
        marray_day.to_netcdf(path+'clim/'+str(yy[0])+'-'+str(yy[-1])+'_'+dmonth+'_day.nc')


        marray_night = xr.DataArray((fill_array_night)[np.newaxis], coords=[('month', np.array([sm])), ('lat', np.arange(-90,90,0.05)+0.05), ('lon', np.arange(-180,180,0.05)+0.05)], dims=['month' ,'lat', 'lon'] )
        marray_night.name = 'lst'
        marray_night = marray_night.sel(lon=slice(-15,0), lat=slice(3.5, 9))
        marray_night.to_netcdf(path + 'clim/' + str(yy[0]) + '-' + str(yy[-1]) + '_' + dmonth + '_night.nc')

