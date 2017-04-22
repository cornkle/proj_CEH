import numpy as np
from wavelet import util
import xarray as xr
import matplotlib.pyplot  as plt
import pdb
from scipy import ndimage
from utils import u_arrays as ua
import multiprocessing


def read_grads(file):
    # if not file:
    #     print('I found no msg files in my fpath')
    #     return False
    #
    # rrShape = (2422, 2316)  # msg shape
    # rrMDI = np.uint8(255)
    #
    # rr = np.fromfile(file, dtype=rrMDI.dtype)
    # rr.shape = rrShape
    # rr = rr.astype(np.int32) - 173
    #
    # if llbox:
    #     i, j = np.where(
    #         (self.lon > llbox[0]) & (self.lon < llbox[2]) & (self.lat > llbox[1]) & (self.lat < llbox[3]))
    #     blat = self.lat[i.min():i.max() + 1, j.min():j.max() + 1]
    #     blon = self.lon[i.min():i.max() + 1, j.min():j.max() + 1]
    #     rr = rr[i.min():i.max() + 1, j.min():j.max() + 1]
    # else:
    #     blat = self.lat
    #     blon = self.lon
    #     rr = rr
    #
    # str = file.split(os.sep)[-1]
    # curr_date = [
    #     pd.datetime(np.int(str[0:4]), np.int(str[4:6]), np.int(str[6:8]), np.int(str[8:10]), np.int(str[10:12]))]
    # date = curr_date  # or np.atleast_1d(dt.datetime())
    #
    # da = xr.DataArray(rr[None, ...], coords={'time': (('time'), date),
    #                                          'lat': (('y', 'x'), blat),
    #                                          'lon': (('y', 'x'), blon)},
    #                   dims=['time', 'y', 'x']).isel(time=0)
    #
    # ds = xr.Dataset({'t': da})
    #
    # if netcdf_path:
    #     savefile = netcdf_path
    #
    #     try:
    #         os.remove(savefile)
    #     except OSError:
    #         pass
    #     try:
    #         ds.to_netcdf(path=savefile, mode='w')
    #     except OSError:
    #         print('File cannot be saved. Maybe directory wrong. Check your permissions')
    #         raise
    #
    #     print('Saved ' + savefile)

    pass


def run(array):

    forest = array['nclass'][0, 1, :, :].values

    sum = np.sum(array['nclass'][0, :, :, :].values, 0)

    mask = sum == 0

    forest[mask] = 100



    wav = util.waveletSurface(forest, 285)

    forest = np.ma.array(forest, mask=mask)

    wl = wav['power']
    scales = wav['scales']

    for id, s in enumerate(scales):
        wl[id, :,:][wl[id,:,:] <= s**.5] =0

    return wl, scales, forest


def save():

    files = ua.locate(".nc", '/users/global/cornkle/data/Amazon')
   # files = files[6:7]

    years = len(files)
    yy = np.arange(1984, 2009, 4)

    scale_id = 7

    ylist = []
    vlist = []

    for i, f in enumerate(files):
        y = yy[i]

        print('Doing ' +f)

        array = xr.open_dataset(f)

        wl, scales, forest = run(array)

        print(scales)

        wwl = wl[scale_id]

        yarr = xr.DataArray(wwl[np.newaxis,...], coords=[array.time, array.lat, array.lon], dims=['time', 'lat', 'lon'])
        veg = xr.DataArray(forest[np.newaxis,...], coords=[array.time, array.lat, array.lon], dims=['time', 'lat', 'lon'])

        ylist.append(yarr)
        vlist.append(veg)

    yarr = xr.concat(ylist, dim='time')
    veg = xr.concat(vlist, dim='time')


    sc = int(scales[scale_id]/1000.)

    xarr = xr.Dataset()
    xarr['vegfra']=veg
    xarr['wav']=yarr

    xarr.to_netcdf('/users/global/cornkle/amazon/nc/rhod_'+str(sc)+'kmt.nc')
    print('Saved '+'/users/global/cornkle/amazon/nc/rhod_'+str(sc)+'kmt.nc')


def plot():
    sc = 3
    f = '/users/global/cornkle/amazon/nc/rhod_'+str(sc)+'kmneg.nc'
    array = xr.open_dataset(f)

    g = array['vegfra'].plot.contourf('lon', 'lat', col='time', col_wrap=4, robust=True)

    for i, ax in enumerate(g.axes.flat):
        try:
            ax.contour(array.lon, array.lat, array['wav'][i,:,:], cmap='Blues')
        except (IndexError, ValueError):
            continue


    plt.savefig('/users/global/cornkle/amazon/rhod_'+str(sc)+'kmneg.png', dpi=400)


    # for i, f in enumerate(files):
    #     y = yy[i]
    #
    #     print('Doing ' +f)
    #
    #     array = xr.open_dataset(f)
    #
    #     wl, scales, forest = run(array)
    #
    #     print(scales)
    #
    #     scale_ind = range(len(scales))
    #
    #     ax = fig.add_subplot(4,2, i+1)
    #
    #     wll = wl[scale_id,:,:]
    #
    #     m = ax.contourf(array['lon'], array['lat'], forest, levels=np.arange(0, 101, 10), cmap='viridis')
    #     m1 = ax.contour(array['lon'], array['lat'],  wll/np.max(wll)*100, levels=np.arange(0,101,10), cmap='Blues' )
    #     plt.colorbar(m)
    #     ax.set_title(str(scales[scale_id])+' m '+str(y) )
    #
    #     xr.Dataset.close(array)
    #
    #
    # plt.savefig('/users/global/cornkle/rhod_30km.png', dpi=400)
    # plt.show()





if __name__ == "__main__":
    plot()