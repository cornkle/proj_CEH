import numpy as np
from wavelet import util
import xarray as xr
import matplotlib.pyplot  as plt
import pdb
from scipy import ndimage


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


def run():
    array = xr.open_dataset('/users/global/cornkle/data/Amazon/fract_cover.nc')

    forest = array['nclass'][0,1,:,:].values




    sum = np.sum(array['nclass'][0,:,:,:].values, 0)

    mask = sum == 0

    forest[mask]= 100



    wav = util.waveletSurface(forest, 285)

    forest = np.ma.array(forest, mask=mask)

    wl = wav['power']
    scales = wav['scales']

    for id, s in enumerate(scales):
        wl[id, :,:][wl[id,:,:] <= s**.5] =0
    #
    # maxout = (
    #     wl == ndimage.maximum_filter(wl, (3, 3,3), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))
    #
    #
    # zz, yy, xx = np.where((maxout == 1))
    #
    # print(zz,yy,xx)
    #
    # dummy = np.zeros_like(forest)
    #
    # for z, y, x in zip(zz,yy,xx):
    #
    #     dummy[y,x] = np.round(scales[z])
    #
    # f = plt.figure()
    # plt.imshow(dummy)
    #

    scale_ind = range(len(scales))

    fcnt = 0
    vv = 5
    f = plt.figure()
    for s in scale_ind:

        fcnt+=1
        ax = f.add_subplot(vv,vv,fcnt)

        wll = wl[s,:,:]

        m1 = ax.contour(array['lon'], array['lat'],  wll/np.max(wll)*100, levels=np.arange(0,101,10), cmap='viridis' )
        plt.colorbar(m1)
        ax.set_title(str(wav['scales'][s])+' m')

    ax = f.add_subplot(vv,vv,fcnt+1)
    m = ax.contourf(array['lon'], array['lat'] ,forest, levels=np.arange(0,101,10) )
    plt.colorbar(m)
    #
    # ax = f.add_subplot(vv, vv, fcnt + 2)
    # m = ax.pcolormesh(array['lon'], array['lat'], dummy)
    # plt.colorbar(m)
    for s in scale_ind:
        f = plt.figure()
        plt.contourf(array['lon'], array['lat'],forest)
        plt.contour(array['lon'], array['lat'],wl[s,:,:], cmap='inferno')
        plt.title(str(np.round(wav['scales'][s])))


    f = plt.figure()
    plt.contourf(array['lon'], array['lat'], forest)


    print('Wavelet finished')
    return wav