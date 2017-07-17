import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy



class RainySeason(pd.Series):
    """ Class doc"""

    def __init__(self, data):
        """ Class initialiser """
        assert isinstance(data, pd.Series), "Exp series got {0.__class__}".format(data)
        if data.index.dtype == int:
            idx = data.index
        else:
            idx = np.arange(data.index.size) + 1

        print(idx.min())
        print(idx.max())

        assert idx.min() == 1 and idx.max() in (365, 366) #testing for right format, else AssertionError
        ts = data.fillna(0)
        mrr = ts.mean()
        diff = ts - mrr
        cumdif = diff.cumsum()

        pd.Series.__init__(self, cumdif.values, index=idx)



def y_loop():

    path = '/users/global/cornkle/data/ARC2/'

    yy1 = ['1983', '1990', '2000', '2010' ]
    yy2 = ['1990', '2000', '2010', '2017']

    area = np.array([[-17,-10,10,16],[-17,-14,12,15], [-10,-2,10,17], [-2,9,10,17], [-3,3,4,12]])

    set = area[1]

    # regions: west, central, east, Ghana

    start = []
    end = []

    for y1, y2 in zip(yy1,yy2):

        file = y1+'-'+y2+'WA.nc'

        data = xr.open_dataarray(path+file)

        ds = data.sel(X=slice(set[0],set[1]), Y=slice(set[2], set[3]))

        if y1 == '1983':
            f = plt.figure(figsize=(15, 7), dpi=400)
            ax = plt.axes(projection=ccrs.PlateCarree())
            ds[2,:,:].plot.contourf('X', 'Y', projection=ccrs.PlateCarree(), vmax=8)
            ax.coastlines()
            # Gridlines
            xl = ax.gridlines(draw_labels=True);
            xl.xlabels_top = False
            xl.ylabels_right = False
            # Countries
            ax.add_feature(cartopy.feature.BORDERS, linestyle='--');

        mean = ds.mean(dim=['X', 'Y'])

        years = np.unique(mean['T.year'])



        for y in years:
            print(y)
            md =mean.sel(T=str(y)).to_pandas()

            diff = RainySeason(md)

            start.append(np.argmin(diff))
            end.append(np.argmax(diff))

    pstart = pd.Series(start, index=np.arange(1983, 2017, 1))
    pend = pd.Series(end, index=np.arange(1983, 2017, 1))

    diff = pd.Series(np.array(end)-np.array(start), index=np.arange(1983, 2017, 1))

    strii = ' [17-14W, 12-15N]'

    f = plt.figure()
    f.add_subplot(3,1,1)

    plt.plot(np.arange(1983,2017,1),pstart)
    plt.plot(np.arange(1983, 2017, 1), pd.rolling_mean(pstart,10, center=True), label='10 yr rolling mean')
    plt.title('ARC2: Start of rainy season |'+strii)
    plt.legend()

    f.add_subplot(3, 1, 2)
    plt.plot(np.arange(1983, 2017, 1), pend)
    plt.plot(np.arange(1983, 2017, 1), pd.rolling_mean(pend, 10, center=True))
    plt.title('ARC2: End of rainy season |'+strii)

    f.add_subplot(3, 1, 3)
    plt.plot(np.arange(1983, 2017, 1), diff)
    plt.plot(np.arange(1983, 2017, 1), pd.rolling_mean(diff, 10, center=True))
    plt.title('ARC2: Length of rainy season |'+strii)




if __name__ == "__main__":
    y_loop()



