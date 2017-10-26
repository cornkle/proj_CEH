import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import multiprocessing



class RainySeason(pd.Series):
    """ Class doc"""

    def __init__(self, data):
        """ Class initialiser """
        assert isinstance(data, pd.Series), "Exp series got {0.__class__}".format(data)
        if data.index.dtype == int:
            idx = data.index
        else:
            idx = np.arange(data.index.size) + 1


        #print(idx.min())
       # print(idx.max())
        if idx.max() < 365:
            print('Problem, not enough days')
            pdb.set_trace()

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

    area = np.array([[-17,-10,10,16], [-10,-2,10,17], [-2,9,10,17], [-3,3,4,12], [-17,-10,13,17]])

    set = area[4]

    # regions: west, central, east, Ghana, Sahel

    start = []
    end = []

    for y1, y2 in zip(yy1,yy2):

        file = y1+'-'+y2+'WA.nc'

        data = xr.open_dataarray(path+file)

        ds = data.sel(X=slice(set[0],set[1]), Y=slice(set[2], set[3]))

        # if y1 == '1983':
        #     f = plt.figure(figsize=(15, 7), dpi=400)
        #     ax = plt.axes(projection=ccrs.PlateCarree())
        #     ds[2,:,:].plot.contourf('X', 'Y', projection=ccrs.PlateCarree(), vmax=8)
        #     ax.coastlines()
        #     # Gridlines
        #     xl = ax.gridlines(draw_labels=True);
        #     xl.xlabels_top = False
        #     xl.ylabels_right = False
        #     # Countries
        #     ax.add_feature(cartopy.feature.BORDERS, linestyle='--');

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

    strii = ' [17-10W, 13-17N]'

    f = plt.figure(figsize=(5,9))
    f.add_subplot(3,1,1)

    plt.plot(np.arange(1983,2017,1),pstart, marker='o')
    plt.plot(np.arange(1983, 2017, 1), pd.rolling_mean(pstart,5, center=True), label='5 yr rolling mean')
    plt.axhline(pstart.mean(), color='k', linestyle='dashed')
    plt.title('west: Start of rainy season (ARC2) |'+strii, fontsize=9)
    plt.legend()

    f.add_subplot(3, 1, 2)
    plt.plot(np.arange(1983, 2017, 1), pend, marker='o')
    plt.plot(np.arange(1983, 2017, 1), pd.rolling_mean(pend, 5, center=True))
    plt.axhline(pend.mean(), color='k', linestyle='dashed')
    plt.title('west: End of rainy season (ARC2) |'+strii, fontsize=9)

    f.add_subplot(3, 1, 3)
    plt.plot(np.arange(1983, 2017, 1), diff, marker='o')
    plt.plot(np.arange(1983, 2017, 1), pd.rolling_mean(diff, 5, center=True))
    plt.axhline(diff.mean(), color='k', linestyle='dashed')
    plt.title('west: Length of rainy season (ARC2) |'+strii, fontsize=9)


def find_start_end(diff):

    dstart= np.argmin(diff)
    dend=np.argmax(diff)



def tamsat(y):

    path = '/users/global/cornkle/data/OBS/TAMSATv3/'

    print('Doing '+str(y))

    area = np.array([[-17,-10,10,17], [-10,-2,10,17], [-2,9,10,17], [-3,3,4,12], [-17,-15,13.5,15.5]]) # [-17,-10,13,17]

    coord = area[4]

    # regions: west, central, east, Ghana, Sahel


    data = xr.open_mfdataset(path + 'rfe'+str(y)+'*.nc')
    data = data.sel(lon=slice(coord[0], coord[1]), lat=slice(coord[3], coord[2]))
    print('Opened data')

    data = data['rfe']

    tstart=[]
    tend = []

    # for yy in np.arange(data.shape[0]):
    #     for xx in np.arange(data.shape[1]):
    mean = data.mean(dim=['lat', 'lon'])
   # mean = data.isel(lat=yy, lon=xx)
    md = mean.to_pandas()
    md = md.reindex(pd.date_range(str(y)+'-01-01', str(y)+'-12-31', freq='D'))

    diff = RainySeason(md)

    dstart= np.argmin(diff)
    dend=np.argmax(diff)

    # f = plt.figure()
    # ax = f.add_subplot(111)
    # plt.plot(diff.index, diff)
    # plt.axvline(dstart, color='k')
    # plt.axvline(dend, color='k')
    # plt.text(dstart-1, -50, str(dstart))
    # plt.text(dend-1, -50, str(dend))
    # plt.minorticks_on()


    # tstart.append(dstart)
    # tend.append(dend)
    # start = np.median(tstart)
    # end = np.median(tend)
    #
    # print('Done ' + str(y))

    return dstart, dend




def multi_loop():

    psave = '/users/global/cornkle/figs/dry_spell/'
    pool = multiprocessing.Pool(processes=1)

    years = np.arange(1983, 2017,1)#2017, 1)
    start=[]
    end = []

    for y in years:
        sstart, send = tamsat(y)
        start.append(sstart)
        end.append(send)

    np.save(psave+'start_season_west_small.npy',np.array(start))
    np.save(psave + 'end_season_west_small.npy', np.array(end))

def plot(str):

    psave = '/users/global/cornkle/figs/dry_spell/'
    start= np.load(psave+'start_season_'+str[0]+'.npy')
    end = np.load(psave + 'end_season_'+str[0]+'.npy')

    pstart = pd.Series(start, index=np.arange(1983, 2017, 1))
    pend = pd.Series(end, index=np.arange(1983, 2017, 1))

    diff = pd.Series(end - start, index=np.arange(1983, 2017, 1))

    strii = str[1]

    f = plt.figure(figsize=(5,9))
    f.add_subplot(3, 1, 1)

    plt.plot(np.arange(1983, 2017, 1), pstart, marker='o')
    plt.plot(np.arange(1983, 2017, 1), pstart.rolling(center=True,window=5).mean(), label='5 yr rolling mean')
    plt.axhline(pstart.mean(), color='k', linestyle='dashed')
    plt.title(str[0]+': Start of rainy season (TAMSATv3)|' + strii, fontsize=9)
    plt.legend()

    f.add_subplot(3, 1, 2)
    plt.plot(np.arange(1983, 2017, 1), pend, marker='o')
    plt.plot(np.arange(1983, 2017, 1), pend.rolling(center=True,window=5).mean())
    plt.axhline(pend.mean(), color='k', linestyle='dashed')
    plt.title(str[0]+': End of rainy season (TAMSATv3)|' + strii, fontsize=9)

    f.add_subplot(3, 1, 3)
    plt.plot(np.arange(1983, 2017, 1), diff, marker='o')
    plt.plot(np.arange(1983, 2017, 1), diff.rolling(center=True,window=5).mean())
    plt.axhline(diff.mean(), color='k', linestyle='dashed')
    plt.title(str[0]+': Length of rainy season (TAMSATv3)|' + strii, fontsize=9)


def loop():

    regions = [('west', '[17-10W, 10-17N]'),
               ('central','[10-2W,10-17N]') ,
               ('east', '[2W-9E,10-17N]'),
               ('west_small', '[17-15W,13.5-15.5N]') ]

    for r in regions:
        plot(r)

#
#
#
#
# if __name__ == "__main__":
#     y_loop()
#
#
# if __name__ == "__main__":
#     y_loop()



