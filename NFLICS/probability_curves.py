import ipdb
import numpy as np
import xarray as xr
from utils import constants as cnst, u_arrays as ua
import pandas as pd
import multiprocessing
import pickle as pkl


def create_tab():
    coref = cnst.network_data + '/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant_daily.nc'
    dat = xr.open_dataarray(coref)

    dic = {'date' : [], 'lon' : [], 'lat': [], 'hour': [], 'month': [], 'year': []}
    for t in (dat.time):
        day = dat.sel(time=t)

        pos = np.where(np.isfinite(day.values))

        for y, x in zip(pos[0], pos[1]):

            point = day.isel(lon=x, lat=y)

            hour = point.values

            date =  str(point['time.year'].values) +'-'+ str(point['time.month'].values).zfill(2) +'-'+ str(point['time.day'].values).zfill(2)+' ' + str(int(hour)).zfill(2) + ':00:00'
            dic['date'].append(date)
            dic['month'].append(int(day['time.month'].values))
            dic['hour'].append(int(hour))
            dic['year'].append(int(day['time.year'].values))
            dic['lon'].append(float(point.lon.values))
            dic['lat'].append(float(point.lat.values))
    pdic = pd.DataFrame.from_dict(dic)

    pdic.to_csv(cnst.network_data + 'data/NFLICS/tables/blob_map_allscales_-50_JJAS_points_dominant_daily.csv',
                na_rep=-999, index_label='id')


def read_SahelStorm(h):

    msgopen = pd.read_csv(
        cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/init_merged2/cores_gt15000km2_table_AMSRE_tracking2_' + str(
            h) + '_init.csv', na_values=[-999, -99])
    msg = pd.DataFrame.from_dict(msgopen)  # &  &
    msg = msg[(msg['topo'] <= 450) & (msg['dtime'] <= 2)]
    date = pd.to_datetime(msg[['year', 'month', 'day', 'hour']])
    msg['date'] = date

    return msg


def read_oldSahel(h):

    pdic = pd.read_csv(cnst.network_data + 'data/NFLICS/tables/blob_map_allscales_-50_JJAS_points_dominant_daily.csv',
                       na_values=[-999])
    pdic = pd.DataFrame.from_dict(pdic)

    tab = pdic[pdic.hour == h]
    return tab



def read_NFLICSstorm(h):  # from core files
    pdic = pd.read_csv(
        '/home/ck/DIR/cornkle/data/NFLICS/tables/coresPower_MSG_-40_700km2_-50points_dominant_vcores_14-4W_'+str(h).zfill(2)+'UTC.csv',
        na_values=[-999])
    pdic = pd.DataFrame.from_dict(pdic)

    return pdic




def run():

    hours = np.arange(15,21)

    pool = multiprocessing.Pool(processes=4)
    hours = [15,16,17,18,19,20,21,22,23]#,0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    res = pool.map(read_lsta, hours)
    pool.close()



def cut_kernel_lsta(xpos, ypos, arr, h):
    dist = 200

    kernel = ua.cut_kernel(arr, xpos, ypos, dist)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    kmean = kernel

    if kernel.shape != (2 * dist + 1, 2 * dist + 1):
        print('WRONG SHAPE')
        return
    ycirc30, xcirc30 = ua.draw_circle(dist + 1, dist + 1, 6)  # 15km radius
    k30 = np.nanmean(kmean[ycirc30, xcirc30])

    ycirc100e, xcirc100e = ua.draw_circle(dist + 51, dist + 1, 25)  # at - 150km, draw 50km radius circle
    data = kmean[ycirc100e, xcirc100e]

    #     if h > 17:
    #         x2 = 95 # 285km upstream
    #     else:
    #         x2 = 80 # 240km
    #     data = kernel[dist - 10:dist + 10, dist+8:dist + x2]

    e100 = np.nanmean(data)

    #     if (np.sum(np.isfinite(data)) / data.size) < 0.01:
    #         print('Too small valid area')
    #         return

    if np.isnan(e100):
        print('Is nan')
        return

    return k30, e100, kernel


def read_lsta(h):

    path = cnst.network_data + '/data/NFLICS/tables/prob_dictionary/'
    dic = read_NFLICSstorm(h)

    tab = dic[dic.hour == h]
    # tab = tab[::400]

    amsrk30 = []
    amsre100 = []

    ramsrk30 = []
    ramsre100 = []

    amskern = None
    ramskern = None

    lonl = []
    latl = []
    tlon = []
    tlat = []

    for date in np.unique(tab.date):

        dt = pd.to_datetime(date)

        # ipdb.set_trace()
        fdate = str(dt.year) + str(dt.month).zfill(2) + str(dt.day).zfill(2)

        lpath = '/home/ck/DIR/cornkle/data/NFLICS/LSTA/netcdf/'

        try:
            # ipdb.set_trace()
            lsta = xr.open_dataset(
                lpath + 'HDF5_LSASAF_ANOM_MSG_LST_MSG-Disk_' + fdate + '1700.nc')  # sma  #'AMSR_L3_LPRMv05_A_'

        except:
            print('Could not find ' + fdate)
            continue

        lsta_da = lsta['lsta']  # .values.squeeze()

        print('Doing ' + 'AMSR_' + str(dt.year) + str(dt.month).zfill(2) + str(
            dt.day).zfill(2) + '.nc')

        # ax.pcolormesh(lsta_da)

        ###############################Blob loop
        cores = 0

        points = np.array(list(zip(lsta.lon.values.flat, (lsta.lat.values).flat)))

        ttab = tab[tab['date'] == date]

        # ipdb.set_trace()

        for lat, lon in zip(ttab.lat, ttab.lon):  # zip([14.7],[ -16.5]):

            # ipdb.set_trace()

            print('Doing ', lat, lon)

            # ax.set_title(str(lat)+' '+ str(lon))

            pout = ua.closest_point((lon, lat), points)
            spos = np.unravel_index(pout, lsta.lon.values.shape)

            ypos = spos[0]
            xpos = spos[1]


            try:
                ak30, ae100, kernel = cut_kernel_lsta(xpos, ypos, lsta_da.values, h)
            except TypeError:
                print('AMSR kernel error')
                continue

            amsrk30.append(ak30)
            amsre100.append(ae100)

            if amskern is not None:
                amskern = np.nansum(np.stack([amskern, kernel]), axis=0)
                amskernc = np.nansum(np.stack([amskernc, np.isfinite(kernel).astype(int)]), axis=0)
            else:
                amskern = kernel
                amskernc = np.isfinite(kernel).astype(int)

            lonl.append(lon)
            latl.append(lat)
            tlon.append(lsta.lon.values[ypos, xpos])
            tlat.append(lsta.lat.values[ypos, xpos])

            cores += 1

            ##### random

            y = ypos
            x = xpos

            # rdist = 50
            # randy50 = [y - rdist, y - rdist, y - rdist, y, y, y + rdist, y + rdist, y + rdist]
            # randx50 = [x - rdist, x, x + rdist, x - rdist, x + rdist, x - rdist, x, x + rdist]
            # randy50_100 = [y - rdist, y - rdist, y, y, y + rdist, y + rdist]
            #
            # rdist = 100
            # randx100 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]
            #
            # rdist = 150
            # randx150 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]
            #
            # rdist = 200
            # randx200 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]
            #
            # randy = np.array(randy50 + randy50_100 + randy50_100 + randy50_100)
            # randx = np.array(randx50 + randx100 + randx150 + randx200)

            # pout = ua.closest_point((lon, 10), points)
            # spos = np.unravel_index(pout, lsta.lon.values.shape)

            xfi = lsta_da.values.shape[1]
            randx = np.random.randint(140, xfi, 30)  # point 140 starts at coast
            if np.min(ypos) == 0:
                ymin = np.min(ypos)
            else:
                ymin = np.min(ypos) - 20
            if np.max(ypos) == lsta_da.values.shape[0] - 1:
                ymax = np.max(ypos)
            else:
                ymax = np.max(ypos) + 20
            randy = np.random.randint(ymin, ymax, 30)
            #posr = (randy, randx)

            rcnt = 0
            for ry, rx in zip(randy, randx):

                # ax.plot(rx, ry, 'ro')
                #print('Doing random', rcnt)

                if ry < 0:
                    continue
                if ry > lsta_da.shape[0] - 1:
                    continue

                if rx < 0:
                    continue
                if rx > lsta_da.shape[1] - 1:
                    continue

                try:
                    arc30, arce100, rkernel = cut_kernel_lsta(rx, ry, lsta_da.values, h)
                except TypeError:
                    print('Random type error')
                    continue
                rcnt +=1
                ramsrk30.append(arc30)
                ramsre100.append(arce100)

                if ramskern is not None:
                   # print('stack')
                    ramskern = np.nansum(np.stack([ramskern, rkernel]), axis=0)
                    ramskernc = np.nansum(np.stack([ramskernc, np.isfinite(rkernel).astype(int)]), axis=0)
                else:
                    ramskern = rkernel
                    ramskernc = np.isfinite(rkernel).astype(int)

            #ipdb.set_trace()

            print('Core lonlat', xpos, ypos)
            print('Random lonlat', rx, ry)

        # ipdb.set_trace()
        del lsta
        del lsta_da

    ams_sum = amskern  # np.nansum(np.stack(amskern, axis=0), axis=0)[np.newaxis,...]
    ams_cnt = amskernc  # np.sum(np.isfinite(np.stack(amskern, axis=0)), axis=0)[np.newaxis,...]

    rams_sum = ramskern  # np.nansum(np.stack(ramskern, axis=0), axis=0)[np.newaxis,...]
    rams_cnt = ramskernc  # np.sum(np.isfinite(np.stack(ramskern, axis=0)), axis=0)[np.newaxis,...]

    dicout = {'amsr': [amsrk30, amsre100],
           'ramsr': [ramsrk30, ramsre100],
           'akern': [ams_sum, ams_cnt],
           'rkern': [rams_sum, rams_cnt],
           'ccord': [lonl, latl],
           'tcord': [tlon, tlat]
           }

    pkl.dump(dicout, open(path + "/Sahel_MCS_Senegal_" + str(h).zfill(2) + ".p", "wb"))



