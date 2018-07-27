import datetime as dt
from eod import trmm_clover
import pickle as pkl
import pdb
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from utils import u_grid

HOD = range(24)  # hours of day
YRANGE = range(2004, 2016)


def save_values():
    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_tropWA' #meteosat_WA30'

    t = trmm_clover.ReadWA(trmm_folder, yrange=YRANGE, area=[-14, 12, 4, 9])   # [-15, 15, 4, 21], [-10, 10, 10, 20]
    value_list = []
    for _y, _m, _d, _h, _mi in zip(t.dates.dt.year,  t.dates.dt.month, t.dates.dt.day, t.dates.dt.hour, t.dates.dt.minute):

        if (_h <10) | (_h>19):
            continue

        date = dt.datetime(_y, _m, _d, _h, _mi)

        tdic = t.get_ddata(date, cut=[4.3, 8])

        print('TRMM:', date)

        value_list.extend((tdic['p'].values)[tdic['p'].values>=0.1].flatten())

    pkl.dump(np.array(value_list), open('/users/global/cornkle/data/CLOVER/saves/trmm_values_2004-2016MAM_10-19UTC_-14W12E_4-9N_0.1mm.p', 'wb'))


def save_values_regrid():
    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_tropWA' #meteosat_WA30'

    t = trmm_clover.ReadWA(trmm_folder, yrange=YRANGE, area=[-14, 12, 4, 9])   # [-15, 15, 4, 21], [-10, 10, 10, 20]
    value_list = []
    for _y, _m, _d, _h, _mi in zip(t.dates.dt.year,  t.dates.dt.month, t.dates.dt.day, t.dates.dt.hour, t.dates.dt.minute):

        if (_h <10) | (_h>19):
            continue

        date = dt.datetime(_y, _m, _d, _h, _mi)

        tdic = t.get_ddata(date, cut=[4.3, 8])
        lon = np.arange(-14,13,1)
        lat = np.arange(4,9,1)
        # make salem grid

        grid = u_grid.make(lon, lat, 5000)
        lon, lat = grid.ll_coordinates

        # interpolate TRM and MSG to salem grid
        inter, tpoints = u_grid.griddata_input(tdic['lon'].values, tdic['lat'].values, grid)

        # Interpolate TRMM using delaunay triangularization
        try:
            dummyt = griddata(tpoints, tdic['p'].values.flatten(), inter, method='linear')
        except ValueError:
            continue
        outt = dummyt.reshape((grid.ny, grid.nx))

        ##remove edges of interpolated TRMM
        for nb in range(5):
            boole = np.isnan(outt)
            outt[boole] = -1000
            grad = np.gradient(outt)
            outt[boole] = np.nan
            outt[abs(grad[1]) > 300] = np.nan
            outt[abs(grad[0]) > 300] = np.nan

        print('TRMM:', date)

        value_list.extend((outt)[outt>=0.1].flatten())

    pkl.dump(np.array(value_list), open('/users/global/cornkle/data/CLOVER/saves/trmm_values_2004-2016MAM_10-19UTC_-14W12E_4-9N_0.1mm_regrid.p', 'wb'))




def get_extreme():

    dic = pkl.load(open('/users/global/cornkle/data/CLOVER/saves/trmm_values_2004-2016MAM_10-19UTC_-14W12E_4-9N_0.1mm_regrid.p',
                        'rb'))  # MSG_TRMM_temp_pcp_300px2004-2013_new.p', 'rb'))

    print(np.percentile(dic,99))
    plt.figure()
    plt.hist(dic)

    #not regridded: ~ 39mm at > 0.1mm
    #regridded: ~ 31mm at >0.1
    #not regridded: ~ 43mm at > 1mm

    pass