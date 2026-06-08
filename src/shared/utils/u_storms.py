import sortedcontainers as sc



def era5_dic():
        dic = {}
        vars = ['q_col', 'u_col', 'r_col', 'v_col', 'q_col_s', 'u_col_s', 'r_col_s', 'v_col_s', 'u925_s', 'u650_s',
                'q925_s','q700_s', 'tcwv_s', 'CAPE_s', 'tcwv', 'CAPE', 'dates', 'tmin', 'tmean', 't10', 'area',
                'area70', 'lat', 'lon', 'u925', 'u650', 'q925', 'q700']

        for v in vars:
            dic[v] = []
        return dic


class storm(object):

    def __init__(self, tmin, tmean, area, area70, time, lon, lat):

        self.surface = {}
        self.erai = {}
        self.era5 = era5_dic()
        self.tmin = tmin  # exponential factor for calculation of distance between decomposition scales, check resulting scales!
        self.tmean = tmean # the number of scales the data is decomposed into
        self.area = area # pixel resolution (e.g. in km)
        self.area70 = area70
        self.time = time # scales in unit of given pixel resolution
        self.lon = lon # wavelet scales for normalising power spectrum
        self.lat = lat # wavelet scales for normalising power spectrum


    def get_era5(self, era_in):

        time = str(self.time.year)+'-'+str(self.time.month)+'-'+'12'
        stormtime = str(self.time.year)+'-'+str(self.time.month)+'-'+'18'

        try:
            era_day_pl = era_in.sel(ymonth=time)
        except (TypeError, IndexError, KeyError):
            print('Era missing:', time)
            return



        return wav_coeffs_pure, norm_power



class storm_container(storm):

    def __init__(self):

        self.stormlist = sc.SortedDict()

    def add_storm(self):



