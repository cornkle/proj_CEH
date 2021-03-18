import ipdb
import numpy as np
import xarray as xr
from utils import constants as cnst, u_arrays as ua
import pandas as pd
import multiprocessing
import pickle as pkl
import glob
from scipy.stats import norm
from utils import u_statistics as u_stat
import matplotlib.pyplot as plt




def read_bigDomain(h):  # from core files
    tables = glob.glob(
        '/media/ck/Elements/Africa/WestAfrica/NFLICS/tables/blob_bigDomain_allscales_-40_9-130km_points_dominant_daily_NFLICSdomain_*.csv')
    tab = pd.read_table(tables[0], parse_dates=True, delimiter=',')
    for t in tables[1::]:
        dummy = pd.read_table(t, parse_dates=True, delimiter=',')
        tab = pd.concat([tab, dummy], axis=0)
    tab = tab[tab.hour == h]
    return tab




def run():

    hours = np.arange(15,21)

    pool = multiprocessing.Pool(processes=3)
    hours = [10,11]#[0,1,2,3,4,5,6,7,8,9,10,11,12]#[13,14,15,16,17,18,19,20,21,22,23]
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

    path = cnst.elements_drive + '/Africa/WestAfrica/NFLICS/tables/prob_dictionary'
    tab = read_bigDomain(h)

    amsrk30 = []
    amsre100 = []

    ramsrk30 = []
    ramsre100 = []

    amskern = None
    ramskern = None

    amskernc = None
    ramskernc = None

    lonl = []
    latl = []
    tlon = []
    tlat = []

    #ipdb.set_trace()

    for date in np.unique(tab.date):

        single = tab[tab.date == date]

        dt = pd.to_datetime(date)
        if h <= 10:
            shift=-1
        else:
            shift=0
        daystring = str(abs(shift))
        dayd = pd.Timedelta(daystring + 'days')
        if shift < 0:
            dt = dt - dayd
        if shift >=0:
            dt = dt + dayd

        fdate = str(dt.year) + str(dt.month).zfill(2) + str(dt.day).zfill(2)

        lpath = cnst.elements_drive + '/Africa/WestAfrica/NFLICS/LSTA_2004-2015/netcdf_onCores/'

        try:
            # ipdb.set_trace()
            lsta = xr.open_dataset(
                lpath + 'HDF5_LSASAF_ANOM_MSG_LST_MSG-Disk_' + fdate + '1700.nc')  # sma  #'AMSR_L3_LPRMv05_A_'

        except:
            print('Could not find ' + fdate)
            continue

        print('Doing', date)
        for index, row in single.iterrows():

            lsta_da = (lsta['lsta'].squeeze()) / 100
            lsta_da.values[lsta_da.values == 0] = np.nan

            print('Doing ' + 'AMSR_' + str(dt.year) + str(dt.month).zfill(2) + str(
                dt.day).zfill(2) + '.nc')

            point = lsta_da.sel(lat=row.lat, lon=row.lon, method='nearest')

            plat = point['lat'].values
            plon = point['lon'].values

            xpos = np.where(lsta_da['lon'].values == plon)
            xpos = int(xpos[0])
            ypos = np.where(lsta_da['lat'].values == plat)
            ypos = int(ypos[0])

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

            lonl.append(row.lon)
            latl.append(row.lat)
            tlon.append(point.lon.values)
            tlat.append(point.lat.values)


            ##### random

            y = ypos
            x = xpos

            rdist = 50
            randy50 = [y - rdist, y - rdist, y - rdist, y, y, y + rdist, y + rdist, y + rdist]
            randx50 = [x - rdist, x, x + rdist, x - rdist, x + rdist, x - rdist, x, x + rdist]
            randy50_100 = [y - rdist, y - rdist, y, y, y + rdist, y + rdist]

            rdist = 100
            randx100 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]

            rdist = 150
            randx150 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]

            rdist = 200
            randx200 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]

            randy = np.array(randy50 + randy50_100 + randy50_100 + randy50_100)
            randx = np.array(randx50 + randx100 + randx150 + randx200)

            # pout = ua.closest_point((lon, 10), points)
            # spos = np.unravel_index(pout, lsta.lon.values.shape)
######################################
            # xfi = lsta_da.values.shape[1]
            # randx = np.random.randint(140, xfi, 30)  # point 140 starts at coast
            # if np.min(ypos) == 0:
            #     ymin = np.min(ypos)
            # else:
            #     ymin = np.min(ypos) - 20
            # if np.max(ypos) == lsta_da.values.shape[0] - 1:
            #     ymax = np.max(ypos)
            # else:
            #     ymax = np.max(ypos) + 20
            # randy = np.random.randint(ymin, ymax, 30)
            ##posr = (randy, randx)
#############################################
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
            #print('Random lonlat', rx, ry)

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

    pkl.dump(dicout, open(path + "/BigDomain_LSTAfromBlob_17W-1E-9-20N_randomFixed_" + str(h).zfill(2) + ".p", "wb"))



def test_file(savepath):


    h=17
    path = cnst.elements_drive + '/Africa/WestAfrica/NFLICS/tables/prob_dictionary/'
    dic = pkl.load(open(path + savepath, "rb"))

    ###### Kernel
    plt.figure()
    plt.pcolormesh(dic['akern'][0] / dic['akern'][1], cmap='RdBu_r', vmin=-2, vmax=2)
    plt.colorbar()


    ###### Histogram comparison
    cinput = np.array(dic['amsr'][0])
    rinput = np.array(dic['ramsr'][0])
    cinput = cinput[np.isfinite(cinput)]
    rinput = rinput[np.isfinite(rinput)]

    nbpoint, pointcount, bins = u_stat.histo_frequency(cinput, bins=np.arange(-10, 10, 0.5))
    nball, allcount, bins = u_stat.histo_frequency(rinput, bins=np.arange(-10, 10, 0.5))

    # nbpoint, bins, pointcount = plt.hist(cinput, bins=np.arange(-12, 12, 0.1), density=1)
    # nball, bins, allcount = plt.hist(rinput, bins=np.arange(-12, 12, 0.1), density=1)
    print(bins)
    bin_centre = bins[0:-1] + ((bins[1::] - bins[0:-1]) / 2)
    bin_edge = bins[0:-1]
    width = bins[1::] - bins[0:-1]

    f = plt.figure(figsize=(9, 5), dpi=200)
    ax = f.add_subplot(111)

    ax.bar(bin_edge, nbpoint, label='Core', edgecolor='k', alpha=0.5, align='edge', width=width)
    ax.bar(bin_edge, nball, label='Random', edgecolor='k', alpha=0.5, align='edge', width=width)
    plt.ylabel('Frequency')
    stri = (np.sum(cinput >= np.percentile(rinput, 75)) / cinput.size * 100).round(2)
    # ipdb.set_trace()
    plt.title(str(stri) + '% of Cells occur in warmest 75% | ' + str(np.round(stri - 25, 2)) + '% more')
    plt.legend()

    stri = (np.sum(cinput >= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
    print(str(stri) + '% of Cells occur in warmest half')

    plt.axvline(0, color='k')

    ########### Curve fit on histogram

    # # from scipy.stats import norm
    # # best fit of data
    (mupoint, sigmapoint) = norm.fit(cinput)
    (muall, sigmaall) = norm.fit(rinput)
    #
    # # add a 'best fit' line
    ypoint = norm.pdf(bins, mupoint, sigmapoint)
    yall = norm.pdf(bins, muall, sigmaall)


    ################ Other curve fit
    # fitting function:
    # Gaussian function

    def fit_curve_scipy(x,y):
        from scipy.optimize import curve_fit
        def gauss_function(x, a, x0, sigma):
            return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

        # estimate mean and standard deviation
        mean = sum(x * y)
        sigma = sum(y * (x - mean) ** 2)
        # do the fit!
        popt, pcov = curve_fit(gauss_function, x, y, p0=[1, mean, sigma])

        new_y = gauss_function(x, *popt)
        return new_y

    def fit_curve_lm(x,y):
        from lmfit import Model
        from numpy import exp, linspace, random

        def gaussian(x, amp, cen, wid):
            """1-d gaussian: gaussian(x, amp, cen, wid)"""
            return (amp / (np.sqrt(2 * np.pi) * wid)) * exp(-(x - cen) ** 2 / (2 * wid ** 2))

        gmodel = Model(gaussian)
        result = gmodel.fit(y, x=x, amp=5, cen=5, wid=1)

        return result.best_fit

    def fit_curve_poly(x,y):
        z = np.polyfit(x, y, 5)
        f = np.poly1d(z)

        return f





    # ###### Calculate probability from fitted curve
    #
    def calc_prob_pdf(x, mu_core, mu_random, sigma_core, sigma_random):
        corep = norm.pdf(x, mu_core, sigma_core)
        randomp = norm.pdf(x, mu_random, sigma_random)

        corec = norm.cdf(x, mu_core, sigma_core)
        randomc = norm.cdf(x, mu_random, sigma_random)

        out_pdf = (corep - randomp) / randomp

        prob_factor = 1 + out_pdf

        return out_pdf, prob_factor

    # pdf based
    prob_pdf = []
    factor_pdf = []
    for b in bins:
        tb, fac = calc_prob_pdf(b, mupoint, muall, sigmapoint, sigmaall)
        prob_pdf.append(tb)
        factor_pdf.append(fac)

    #histogram based
    prob_histo = (nbpoint-nball)/nball
    factor_histo = 1 + prob_histo

    #curve fit scipy
    nball_fit = fit_curve_scipy(bin_centre, nball)
    nbpoint_fit = fit_curve_scipy(bin_centre, nbpoint)
    #
    prob_scipy = (nbpoint_fit-nball_fit)/nball_fit
    factor_scipy = 1 + prob_scipy


    #curve fit lm
    nball_fit = fit_curve_lm(bin_centre, nball)
    nbpoint_fit = fit_curve_lm(bin_centre, nbpoint)
    #
    prob_lm = (nbpoint_fit-nball_fit)/nball_fit
    factor_lm = 1 + prob_scipy


    #polyfit
    # calculate polynomial
    probf = fit_curve_poly(bin_centre,prob_histo)
    factorf = fit_curve_poly(bin_centre, factor_histo)

    x_new = np.linspace(-10,10.1,1000)
    prob_poly = probf(x_new)
    factor_poly = factorf(x_new)


    ########################

    f = plt.figure(figsize=(9, 5), dpi=200)
    ax = f.add_subplot(111)

    ax.bar(bin_edge, nbpoint, label='Core', edgecolor='k', alpha=0.5, align='edge', width=width)
    ax.bar(bin_edge, nball, label='Random', edgecolor='k', alpha=0.5, align='edge', width=width)

    # plt.plot(bins, ypoint, 'b--', linewidth=2)
    # plt.plot(bins, yall, 'r--', linewidth=2)

    # plt.plot(x_new, y_point, 'b--', linewidth=2)
    # plt.plot(x_new, y_all, 'r--', linewidth=2)

    plt.ylabel('Probability')
    stri = (np.sum(cinput >= np.percentile(rinput, 75)) / cinput.size * 100).round(2)
    enh = (stri - 25) / 25 * 100
    plt.title(str(stri) + '% of Cells occur in warmest 25% | ' + str(np.round(enh)) + '% enhancement')
    plt.legend()

    stri = (np.sum(cinput >= np.percentile(rinput, 50)) / cinput.size * 100).round(2)
    enh = (stri - 50) / 50 * 100
    print(str(stri) + '% of Cells occur in warmest half', 'Enhancement:' + str(np.round(enh, 2)) + '%')

    plt.axvline(0, color='k')







    ###### Plot probability curves
   # ipdb.set_trace()
    f = plt.figure(figsize=(12, 6), dpi=300)
    lw = 1
    ax = f.add_subplot(121)
    ax.plot(bins, np.array(prob_pdf), label='PDF-based', lw=lw)
    ax.plot(bin_centre, np.array(prob_histo), marker='o', label='Raw-Histo-based', lw=lw)
    ax.plot(bin_centre, np.array(prob_lm), label='LM-based', lw=2)
    ax.plot(bin_centre, np.array(prob_scipy),  label='Scipy-based', lw=lw)
    ax.plot(x_new, np.array(prob_poly), label='Polyfit-raw', lw=lw)
    plt.ylabel('Probability ')
    plt.axhline(0, color='k', linestyle='dashed')
    plt.axvline(0, color='k', linestyle='dashed')
    plt.xlim(-10, 10)
    plt.xlabel('LSTA (K)')
    plt.title('Probability difference: (CoreP-RandomP)/RandomP')
    plt.legend()

    ax = f.add_subplot(122)
    ax.plot(bins, np.array(factor_pdf),  label='PDF-based', lw=lw)
    ax.plot(bin_centre, np.array(factor_histo), marker='o', label='Raw-Histo-based', lw=lw)
    ax.plot(bin_centre, np.array(factor_lm),  label='LM-based', lw=2)
    ax.plot(bin_centre, np.array(factor_scipy),  label='Scipy-based', lw=lw)
    ax.plot(x_new, np.array(factor_poly), label='Polyfit-raw', lw=lw)
    plt.ylabel('-')
    plt.axhline(1, color='k', linestyle='dashed')
    plt.axvline(0, color='k', linestyle='dashed')
    plt.xlim(-10, 10)
    plt.xlabel('LSTA (K)')
    plt.title('Probability factor (1-P)')
    plt.savefig('/home/ck/DIR/cornkle/figs/NFLICS/LSTA_probability_tests.jpg')

    #plt.close('all')


def calculate_hourly_ptables():

    def calc_prob_pdf(x, mu_core, mu_random, sigma_core, sigma_random):
        corep = norm.pdf(x, mu_core, sigma_core)
        randomp = norm.pdf(x, mu_random, sigma_random)

        out_pdf = (corep - randomp) / randomp

        prob_factor = 1 + out_pdf

        return out_pdf, prob_factor

    def fit_curve_poly(x,y, h):
        if (h <=10) | (h>=23):
            order = 2
        else:
            order = 6
        z = np.polyfit(x, y, order)
        f = np.poly1d(z)

        return f


    hours = np.arange(0,24)

    for h in hours:

        f = plt.figure(figsize=(15, 5), dpi=300)
        lw = 2
        ax = f.add_subplot(131)
        ax1 = f.add_subplot(132)
        ax2 = f.add_subplot(133)

        path = cnst.elements_drive + '/Africa/WestAfrica/NFLICS/tables/prob_dictionary/BigDomain_LSTAfromBlob_17W-1E-9-20N_randomFixed_'+str(h).zfill(2)+'.p'
        dic = pkl.load(open(path , "rb"))

        ###### Histogram comparison
        cinput = np.array(dic['amsr'][0])
        rinput = np.array(dic['ramsr'][0])
        cinput = cinput[np.isfinite(cinput)]
        rinput = rinput[np.isfinite(rinput)]

        print(np.min(cinput), np.max(cinput))
        print(np.min(rinput), np.max(rinput))

        # nbpoint, pointcount, bins = u_stat.histo_frequency(cinput, bins=np.arange(-10, 10.1, 1))
        # nball, allcount, bins = u_stat.histo_frequency(rinput, bins=np.arange(-10, 10.1, 1))

        nbpoint, pointcount, bins = u_stat.histo_frequency(cinput, bins=np.arange(np.percentile(rinput,0.1), np.percentile(rinput,99.9)+0.1,1))
        nball, allcount, bins = u_stat.histo_frequency(rinput, bins=np.arange(np.percentile(rinput,0.1), np.percentile(rinput,99.9)+0.1,1))

        # nbpoint, bins, pointcount = plt.hist(cinput, bins=np.arange(-12, 12, 0.1), density=1)
        # nball, bins, allcount = plt.hist(rinput, bins=np.arange(-12, 12, 0.1), density=1)
        print(bins)
        bin_centre = bins[0:-1] + ((bins[1::] - bins[0:-1]) / 2)
        bin_edge = bins[0:-1]
        width = bins[1::] - bins[0:-1]

        x_new = np.linspace(np.percentile(bin_centre, 10), np.percentile(bin_centre, 90) + 0.1, len(bin_centre))

        ########### Curve fit on histogram

        # # from scipy.stats import norm
        # # best fit of data

        # pdf based
        (mupoint, sigmapoint) = norm.fit(cinput)
        (muall, sigmaall) = norm.fit(rinput)
        prob_pdf, factor_pdf = calc_prob_pdf(x_new, mupoint, muall, sigmapoint, sigmaall)

        # histogram based polyfit
        prob_histo = (nbpoint - nball) / nball
        factor_histo = 1 + prob_histo
        # calculate polynomial function
        probf = fit_curve_poly(bin_centre, prob_histo, h)
        factorf = fit_curve_poly(bin_centre, factor_histo, h)
        #insert new x
        prob_poly = probf(x_new)
        factor_poly = factorf(x_new)

        ax.plot(bin_centre, np.array(factor_histo), marker='o', label='Raw-Histo-based', lw=1, ms=0.5, color='k')
        ax.plot(x_new, np.array(factor_pdf), lw=lw, label='PDF-based')


        ax1.plot(bin_centre, np.array(factor_histo), marker='o', label='Raw-Histo-based', lw=1, ms=0.5, color='k')
        ax1.plot(x_new, np.array(factor_poly), label='Polyfit-raw', lw=lw)

        ax2.bar(bin_edge, nbpoint, label='Core', edgecolor='k', alpha=0.5, align='edge', width=width)
        ax2.bar(bin_edge, nball, label='Random', edgecolor='k', alpha=0.5, align='edge', width=width)

        plt.ylabel('Probability')
        stri = (np.sum(cinput >= np.percentile(rinput, 75)) / cinput.size * 100).round(2)
        enh = (stri - 25) / 25 * 100
        plt.title(str(stri) + '% of Cells in warmest 25% | ' + str(np.round(enh)) + '% enhancement')
        ax2.legend()
        plt.axvline(0, color='k')

        ax.set_ylabel('-')
        ax.axhline(1, color='k', linestyle='dashed')
        ax.axvline(0, color='k', linestyle='dashed')
        ax.set_xlim(-15, 15)
        ax.set_xlabel('LSTA (K)')
        ax.set_title('Probability factor (1-P): PDF-based')


        ax1.set_ylabel('-')
        ax1.axhline(1, color='k', linestyle='dashed')
        ax1.axvline(0, color='k', linestyle='dashed')
        ax1.set_xlim(-15, 15)
        ax1.set_xlabel('LSTA (K)')
        ax1.set_title('Probability factor (1-P): Polyfit-based')
        plt.tight_layout()
        plt.savefig('/home/ck/DIR/cornkle/figs/NFLICS/probability_test_hourly/LSTA_probability_tests_'+str(h).zfill(2)+'_30km.jpg')

        plt.close('all')

        Pdic_local ={}
        Pdic_local['bins'] = bin_centre
        Pdic_local['factor'] = factor_histo

        df = pd.DataFrame.from_dict(Pdic_local)
        dpath = cnst.elements_drive + '/Africa/WestAfrica/NFLICS/tables/prob_dictionary/BigDomain_LSTAprobability_17-0W_9-20N_'+str(h).zfill(2)+'.csv'
        df.to_csv(dpath)

