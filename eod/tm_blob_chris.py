import numpy as np
import pickle as pkl
import datetime as dt
from eod import trmm, msg, tm_utils
import ipdb
from collections import defaultdict
import matplotlib.pyplot as plt
from utils import u_arrays as ua


import cartopy.crs as ccrs




YRANGE=range(2004,2014)
#AREA=[-16,4,16,20]
AREA=None

def tm_overlap_blobs():

    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_WA30'

    tObj = trmm.ReadWA(trmm_folder, area=AREA, yrange=YRANGE)
    mObj = msg.ReadMsg(msg_folder)

    files = tObj.fpaths
    dates = tObj.dates

    mdic = defaultdict(list)
    mdic_f = defaultdict(list)

    mlon = mObj.lon
    mlat = mObj.lat

    mll = tm_utils.ll_toMSG(mlon, mlat)
    mxy = ua.unique_of_pair(mll['x'], mll['y'])

    cnt = 0
    datess = []

    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered
    for _y, _m, _d, _h, _mi in zip(dates.y, dates.m, dates.d, dates.h, dates.mi):

        # set zero shift time for msg
        date = dt.datetime(_y, _m, _d, _h, _mi)

        dt0 = tm_utils.minute_delta(_mi, 30)
        print('TRMM', date, 'dt', dt0, 'MSG', date + dt.timedelta(minutes=int(dt0)) )
        #time difference max
        # if abs(dt0) > 4:
        #     continue

        ndate = date + dt.timedelta(minutes=int(dt0))
        print('TRMM', date, 'MSG', ndate)

        mObj.set_date(ndate.year, ndate.month, ndate.day , ndate.hour, ndate.minute)

        if not (mObj.tpath or mObj.bpath):
            print('No table or blob file, continue')
            continue

        dff = mObj.get_table()
        dstring = str(ndate.year) + '-' + str(ndate.month).zfill(2) + '-' + str(ndate.day).zfill(2) + ' ' + str(ndate.hour).zfill(2) + ':' + str(
            ndate.minute).zfill(2) + ':' + str(00).zfill(2)
        if not dstring in dff['Date'].as_matrix():
            continue

        sel = dff.loc[dff['Date'] == dstring]
        big = sel.loc[sel['Area'] >= 25000]  # only mcs over 25.000km2
        print('big area', big['Area'].values)
        if big.empty:
            continue

        td = tObj.get_ddata(_y, _m, _d, _h, _mi, cut=[0,22])
        try:
            if not td:
                print('TRMM problem')
                continue
        except:
            pass

        md = mObj.get_data(llbox=AREA)

        md_blob = mObj.get_blob(llbox=AREA)

        blobs = md_blob.values

        blat = big['Lat'].values.tolist()
        blon = big['Lon'].values.tolist()
        barea = big['Area'].values.tolist()
        btemp = big['Temp'].values.tolist()

        for lon, lat, bt, ba in zip(blon, blat, btemp, barea):

            mcs = tm_utils.ll_toMSG(lon, lat)
            point = np.where((mll['x'] == mcs['x']) & (mll['y'] == mcs['y']))

            # if not all(point):
            #     if mcs['x'] > mll['x'].max() or mcs['x'] < mll['x'].min() or mcs['y'] > mll['y'].max() or mcs['y'] < mll['y'].min():
            #         continue
            #     else:
            #         print('Point not found but should be in!')
            #         continue

            # blob number

            nb = blobs[point]
            # if we find a 0 instead of a blob, continue
            if not nb[0]:
                continue

            isblob = np.where(blobs == nb)

            if isblob[0].size < 2500:
                print('Ooops blob too small? This should not happen')
                continue

            # lat lons of complete blob
            blats = md['lat'].values[isblob]
            blons = md['lon'].values[isblob]

            # msg indices of complete blob
            my = mll['y'][isblob]
            mx = mll['x'][isblob]

            blatmin, blatmax = blats.min(), blats.max()
            blonmin, blonmax = blons.min(), blons.max()


            # whole blob must be inside TRMM. Attention: This draws a rectangle.
            # There is still a chance that blob is not in TRMM. Checked later!

            if not (td['lon'].values.min() < blonmin) & (td['lon'].values.max() > blonmax):
                continue
            if not (td['lat'].values.min() < blatmin) & (td['lat'].values.max() > blatmax):
                continue

            ll_trmm = tm_utils.ll_toMSG(td['lon'].values, td['lat'].values)

            tx = ll_trmm['x']
            ty = ll_trmm['y']

            mpair = ua.unique_of_pair(mx, my)
            tpair = ua.unique_of_pair(tx, ty)
            #Do we need to do it that way?

            inter = np.in1d(tpair, mpair)  # returns false and true, whole grid
            inter_rev = np.in1d(mpair, tpair.flat[inter])  # Attention: this leaves out meteosat cells where no closest TRMM cell (since TRMM is coarser!)

            # have at least 500 pixels shared for MCS between TRMM and MSG
            if sum(inter) < 500:
                continue


            print(_y, _m, _d, _h, _mi)

            bprcp = td['p'].values.flat[inter]
            bflags = td['flags'].values.flat[inter]
            mtt = md['t'].values[isblob].flat[inter_rev]

            # we need same number of TRMM and MSG per plot to do the masking
            if not bprcp.size == mtt.size:
                print('Tprcp and MSGT not same, someting wrong!')
                continue


            # rtest = np.copy(td['p'].values)  # check the TRMM pixels identified
            # rtest.flat[inter] = 1500  # np.where(inter)
            #
            #
            # maskr = np.zeros_like(md['t'].values)
            # maskr[isblob] = 1000
            # # np.where(maskr>999)
            #
            # mxinter = np.in1d(mxy, mpair[inter_rev])
            # maskrr = np.zeros_like(md['t'].values)
            # maskrr.flat[mxinter] = 1100
            #
            # plt.figure()
            # ax = plt.axes(projection=ccrs.PlateCarree())
            #
            # plt.contourf(mlon, mlat, maskr,
            #              transform=ccrs.PlateCarree())  # green, MSG blob
            # plt.contourf(td['lon'].values, td['lat'].values, rtest, levels=np.arange(1300, 1600, 100),
            #              transform=ccrs.PlateCarree())  # identified TRMM pixel
            # #Identified MSG temperatures, problem: only nearest to TRMM, omits MSG pixels
            # plt.contourf(mlon, mlat, maskrr, levels=np.arange(1097, 1099, 1),
            #              transform=ccrs.PlateCarree())  # green, MSG blob
            # ax.coastlines()

            if np.count_nonzero(bprcp) < 50:
                continue

            mask = tm_utils.getTRMMconv(bflags)  # filter for convective rain
            mask = np.array(mask)

            smask = tm_utils.getTRMMstrat(bflags)  # filter for convective rain
            smask = np.array(smask)

            nz_bprcp = np.sum(bprcp>0.1)

            tall = np.nanmean(mtt[np.isfinite(bprcp)])

            # remove all these zero rainfall from blob
            bprcpNZ = bprcp[bprcp>0.1]
            mttNZ = mtt[bprcp>0.1]
            flagsNZ = bflags[bprcp>0.1]
            maskNZ = tm_utils.getTRMMconv(flagsNZ)  # list of 0 and 1, flattened!
            smaskNZ = tm_utils.getTRMMstrat(flagsNZ)  # list of 0 and 1, flattened!

            if sum(maskNZ) < 2:
                continue
            datess.append((_y, _m, _d, _h, _mi, ba, td['lon'].values.min(), td['lon'].values.max(),
                           td['lat'].values.min(), td['lat'].values.max(), blonmin, blonmax, blatmin, blatmax))

            pm = np.nanmean(bprcpNZ)
            tm = np.nanmean(mttNZ)

            ppm = np.percentile(bprcpNZ, 98)
            pmax = np.nanmax(bprcp)
            pi = float(np.sum(bprcpNZ>30)) / float(bprcpNZ.size)

            mdic['p'].append(pm)  # prcp mean of every MCS (no zero)
            mdic['pp'].append(ppm)  # value of 98 percentile of MCS (no zero)
            mdic['rain'].append(bprcpNZ) # whole rainfall field, no sum
            mdic['pmax'].append(pmax)  # maximum pcp in MCS
            mdic['pi'].append(pi)  # share of > 30mmh pixel of > 0 pixel
            mdic['t'].append(tm)  # T where PCP > 0 and overlap
            mdic['tall'].append(tall)  # T where cloud and TRMM valid (incl 0 rain)
            mdic['hod'].append(_h)  # hour of day for image
            mdic['yr'].append(_y)  # year for image
            mdic['mon'].append(_m)  # month for image
            mdic['lat'].append(lat)
            mdic['lon'].append(lon)
            mdic['tpixel_nzero'].append(np.count_nonzero(bprcp))  # nb pixel of MCS for PCP > 0
            mdic['tpixel'].append(bprcp.size)  # nb pixel of MCS including 0
            mdic['tpixel_conv'].append(sum(mask))  # number convective pixel
            mdic['tpixel_strat'].append(sum(smask))  # number stratiform pixel
            mdic['tpixel_zero'].append(np.size(bprcp) - np.count_nonzero(bprcp))  # number zero pixel
            mdic['tpixel_derived'].append(sum(mask) + sum(smask) + (np.size(bprcp) - np.count_nonzero(bprcp)))
            mdic['twhole'].append(bt)
            mdic['area'].append(isblob[0].size)

            print('Passed flag filter')

            # check for at least 500 TRMM pixels in MSG above 0 rain
            # if np.count_nonzero(bprcp) < 500:
            #    continue


            pc = np.nanmean(bprcpNZ.flat[np.where(maskNZ)])
            tc = np.nanmean(mttNZ.flat[np.where(maskNZ)])
            pic = float(np.greater(bprcpNZ.flat[np.where(maskNZ)], 30.).sum()) / float(sum(maskNZ))
            ppc = np.percentile(bprcpNZ.flat[np.where(maskNZ)], 98)
            pmaxc = bprcpNZ.flat[np.where(maskNZ)].max()

            #  print 'Nb', nb
            mdic_f['pconv'].append(pc)
            mdic_f['piconv'].append(pic)
            mdic_f['ppconv'].append(ppc)
            mdic_f['pmaxconv'].append(pmaxc)
            mdic_f['tconv'].append(tc)
            mdic_f['tnfconv'].append(tm)
            mdic_f['hod'].append(_h)
            mdic_f['yr'].append(_y)
            mdic_f['mon'].append(_m)
            mdic_f['lat'].append(lat)
            mdic_f['lon'].append(lon)
            mdic_f['tpixel_convNZ'].append(sum(maskNZ))
            mdic_f['tpixel_stratNZ'].append(sum(smaskNZ))
            cnt = cnt + 1
            print(cnt)


    myDicts = [mdic, mdic_f]
    for d in datess: print(d)

    pkl.dump(myDicts, open('/users/global/cornkle/data/OBS/test/c_paper_rainfield.p',
                           'wb'))  # MSG_TRMM_temp_pcp_300px'+str(yrange[0])+'-'+str(yrange[-1])+'_new.p', 'wb'))
    print(
        'Saved ' + 'MSG_TRMM_temp_pcp_' + str(YRANGE[0]) + '-' + str(YRANGE[-1]) + '_new.p with ' + str(
            cnt) + ' MCSs')


if __name__ == "__main__":
    tm_overlap_blobs()










