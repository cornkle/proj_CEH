import numpy as np
import pandas as pd
import datetime as dt
from eod import trmm, msg, tm_utils
import ipdb
from collections import defaultdict



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
    dic = defaultdict(list)

    mlon = mObj.lon
    mlat = mObj.lat

    mll = tm_utils.ll_toMSG(mlon, mlat)

    cnt = 0

    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered
    for _y, _m, _d, _h, _mi in zip(dates.y, dates.m, dates.d, dates.h, dates.mi):

        # set zero shift time for msg
        date = dt.datetime(_y, _m, _d, _h, _mi)

        dt0 = tm_utils.minute_delta(_mi, 30)

        #time difference max + 10min
        if abs(dt0) > 5:
            continue

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
        if big.empty:
            continue

        td = tObj.get_ddata(_y, _m, _d, _h, _mi, cut=[4, 20])
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
            #print(mll['x'].min(), mll['y'].min(), mll['x'].max(), mll['y'].max() )
            #print(mcs['x'],mcs['y'])
            if not all(point):
                if mcs['x'] > mll['x'].max() or mcs['x'] < mll['x'].min() or mcs['y'] > mll['y'].max() or mcs['y'] < mll['y'].min():
                    continue
                else:
                    print('Point not found but should be in!')
                    continue

            # blob number

            nb = blobs[point]
            # if we find a 0 instead of a blob, continue
            if not nb[0]:
                yy,xx = point
                #ipdb.set_trace()
                if blobs[yy-1:yy+1, xx-1:xx+1].any():
                    sub = blobs[yy-1:yy+1, xx-1:xx+1]
                    nb = sub[sub>0][0]
                    print('Found new ind in blob ', nb)
                else:
                    print('Found only 0 ind for blob, continue')

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
            mpair = (mx + my) * (mx + my + 1) / 2 + my

            blatmin, blatmax = blats.min(), blats.max()
            blonmin, blonmax = blons.min(), blons.max()

            ll_trmm = tm_utils.ll_toMSG(td['lon'].values, td['lat'].values)

            tx = ll_trmm['x']
            ty = ll_trmm['y']

            tpair = (tx + ty) * (tx + ty + 1) / 2 + ty
            inter = np.in1d(tpair, mpair)  # returns false and true, whole grid
            inter_rev = np.in1d(mpair, tpair.flat[inter])  # Attention: this leaves out meteosat cells where no closest TRMM cell (since TRMM is coarser!)

            # have at least 500 pixels shared for MCS between TRMM and MSG
            if sum(inter) < 500:
                continue

            bprcp = td['p'].values.flat[inter]
            bflags = td['flags'].values.flat[inter]
            mtt = md.values[isblob].flat[inter_rev]

            # we need same number of TRMM and MSG per plot to do the masking
            if not bprcp.size == mtt.size:
                print('Tprcp and MSGT not same, someting wrong!')
                continue

            zmask = tm_utils.getTRMMconv(bflags)  # filter for convective rain
            zmask = np.array(zmask)

            nz_bprcp = np.count_nonzero(bprcp)

            if np.count_nonzero(bprcp) < 50:
                continue

            # remove all these zero rainfall from blob
            bprcp = bprcp.flat[np.where(bprcp)]
            bflags = bflags.flat[np.where(bprcp)]
            mtt = mtt.flat[np.where(bprcp)]

            pm = np.nanmean(bprcp)
            tm = np.nanmean(mtt)
            tmin = np.nanmin(mtt)
            ppm = np.percentile(bprcp, 98)
            pmax = np.nanmax(bprcp)

            mask = tm_utils.getTRMMconv(bflags)  # filter for convective rain
            mask = np.array(mask)

            if sum(mask) < 2:
                continue

            mdic['pmean'].append(pm)
            mdic['pmax'].append(pmax)
            mdic['p98'].append(ppm)
            mdic['tmean'].append(tm)
            mdic['tmean_tab'].append(bt)
            mdic['area_tab'].append(ba)
            mdic['tmin'].append(tmin)
            mdic['hod'].append(_h)
            mdic['yr'].append(_y)
            mdic['mon'].append(_m)
            mdic['lat'].append(lat)
            mdic['lon'].append(lon)
            mdic['p_nb_conv'].append(sum(mask))
            mdic['p_znb_conv'].append(sum(zmask))
            mdic['p_znb'].append(nz_bprcp)
            mdic['p_nb_nz'].append(bprcp.size)
            mdic['p_nb_gt30'].append(np.sum(bprcp>30))
            dic['all_p_nz'].extend(bprcp)
            dic['all_t_nz'].extend(mtt)

            df = pd.DataFrame(data=mdic)
            df2 = pd.DataFrame(data=dic)

            cnt += 1

    df.to_pickle('/users/global/cornkle/C_paper/MCS_blobs/df_MCS.pkl')
    df2.to_pickle('/users/global/cornkle/C_paper/MCS_blobs/df_allPT_inMCS.pkl')
    print('Saved '+ str(cnt) + ' MCSs')


if __name__ == "__main__":
    tm_overlap_blobs()










