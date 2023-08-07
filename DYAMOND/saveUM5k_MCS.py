# -*- coding: utf-8 -*-


import numpy as np
from scipy.ndimage.measurements import label
import xarray as xr
import os
import ipdb

import glob
from scipy import ndimage
from utils import u_interpolate as u_int


HOD = range(24)  # hours of day
#YRANGE = range(2004, 2015)

def olr_to_bt(olr):
    sigma = 5.670373e-8
    return ((olr/sigma)**0.25)-273.15


def saveMCS():

    ffiles = '/home/ck/DIR/cornkle/data/DYAMOND/UM-5km/SAm/'
    vars = ['rlut', 'pr']

    metum_res = 5

    fnames = glob.glob(ffiles+'/rlut/*.nc')

    windgrid = glob.glob(ffiles+'/ua/*.nc')
    vgrid = glob.glob(ffiles + '/va/*.nc')

    srfc_dummy = xr.open_dataset(fnames[0])
    pl_dummy = xr.open_dataset(windgrid[5])
    plv_dummy = xr.open_dataset(vgrid[5])

    inds, weights, shape = u_int.interpolation_weights(pl_dummy.longitude, pl_dummy.latitude, srfc_dummy.longitude,
                                                       srfc_dummy.latitude)

    indsv, weightsv, shapev = u_int.interpolation_weights(plv_dummy.longitude, plv_dummy.latitude, srfc_dummy.longitude,
                                                       srfc_dummy.latitude)

    cnt = 0
    for f in fnames:

        orl = (xr.open_dataset(f)['rlut']).load()
        vals = orl.values
        nvals = olr_to_bt(vals)
        orl.values = nvals


        bname = os.path.basename(f)
        strdate = bname[-17:-9]
        #ipdb.set_trace()
        pfile = glob.glob(ffiles+'/pr/*-'+strdate+'*.nc')
        pcp = (xr.open_dataset(pfile[0])['pr'])*3600

        pfile = glob.glob(ffiles + '/prw/*-'+strdate+'*.nc')
        try:
            tcw = (xr.open_dataset(pfile[0])['prw'])
        except:
            continue
        pfile = glob.glob(ffiles + '/ua/*-'+strdate+'*.nc')
        uu = (xr.open_dataset(pfile[0])['ua'])
        udiff = uu.sel(model_level_number=24) - uu.sel(model_level_number=3)
        pfile = glob.glob(ffiles + '/va/*-'+strdate+'*.nc')
        vv = (xr.open_dataset(pfile[0])['va'])
        vdiff = vv.sel(model_level_number=24) - vv.sel(model_level_number=3)

        pfile = glob.glob(ffiles + '/hfss/*-'+strdate+'*.nc')
        sh = (xr.open_dataset(pfile[0])['hfss'])
        pfile = glob.glob(ffiles + '/hfls/*-'+strdate+'*.nc')
        lh = (xr.open_dataset(pfile[0])['hfls'])
        ef = lh / (sh+lh)
        ef.values[ef.values>1] = np.nan
        ef.values[ef.values<0] = np.nan
        w = (tcw[tcw['time.hour'] == 15]).mean('time').squeeze()
        eff = (ef[ef['time.hour'] == 15]).mean('time').squeeze()
        u12 = (udiff[udiff['time.hour'] == 15]).squeeze()
        v12 = (vdiff[vdiff['time.hour'] == 15]).squeeze()
        inds, weights, shape = u_int.interpolation_weights(u12.longitude, u12.latitude, srfc_dummy.longitude,
                                                           srfc_dummy.latitude)
        indsv, weightsv, shapev = u_int.interpolation_weights(v12.longitude, v12.latitude, srfc_dummy.longitude,
                                                           srfc_dummy.latitude)
        try:
            u12_regrid = u_int.interpolate_data(u12.values, inds, weights, shape)
        except:
            ipdb.set_trace()
            print('Regrid failed, continue')
            continue
        v12_regrid = u_int.interpolate_data(v12.values, indsv, weightsv, shapev)

        shear_wg = np.sqrt(u12_regrid**2+v12_regrid**2)


        #ipdb.set_trace()
        shear = w.copy(deep=True)
        shear.values = shear_wg
        shear.name = 'ua'

        lon, lat = np.meshgrid(orl.longitude, orl.latitude)

        for p,t in zip(pcp, orl):


            # ipdb.set_trace()
            #
            # _y = t['time.year'].values[0]
            # _m = t['time.month'].values[0]
            _h = int(t['time.hour'].values)
            #
            # _d = t['time.day'].values[0]
            _mi = int(t['time.minute'].values)
            #

            # date = dt.datetime(_y, _m, _d, _h, 0)
            t.values = ndimage.gaussian_filter(t.values, 2, mode='nearest')
            #ipdb.set_trace()
            t.values[t.values >= -40] = 0  # T threshold -10
            labels, numL = label(t.values)

            u, inv = np.unique(labels, return_inverse=True)
            n = np.bincount(inv)

            goodinds = u[n >= 200]  # defines minimum MCS size e.g. 350 km2 = 39 pix at 3x3km res , 5000km2 = 556 pixel
            print(goodinds)
            if not sum(goodinds) > 0:
                continue

            for gi in goodinds:
                if gi == 0:  # index 0 is always background, ignore!
                    continue

                inds = np.where(labels == gi)

                tmin = np.min(t.values[inds])
                if tmin > -50:
                    continue

                # cut a box for every single blob from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
                try:
                    latmax, latmin = np.max(lat[inds]), np.min(lat[inds])
                except:
                    ipdb.set_trace()
                lonmax, lonmin = np.max(lon[inds]), np.min(lon[inds])

                tout = t.sel(latitude=slice(latmin,latmax), longitude=slice(lonmin,lonmax))
                pout = p.sel(latitude=slice(latmin,latmax), longitude=slice(lonmin,lonmax))
                wout = w.sel(latitude=slice(latmin, latmax), longitude=slice(lonmin, lonmax))
                efout = eff.sel(latitude=slice(latmin, latmax), longitude=slice(lonmin, lonmax))
                shearout = shear.sel(latitude=slice(latmin, latmax), longitude=slice(lonmin, lonmax))

                ds = xr.Dataset()
                ds['tir'] = tout
                ds['tcw'] = wout
                ds['prcp'] = pout
                ds['ef'] = efout
                ds['shear'] = shearout

                #ipdb.set_trace()


                #
                # if (sum(mmask.flatten())*25 < 350) | (outt.max()>250):# or (sum(mmask.flatten())*25 > 1500000): #or (outt.max()<0.1)
                #     continue
                #
                # if sum(mask2.flatten()) < 5:  # sum(mmask.flatten())*0.3:
                #     print('Kickout: TRMM MSG overlap less than 3pix of cloud area')
                #     continue
                #
                # print('Hit:', gi)
                #
                # da = xr.Dataset({'p': (['x', 'y'], outt),
                #                  't_lag0': (['x', 'y'], dummy),
                #                  'tc_lag0': (['x', 'y'], outl),
                #                  },
                #                 coords={'lon': (['x', 'y'], lon),
                #                         'lat': (['x', 'y'], lat),
                #                         'time': date})

                ds.attrs['meanT'] = np.mean(t.values[inds])
                ds.attrs['minT'] = np.min(t.values[inds])
                ds.attrs['meanP'] = np.mean(p.values[inds])
                ds.attrs['maxP'] = np.max(p.values[inds])
                ds.attrs['area'] = inds[0].size

                savefile = '/home/ck/DIR/cornkle/data/DYAMOND/UM-5km/MCS/SAm/UM-5km/' + strdate + '_'+str(_h).zfill(2)+str(_mi).zfill(2)+'_' + str(gi) + '.nc'
                # try:
                #     os.remove(savefile)
                # except OSError:
                #     print('OSError, no dir?')
                #     pass
                print('Hour', ds['time.hour'])
                ds.to_netcdf(path=savefile, mode='w')
                print('Saved ' + savefile)


        print('Saved ' + str(cnt) + ' MCSs as netcdf.')
