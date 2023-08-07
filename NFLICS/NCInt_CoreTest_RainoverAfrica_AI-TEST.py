import xarray as xr
from saveCore_standalone_NFLICS import run_powerBlobs, powerBlob_utils as utils, util as wavelet
import matplotlib.pyplot as plt
import pandas as pd
from utils import constants as cnst, u_grid, u_interpolate as u_int
import glob
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
import numpy as np
import os


BOX = [-20,30,3,27]

def draw_map(ax, data, lon, lat, title=None, mask_sig=None, quiver=None, contour=None, cbar_label=None, **kwargs):
    mapp = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), **kwargs)  # this is the actual plot

    ## mask for significance indicator
    if mask_sig is not None:
        plt.contourf(lon, lat, mask_sig, colors='none', hatches='.',
                     levels=[0.5, 1], linewidth=0.1)

    ## quiver list
    if quiver is not None:
        qu = ax.quiver(quiver['x'], quiver['y'], quiver['u'], quiver['v'], scale=quiver['scale'])
    ## additional contour on plot
    if contour is not None:
        ax.contour(contour['x'], contour['y'], contour['data'], levels=contour['levels'], cmap=contour['cmap'])

    ax.coastlines()  ## adds coastlines
    # Gridlines
    xl = ax.gridlines(draw_labels=True);  # adds latlon grid lines
    xl.top_labels = False  ## labels off
    xl.right_labels = False
    plt.title(title)
    # Countries
    #   ax.add_feature(cartopy.feature.BORDERS, linestyle='--'); # adds country borders
    cbar = plt.colorbar(mapp)  # adds colorbar
    cbar.set_label(cbar_label)



def time_to_imergMin(h,mins):
    return str(h*60+mins).zfill(4)

def calc_wavelet(msg, inds, weights, shape):

    data_test = msg['IR108_BT'].squeeze().values
    try:
        data5k = u_int.interpolate_data(data_test, inds, weights, shape)
    except IndexError:
        print('Interpolation problem, continue')


    outt, nogood, t_thresh_size, t_thresh_cut, pix_nb, area_img = utils.filter_img(data5k, 5)
    wav = wavelet.waveletT(outt, dataset='METEOSAT5K_vera')
    power_msg = utils.find_scales_dominant(wav, nogood, area_img, dataset='MSG')

    return power_msg

def calc_5kdummy():
    dummy = xr.open_dataset(glob.glob(cnst.other_drive + '/nflics/core_rainfall_test/IR_108_BT_202305/IR_108_BT_*.nc')[0], decode_times=False)
    data_resolution = 5  # in km
    # make salem grid
    grid5k = u_grid.make(np.arange(BOX[0], BOX[1]), np.arange(BOX[2], BOX[3]), data_resolution * 1000)
    dlon = dummy.lon_2d.values
    dlat = dummy.lat_2d.values
    inds, weights, shape = u_int.interpolation_weights_grid(dlon, dlat, grid5k)
    lonN, latN = grid5k.ll_coordinates

    del dummy

    return inds, weights, shape, lonN, latN

inds, weights, shape, lonN, latN = calc_5kdummy()

for dday in range(1,31):
    for hhour in range(0,24):
        for mmin in [0,30]:

            idate = '202305'+str(dday)
            hour = hhour
            mins = mmin
            msg_time = pd.to_datetime('2000-01-01'+' '+str(hour)+':'+str(mins))#+pd.Timedelta('15min')
            #tag = 'MFG'
            tag = 'MSG'

            if os.path.isfile('/home/ck/Desktop/lmcs/ncint_test/202305'+str(dday).zfill(2)+'_'+str(hour).zfill(2)+str(mmin).zfill(2)+'.jpg'):
                print('file exists, continue')
                continue

            try:
                # msgfile = sorted(glob.glob(cnst.other_drive + '/nflics/core_rainfall_test/IR_108_BT_202305/IR_108_BT_'+idate+'_'+str(msg_time.hour).zfill(2)+str(msg_time.minute).zfill(2)+'.nc'))[0]
                # rainAIfile = sorted(glob.glob(cnst.other_drive + '/nflics/rainoverafricaAI/2023/05/MSG3'+idate+'-S'+str(hour).zfill(2)+str(mins).zfill(2)+'-E*.nc'))[0]
                # #GPM_DPR = glob.glob(cnst.other_drive + '/nflics/core_rainfall_test/GPM_DPR_GMI_L2B_v7/2B.GPM*.'+idate+)
                # IMERG_lateRunfile = sorted(glob.glob(cnst.other_drive + '/nflics/core_rainfall_test/IMERG_lateRun_L3v6/3B-HHR*.'+idate+'-S*.'+time_to_imergMin(hour,mins)+'.V06*.nc4'))[0]
                msg_time = pd.to_datetime('2000-01-01' + ' ' + str(hour) + ':' + str(mins))  # +pd.Timedelta('15min')
                # tag = 'MFG'
                tag = 'MSG'
                msgfile = sorted(glob.glob(cnst.other_drive + '/nflics/core_rainfall_test/IR_108_BT_202305/IR_108_BT_' + idate + '_' + str(hour).zfill(2) + str(mins).zfill(2) + '.nc'))[0]
                rainAIfile = sorted(glob.glob(cnst.other_drive + '/nflics/rainoverafricaAI/2023/05/MSG3' + idate + '-S*-E' + str(msg_time.hour).zfill(2) + str(msg_time.minute).zfill(2) + '.nc'))[0]
                # GPM_DPR = glob.glob(cnst.other_drive + '/nflics/core_rainfall_test/GPM_DPR_GMI_L2B_v7/2B.GPM*.'+idate+)
                IMERG_lateRunfile = sorted(glob.glob(cnst.other_drive + '/nflics/core_rainfall_test/IMERG_lateRun_L3v6/3B-HHR*.' + idate + '-S*.' + time_to_imergMin(hour, mins) + '.V06*.nc4'))[0]
            except:
                continue
            print(msgfile)
            print(rainAIfile)
            print(IMERG_lateRunfile)

            # if np.array(msgfile, rainAIfile, IMERG_lateRunfile).any() == []:
            #     continue
            msg = xr.open_dataset(msgfile).squeeze()
            rainA = xr.open_dataset(rainAIfile).squeeze()
            imerg = xr.open_dataset(IMERG_lateRunfile).sel(lat=slice(BOX[2],BOX[3]), lon=slice(BOX[0],BOX[1])).squeeze()

            posi = np.where(imerg['HQprecipitation'].T.values>1)

            if np.sum(posi) == 0:
                continue

            ilats = imerg.lat[posi[0]]
            ilons = imerg.lon[posi[1]]
            # delt_lats = np.max(ilats)-np.min(ilats)
            # delt_lons = np.max(ilons)-np.min(ilons)

            power_msg = calc_wavelet(msg, inds, weights, shape)

            boxi = [np.min(ilons)-2,np.max(ilons)+2,np.min(ilats)-2,np.max(ilats)+2]
            f=plt.figure(figsize=(15,5), dpi=200)  # this opens a plot window
            ax = f.add_subplot(121, projection=ccrs.PlateCarree())  # this opens a new plot axis
            #draw_map(ax, ds['t'], ds.lon, ds.lat, levels=np.arange(-100,50), cmap='jet')
            tmin = -90
            tmax = 50
            draw_map(ax, msg['IR108_BT'], msg.lon_2d, msg.lat_2d, cmap='jet', levels=np.arange(tmin,tmax,2))
            ax.contour(rainA.longitude, rainA.latitude,  rainA['posterior_mean'], levels=[-1,1,5,10,30], cmap='Reds', linewidths=2)
            try:
                ax.contour(lonN, latN,  power_msg, levels=[-1,1,5,10,30], cmap='Blues', linewidths=2)
            except:
                print('Contour error,pass')
                continue
            ax.set_ylim(boxi[2], boxi[3])
            ax.set_xlim(boxi[0],boxi[1])
            ax.set_title('rainoverAI')

            ax = f.add_subplot(122, projection=ccrs.PlateCarree())  # this opens a new plot axis
            draw_map(ax, msg['IR108_BT'], msg.lon_2d, msg.lat_2d, cmap='jet', levels=np.arange(tmin,tmax,2))
            ax.contour(imerg.lon, imerg.lat,  imerg['HQprecipitation'].T, levels=[-1,1,10,20], cmap='Reds', linewidths=2)
            ax.set_ylim(boxi[2], boxi[3])
            ax.set_xlim(boxi[0],boxi[1])
            ax.set_title('imerg')
            plt.tight_layout()
            f.savefig('/home/ck/Desktop/lmcs/ncint_test/202305'+str(dday).zfill(2)+'_'+str(hour).zfill(2)+str(mmin).zfill(2)+'.jpg')

            plt.close('all')


            del msg
            del rainA
            del imerg

