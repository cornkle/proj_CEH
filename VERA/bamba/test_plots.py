import xarray as xr
import numpy as np
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd

outpath = '/users/global/cornkle/w2018_bamba/mini_forest/test_plots/'
datapath = '/users/global/cornkle/w2018_bamba/mini_forest/pb2014040*.nc'

def t_map_plot(datapath=datapath, outpath=outpath):

    outf = outpath
    ds = xr.open_mfdataset(datapath, concat_dim='TH1')

    ds = ds.sel(grid_latitude_t=slice(5, 6.6), grid_longitude_t=slice(-10, -5.6), TH1='2014-04-06T10:00:00',
                grid_latitude_uv=slice(5, 6.6), grid_longitude_uv=slice(-10, -5.6))

    veg = ds['forest_frac'].squeeze()

    lon = ds.grid_longitude_t
    lat = ds.grid_latitude_t

    lonuv = ds.grid_longitude_uv
    latuv = ds.grid_latitude_uv
    #for dh in range(24):

    f = plt.figure(figsize=(5,6))

    ax1 = f.add_subplot(311, projection=ccrs.PlateCarree())
    mapp = ax1.contourf(lon, lat, ds['lst']-273.15, levels=np.arange(20,36,1), cmap='inferno', transform=ccrs.PlateCarree(), extend='both')
    plt.colorbar(mapp)
    ax1.contour(lon, lat, veg, vmin=0.8, vmax=0.9, nlevels=2,  cmap='Greens_r', transform=ccrs.PlateCarree())
    plt.title('LST '+ str(pd.to_datetime(ds.TH1.values)))


    ax1 = f.add_subplot(312, projection=ccrs.PlateCarree())
    mapp = ax1.contourf(lon, lat, ds['T2'].squeeze()-273.15, levels=np.arange(20,36,1), cmap='inferno', transform=ccrs.PlateCarree(), extend='both')
    plt.colorbar(mapp)
    ax1.contour(lon, lat, veg, vmin=0.8, vmax=0.9, nlevels=2,  cmap='Greens_r', transform=ccrs.PlateCarree())
    plt.title('T2 '+ str(pd.to_datetime(ds.TH1.values)))

    ax1 = f.add_subplot(313, projection=ccrs.PlateCarree())
    mapp = ax1.contourf(lonuv, latuv, ds['T_pl'][0,:,:]-273.15, levels=np.arange(20,36,1), cmap='inferno', transform=ccrs.PlateCarree(), extend='both')
    plt.colorbar(mapp)
    ax1.contour(lon, lat, veg, vmin=0.8, vmax=0.9, nlevels=2,  cmap='Greens_r', transform=ccrs.PlateCarree())
    plt.title('T 950hPa '+ str(pd.to_datetime(ds.TH1.values)))
    plt.tight_layout()


def atmo_map_plot(datapath=datapath, outpath=outpath):

    outf = outpath
    ds = xr.open_mfdataset(datapath, concat_dim='TH1')

    ds = ds.sel(grid_latitude_t=slice(5, 6.6), grid_longitude_t=slice(-10, -5.6), TH1='2014-04-06T15:00:00',
                grid_latitude_uv=slice(5, 6.6), grid_longitude_uv=slice(-10, -5.6))

    veg = ds['forest_frac'].squeeze()

    lon = ds.grid_longitude_t
    lat = ds.grid_latitude_t

    lonuv = ds.grid_longitude_uv
    latuv = ds.grid_latitude_uv
    # for dh in range(24):

    f = plt.figure(figsize=(5, 6))

    ax1 = f.add_subplot(311, projection=ccrs.PlateCarree())
    mapp = ax1.contourf(lon, lat, ds['lst'] - 273.15, levels=np.arange(20, 36, 1), cmap='inferno',
                        transform=ccrs.PlateCarree(), extend='both')
    plt.colorbar(mapp)
    ax1.contour(lon, lat, veg, vmin=0.8, vmax=0.9, nlevels=2, cmap='Greens_r', transform=ccrs.PlateCarree())
    plt.title('LST ' + str(pd.to_datetime(ds.TH1.values)))

    ax1 = f.add_subplot(312, projection=ccrs.PlateCarree())

    mapp = ax1.contourf(lonuv, latuv, ds['u10'][0,:,:].squeeze(), levels=np.arange(-5,6,1), cmap='RdBu',
                        transform=ccrs.PlateCarree(), extend='both')
    plt.colorbar(mapp)
    ax1.contour(lon, lat, veg, vmin=0.8, vmax=0.9, nlevels=2, cmap='Greens_r', transform=ccrs.PlateCarree())
    plt.title('u10 ' + str(pd.to_datetime(ds.TH1.values)))

    ax1 = f.add_subplot(313, projection=ccrs.PlateCarree())
    mapp = ax1.contourf(lon, lat, ds['theta_pl'][0, :, :]-273.15,  levels=np.arange(28,35.1,0.2), cmap='inferno',
                        transform=ccrs.PlateCarree(), extend='both')
    plt.colorbar(mapp)
    ax1.contour(lon, lat, veg, vmin=0.8, vmax=0.9, nlevels=2, cmap='Greens_r', transform=ccrs.PlateCarree())
    plt.title('theta 950hPa ' + str(pd.to_datetime(ds.TH1.values)))
    plt.tight_layout()




def cross_t(datapath=datapath, outpath=outpath):

    outf = outpath
    dsorig = xr.open_mfdataset(datapath, concat_dim='TH1')

    #gap1 5.44-5.67
    #gap2 5.8-6.3
    for date in dsorig.TH1.values:

        ds = dsorig.sel(grid_latitude_t=slice(5.89, 6.25), grid_longitude_t=slice(-10, -5.6), TH1=date,
                    grid_latitude_uv=slice(5.89, 6.25), grid_longitude_uv=slice(-10, -5.6)).copy()

        ### quick test plot to check domain and included vegetation
        plt.figure()
        ds['forest_frac'][0,:,:].plot.contourf()
        ds = ds.mean(dim='grid_latitude_t')
        ds = ds.mean(dim='grid_latitude_uv')

        veg = ds['forest_frac'].squeeze()


        lon = ds.grid_longitude_t
        lonuv = ds.grid_longitude_uv

        f = plt.figure(figsize=(7, 8))

        ax1 = f.add_subplot(411)
        mapp = ax1.contourf(lonuv, ds.P_ECMWF.values, ds['u_pl'].squeeze(), cmap='RdBu', extend='both', levels=np.arange(-8,9,1))
        plt.gca().invert_yaxis()
        plt.title('u-wind ' + str(pd.to_datetime(date)), fontsize=9)
        #ax1.set_xlim(-9,-6)
        ax2 = ax1.twinx()
        ax2.plot(lon, veg, 'k')
        #ax2.set_xlim(-9, -6)

        #plt.title('LST ' + str(pd.to_datetime(veg.TH1.values)))


        ax3 = f.add_subplot(412)

        mapp1 = ax3.contourf(lon, ds.P_ECMWF.values[0:5], ds['theta_pl'].squeeze()[0:5,:]-273.15, cmap='inferno', extend='both', levels=np.arange(28,35.1,0.2))
        plt.gca().invert_yaxis()
        plt.title('Pot. T ' + str(pd.to_datetime(date)), fontsize=9)
        #ax3.set_xlim(-9,-6)
        ax4 = ax3.twinx()
        ax4.plot(lon, veg, 'k')
        #ax4.set_xlim(-9, -6)

        ax7 = f.add_subplot(413)

        mapp3 = ax7.contourf(lonuv, ds.P_ECMWF.values[0:5], ds['rh_pl'].squeeze()[0:5,:], cmap='viridis', extend='both', levels=np.arange(50,101,5))
        plt.gca().invert_yaxis()
        plt.title('RH ' + str(pd.to_datetime(date)), fontsize=9)
        #ax7.set_xlim(-9,-6)
        ax8 = ax7.twinx()
        ax8.plot(lon, veg, 'k')
        #ax8.set_xlim(-9, -6)


        ax5 = f.add_subplot(414)

        mapp2 = ax5.contourf(lonuv, ds.P_ECMWF.values, ds['omega_pl'].squeeze(), cmap='RdBu', extend='both', levels=np.arange(-3,3.1,0.2))
        plt.gca().invert_yaxis()
        plt.title('Omega ' + str(pd.to_datetime(date)), fontsize=9)
        #ax5.set_xlim(-9,-6)
        ax6 = ax5.twinx()
        ax6.plot(lon, veg, 'k')
        #ax6.set_xlim(-9, -6)

        plt.tight_layout()

        f.subplots_adjust(right=0.77)
        cax = f.add_axes([0.84, 0.77, 0.02, 0.18])
        cbar = f.colorbar(mapp, cax)
        cbar.ax.tick_params(labelsize=9)
        cbar.set_label('m s-1', fontsize=9)

        cax = f.add_axes([0.84, 0.54, 0.02, 0.18])
        cbar = f.colorbar(mapp1, cax)
        cbar.ax.tick_params(labelsize=9)
        cbar.set_label('C', fontsize=9)

        cax = f.add_axes([0.84, 0.3, 0.02, 0.18])
        cbar = f.colorbar(mapp3, cax)
        cbar.ax.tick_params(labelsize=9)
        cbar.set_label('C', fontsize=9)

        cax = f.add_axes([0.84, 0.05, 0.02, 0.18])
        cbar = f.colorbar(mapp2, cax)
        cbar.ax.tick_params(labelsize=9)
        cbar.set_label('m s-1', fontsize=9)

        f.savefig(outf+str(pd.to_datetime(date))+'.png')

        plt.close('all')


def timeseries(datapath=datapath, outpath=outpath):

    outf = outpath
    ds = xr.open_mfdataset(datapath, concat_dim='TH1')

    # gap1 5.44-5.67
    # gap2 5.8-6.3


    box = [(5.89, 6.25), (-10, -5.6), (5.89, 6.25), (-10, -5.62) ] # tai park2
    #box= [(5.2, 6.3), (-8.5,-7), (5.2, 6.32), (-8.51, -7) ] # tai park1
    #box = [(6.5, 8), (-10.57, -9.66), (6.5, 7.99), (-10.57, -9.68)] # liberia

    ds = ds.sel(grid_latitude_t=slice(box[0][0], box[0][1]), grid_longitude_t=slice(box[1][0], box[1][1]),
                grid_latitude_uv=slice(box[2][0], box[2][1]), grid_longitude_uv=slice(box[3][0], box[3][1]))  # uv coords slightly off to cheat with staggering


    veg = ds['forest_frac'].squeeze()
    w = ds['w_pl'].sel(P_ECMWF=600)
    rh = ds['rh_pl'].sel(P_ECMWF=950)
    q = ds['q_pl'].sel(P_ECMWF=950)
    u = ds['u_pl'].sel(P_ECMWF=950)
    t = ds['T_pl'].sel(P_ECMWF=950)

    ### quick test plot to check domain and included vegetation
    plt.figure()
    veg[0,:,:].plot.contourf()

    def pick(data, veg):

        forested_pos = np.where((veg>0.7) & (data!=0))
        nonforested_pos = np.where((veg<0.45) & (data!=0))

        forested = data.copy()
        nonforested = data.copy()

        try:
            forested[nonforested_pos] = np.nan
            nonforested[forested_pos] = np.nan
        except IndexError:
            pdb.set_trace()

        fmean = np.nanmean(np.nanmean(forested, axis=1), axis=1)
        nonmean = np.nanmean(np.nanmean(nonforested, axis=1), axis=1)

        return fmean, nonmean

    tveg, tnon = pick(ds['lst'].squeeze().values-273.15, veg.values)
    qveg, qnon = pick(ds['q2'].squeeze().values, veg.values)
    t2veg, t2non = pick(ds['T2'].squeeze().values - 273.15, veg.values)
    u10veg, u10non = pick(ds['u10'].squeeze().values, veg.values)
    w600veg, w600non = pick(w.squeeze().values, veg.values)
    q950veg, q950non = pick(q.squeeze().values, veg.values)
    rh950veg, rh950non = pick(rh.squeeze().values, veg.values)
    u950veg, u950non = pick(u.squeeze().values, veg.values)
    t950veg, t950non = pick(t.squeeze().values, veg.values)

    pdb.set_trace()

    ### proper plot
    f = plt.figure(figsize=(10, 10))

    ax1 = f.add_subplot(711)

    ax1.plot_date(ds.TH1.values, tveg, label='LST forest', linestyle='-')
    ax1.plot_date(ds.TH1.values, tnon, label='LST non-forest', linestyle='-')
    plt.legend()
    plt.ylabel('$\circ$C')

    ax1 = f.add_subplot(712)
    ax1.plot_date(ds.TH1.values, t2veg, label='T2 forest', linestyle='-')
    ax1.plot_date(ds.TH1.values, t2non, label='T2 non-forest', linestyle='-')
    plt.ylabel('$\circ$C')
    plt.legend()

    ax1 = f.add_subplot(714)
    ax1.plot_date(ds.TH1.values, qveg*1000, label='Q2 forest', linestyle='-')
    ax1.plot_date(ds.TH1.values, qnon*1000, label='Q2 non-forest', linestyle='-')
    plt.ylabel('g kg$^-1$')
    plt.legend()

    ax1 = f.add_subplot(715)
    ax1.plot_date(ds.TH1.values, q950veg*1000, label='Q950 forest', linestyle='-')
    ax1.plot_date(ds.TH1.values, q950non*1000, label='Q950 non-forest', linestyle='-')
    plt.ylabel('g kg$^-1$')
    plt.legend()

    ax1 = f.add_subplot(716)
    ax1.plot_date(ds.TH1.values, u950veg, label='u950 forest', linestyle='-')
    ax1.plot_date(ds.TH1.values, u950non, label='u950 non-forest', linestyle='-')
    plt.ylabel('m s$^-1$')
    plt.legend()

    ax1 = f.add_subplot(717)
    ax1.plot_date(ds.TH1.values, w600veg, label='w600 forest', linestyle='-')
    ax1.plot_date(ds.TH1.values, w600non, label='w600 non-forest', linestyle='-')
    plt.ylabel('m s$^-1$')
    plt.legend()

    ax1 = f.add_subplot(713)
    ax1.plot_date(ds.TH1.values, t2veg, label='T950 forest', linestyle='-')
    ax1.plot_date(ds.TH1.values, t2non, label='T950 non-forest', linestyle='-')
    plt.ylabel('$\circ$C')
    plt.legend()

    plt.tight_layout()






