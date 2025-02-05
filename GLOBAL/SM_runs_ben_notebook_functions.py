import iris
import xarray as xr
import pandas as pd
import glob
import os
import datetime
import glob
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import numpy as np
import cartopy.crs as ccrs
import scipy.stats as stats
import statsmodels.api as sm
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
import matplotlib.patches as mpatch
from scipy.ndimage import uniform_filter
import warnings
warnings.filterwarnings("ignore")

def olr_to_bt(olr):
    #Application of Stefan-Boltzmann law
    sigma = 5.670373e-8
    tf = (olr/sigma)**0.25
    #Convert from bb to empirical BT
    a = 1.228
    b = -1.106e-3
    Tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return Tb - 273.15

def empbt(tf):
    tf=tf+273.15
    #Convert Stefan-Boltzmann law derived brightness T to empirical
    a = 1.228
    b = -1.106e-3
    Tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return Tb - 273.15

def get_Sahel(MCS_data):
    try:
        MCS_data=MCS_data[(MCS_data['clat']< 19) & (MCS_data['clat']> 9) & (MCS_data['clon']< 12) & (MCS_data['clon']> -12)]
    except:
        MCS_data=MCS_data[(MCS_data['lat']< 19) & (MCS_data['lat']> 9) & (MCS_data['lon']< 12) & (MCS_data['lon']> -12)]
    return MCS_data

def bin_2d(MCS_data,moist_var="tcw",shear_var="ushear650_850",rain_var="rain_max",t_var="tmin",max_pcle=0.99,min_pcle=0.01,grid_spec=10,grid_spec2=None):
    if grid_spec2 is None:
        grid_specq = grid_spec
    else:
        grid_specq = grid_spec2
        
    q_offset,shr_offset=MCS_data[moist_var].quantile(min_pcle), MCS_data[shear_var].quantile(min_pcle)
    q_bins=np.linspace(q_offset,MCS_data[moist_var].quantile(max_pcle),grid_specq)
    shr_bins=np.linspace(shr_offset,MCS_data[shear_var].quantile(max_pcle),grid_spec)
    q_bin=q_bins[1]-q_bins[0]
    shr_bin=shr_bins[1]-shr_bins[0]
    
    rain_vals=xr.DataArray(np.zeros((grid_specq,grid_spec)),dims=["q","shear"],coords=[q_bins,shr_bins])
    ctt_vals=rain_vals.copy()
    counts=rain_vals.copy()
    area_vals=rain_vals.copy()
    for i in range(0,grid_specq):
        for j in range(0,grid_spec):
            q,shr = q_bins[i],shr_bins[j]
            bin_mcs = MCS_data[(MCS_data[moist_var] > q) & (MCS_data[moist_var] <= q+q_bin) & 
                            (MCS_data[shear_var] > shr) & (MCS_data[shear_var] <= shr+shr_bin)]
            counts[i,j] = len(bin_mcs)
            if len(bin_mcs) > 0:
                rain_vals[i,j] = bin_mcs[rain_var].mean()
                ctt_vals[i,j] = bin_mcs[t_var].mean()
                area_vals[i,j] = bin_mcs['area'].mean()
            else:
                rain_vals[i,j] = np.NaN
                ctt_vals[i,j] = np.NaN
                area_vals[i,j] = np.NaN

    return rain_vals, ctt_vals, area_vals, counts

def plot_hist(MCS_data,region,model,MCS_var,unit="",cmap="viridis",nice_name="",moist_var="tcw",shear_var="ushear650_850",ax=None,vmax=None,cut_off=5,grid_spec=10,save=False):      
    try:
        plev=int(shear_var[-3:])
        try:
            plev=np.max(int(moist_var[-3:]),plev)
        except:
            pass
        print(len(MCS_data))
        MCS_data=MCS_data[MCS_data["psfc"]>plev]
        print(len(MCS_data))
    except:
        pass

    hist, ctt_vals, area_vals, counts = bin_2d(MCS_data,rain_var=MCS_var,moist_var=moist_var,shear_var=shear_var,t_var="tmin",grid_spec=grid_spec)
    hist=hist.where(counts>=cut_off)
    if hist.max() - hist.min() > hist.max() and cmap=="viridis":
        cmap="RdBu"
    
    hours=MCS_data["hour"].unique()
    if len(hours)>1:
        tstr="{:02d}-{:02d}UTC".format(hours[0],hours[-1])
    else:
        tstr="{:02d}UTC".format(hours[0])
    
    ylab=moist_var
    if moist_var[0]=="q":
        ylab="q%s (g kg$^{-1}$)" % (moist_var[-3:])
    elif moist_var=="tcw":
        ylab="Total column water (mm)"
    try:
        xlab="u%s - u%s (m s$^{-1}$)"%(shear_var[6:9],shear_var[-3:])
    except:
        xlab=shear_var

    if len(nice_name)==0:
        nice_name=MCS_var

    if "rain" in MCS_var:
        cmap,unit,nice_name="Blues",r"mm hr$^{-1}$",MCS_var[5:]+" rain"
    elif MCS_var == "tmean":
        cmap,unit,nice_name="Reds",r"$^{\circ}$C","mean BT"
    elif MCS_var == "tmin":
        cmap,unit,nice_name="Oranges",r"$^{\circ}$C","min BT"
    elif "area" in MCS_var:
        cmap,unit,nice_name="Purples",r"km${^2}$","storm area"
    elif "Q"in MCS_var:
        unit=r"Net W m$^{-2}$"

    if ax==None:
        fig,ax=plt.subplots(figsize=(7,5))
    hist.plot(ax=ax,cmap=cmap,cbar_kwargs={"label":unit},robust=True)
    ax.set_xlabel(xlab,fontsize=12)
    ax.set_ylabel(ylab,fontsize=12)
    ax.set_box_aspect(0.95)
    ax.set_facecolor("grey")
    ax.set_title("{}, {}: {}\n{} MCSs ".format(model.replace("_"," "), region, nice_name, tstr) + r"($\bf n=%s$) "% len(MCS_data), fontsize=14)
    return ax

def bin_1d(MCS_data,MCS_var,bin_var,red_var,grid_spec=10,grid_spec2=None,cut_off=5,max_pcle=0.99,min_pcle=0.01):
    bins=np.linspace(MCS_data[bin_var].quantile(min_pcle),MCS_data[bin_var].quantile(max_pcle),grid_spec)
    MCS_data=MCS_data[(MCS_data[red_var]>MCS_data[red_var].quantile(min_pcle)) & (MCS_data[red_var]<MCS_data[red_var].quantile(max_pcle))]
    binw=bins[1]-bins[0]
    
    bin_vals=xr.DataArray(np.zeros(grid_spec),dims=["x"],coords=[bins])
    bin_errs=bin_vals.copy()
    bin_counts=bin_vals.copy()
    for i in range(0,grid_spec):
        q = bins[i]
        bin_mcs = MCS_data[(MCS_data[bin_var] > q) & (MCS_data[bin_var] <= q+binw)]
        bin_counts[i] = len(bin_mcs)
        if len(bin_mcs) >= cut_off:
            bin_vals[i] = bin_mcs[MCS_var].mean()
            bin_errs[i] = bin_mcs[MCS_var].std()/np.sqrt(len(bin_mcs))
        else:
            bin_vals[i] = np.NaN
            bin_errs[i] = np.NaN
    bin_vals=bin_vals.where(~np.isnan(bin_errs.values))
    return bin_vals, bin_errs, bin_counts

def plot_fits(MCS_data,ax,MCS_var,bin_var,red_var,rain_var="rain_max",t_var="tmin",lobf=True,grid_spec=10,cut_off=5,color="maroon"):
    y, y_err, counts = bin_1d(MCS_data,MCS_var,bin_var,red_var,grid_spec=grid_spec,cut_off=cut_off)
    y, y_err = y.dropna(dim="x"), y_err.dropna(dim="x")
    
    if MCS_var[0]=="t":
        color="peru"
    elif "rain" in MCS_var:
        color="steelblue"
    y.plot(ax=ax,marker="o",color=color)
    ax.fill_between(y["x"].values,(y-y_err).values,(y+y_err).values,color=color,alpha=0.2)
    X = y["x"].values
    fit = sm.WLS(y.values, sm.add_constant(X), weights=y_err.values).fit()
    out=fit.params[1]*X + fit.params[0]
    if lobf==True:
        ax.plot(X,out,color="darkgrey")
    r,p=stats.pearsonr(y.values,out)
    ax.set_title("Grad.={:.2f}$\pm${:.2f}; r={:.2f} ({:.2f})".format(fit.params[1],fit.bse[1],r,p),fontsize=10,loc="center")
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_box_aspect(0.75)
    return ax

def data_prep(region,model,com_box=True,paper=False,version="v0",source="jmn"):
    #Handles the ERA5 obs tables Conni has made; note slightly different naming conventions
    if model == "Obs" or model == "Obs_1yr":
        region2 = region
        if region == "Sahel":
            region2 = "WAf"

        if region2=="WAf" and version == "v1":
            MCS_data=pd.read_csv("Tables/MCS_5000km2_-50C_0.1degTIR-IMERG_hourly_WAf_ERA5.csv")
        else:
            files=glob.glob("/gws/nopw/j04/lmcs/global_MCS_tables_OBS/ERA5_added_12LT/{}/????_MCS_5000km2_-50C_0.1degTIR-IMERG_hourly_{}_ERA5.csv".format(region2,region2))
            MCS_data=[]
            for file in files:
                MCS_data.append(pd.read_csv(file))
            MCS_data=pd.concat(MCS_data,axis=0)
        
        if region2 == "WAf" or region2 == "india":
            sdate, edate = pd.Timestamp("2016-08-01"), pd.Timestamp("2016-09-09")
            #MCS_data=MCS_data[(MCS_data["utc_date"]>="2016-08-01") & (MCS_data["utc_date"]<="2016-09-09")]
            if region == "Sahel":
                MCS_data=MCS_data[(MCS_data['clat']< 19) & (MCS_data['clat']> 9) & (MCS_data['clon']< 12) & (MCS_data['clon']> -12)]
        else:
            #MCS_data=pd.read_csv("/gws/nopw/j04/lmcs/global_MCS_tables_OBS/ERA5_added_12LT/{}/2019_MCS_5000km2_-50C_0.1degTIR-IMERG_hourly_{}_ERA5.csv".format(region,region))
            #MCS_data=MCS_data[(MCS_data["utc_date"]>="2019-01-20") & (MCS_data["utc_date"]<="2019-02-28")]
            #Note these are not in the right year! Obs table missing 2020, darn. Will only affect Obs_1yr crop.
            sdate, edate = pd.Timestamp("2019-01-20"), pd.Timestamp("2019-02-28")
        MCS_data=MCS_data.rename(columns={"precipitation_max":"rain_max","precipitation_mean":"rain_mean",
                                          "ushear925_650":"ushear650_925","ushear850_650":"ushear650_850","tcwv":"tcw"})

        MCS_data["utc_date"]=pd.to_datetime(MCS_data["utc_date"])
        #Cut down obs to only cover days of the year within the DYAMOND 40 day periods
        MCS_data=MCS_data[(MCS_data["utc_date"].dt.dayofyear >= sdate.dayofyear) & (MCS_data["utc_date"].dt.dayofyear <= edate.dayofyear)]
        #Further possible restriction to reduce data to only the exact same 40 days
        if model=="Obs_1yr":
            MCS_data=MCS_data[MCS_data["utc_date"].dt.year == sdate.year]
            
    elif model == "CP4conni":
        MCS_data=pd.read_csv("/gws/nopw/j04/lmcs/WestAfrica_MCS_tables_OBS-clover/CP4hist_bulk_-50_5000km2_5-25N_12W15E_15kprecip_hourly.csv")
        MCS_data=MCS_data.rename(columns={"shear":"ushear650_925","prcp":"rain_max","tcwv":"tcw"})
        MCS_data['ushear650_925']=-1*MCS_data['ushear650_925']
        try:
            MCS_data=MCS_data[(MCS_data['year']==2016) & (MCS_data['month'].isin([8,9]))]
        except:
            MCS_data=MCS_data[MCS_data['month'].isin([8,9])]
        if region == "Wafrica":
            MCS_data=MCS_data[(MCS_data['lat']< 18) & (MCS_data['lat']> 5) & (MCS_data['lon']< 13) & (MCS_data['lon']> -17)]
        if region == "Sahel":
            MCS_data=MCS_data[(MCS_data['lat']< 19) & (MCS_data['lat']> 9) & (MCS_data['lon']< 12) & (MCS_data['lon']> -12)]
        MCS_data['q925'] = MCS_data['q925']/1000
        MCS_data[["tmin","tmean"]]=empbt(MCS_data[["tmin","tmean"]])
            
    #Following handles the DYAMOND kscale + CP4 tables I've made
    else:
        region2 = region
        if region == "Sahel":
            region2 = "Wafrica"

        if model.upper() == "CP4":
            MCS_data=pd.read_csv('Tables/CP4/Cp4_{}_MCS_table.csv'.format(region2))
        elif version == "v0":
            MCS_data=pd.read_csv('Tables/kscale_v0/{}_{}_MCS_table_extended.csv'.format(model,region2))
        elif version == "v1":
            MCS_data=pd.read_csv('Tables/kscale_v1/{}_{}_MCS_table_{}.csv'.format(model,region2,source))

        #Regions where LAM and channel domains don't entirely overlap are by default cut to common area; also extract Sahel-only storms from wider Wafrica data.
        if region == "Sahel":
            MCS_data=MCS_data[(MCS_data['clat']< 19) & (MCS_data['clat']> 9) & (MCS_data['clon']< 12) & (MCS_data['clon']> -12)]
        """
        if com_box and region == "Samerica":
            MCS_data=MCS_data[MCS_data['clat']>-30]
        if com_box and region == "India":
            MCS_data=MCS_data[MCS_data['clat']<=26]
        if com_box and region == "Aus":
            MCS_data=MCS_data[MCS_data['clat']>-18]
        """
        try:
            MCS_data=MCS_data[MCS_data['q925']<1]
        except:
            pass
        if model[:3].upper()!="CP4":
        #    MCS_data[["tmin","tmean"]]=empbt(MCS_data[["tmin","tmean"]])
            res=11.1
        else:
            MCS_data["date"]=pd.to_datetime(MCS_data["date"])#,format="%Y%m%d_%H:%M")
            sdate, edate = pd.Timestamp("2016-08-01"), pd.Timestamp("2016-09-10") # end date 1 day later than DYAMOND to account for CP4 360 day calendar
            MCS_data=MCS_data[(MCS_data["date"].dt.dayofyear >= sdate.dayofyear) & (MCS_data["date"].dt.dayofyear <= edate.dayofyear)]
            res=4.5
        MCS_data['area'] = res**2 * MCS_data['area']

    #MCS_data=MCS_data[MCS_data["hour"].isin(np.arange(16,20))]
    try:
        MCS_data['q850'] = 1000*MCS_data['q850']
        MCS_data['q925'] = 1000*MCS_data['q925']
    except:
        pass
            
    #Apply precipitation filter:
    MCS_data = MCS_data[MCS_data["rain_max"] > 1]
    #Apply filter on q925; very, very low values entering likely due to pressure level sitting below surface! Value taken from typical minimum ERA5 values
    #try:
    #    MCS_data = MCS_data[MCS_data["q925"] >= 8]
    #except:
    #    pass
    return MCS_data

def spatial_binning(MCS_data,region="Wafrica",var="count",res=1):
    if region=="Wafrica":
        xbin = np.arange(-17, 26,res) 
        ybin = np.arange(5,25,res) 
    elif region=="Sahel":
        xbin = np.arange(-12,12,res)
        ybin = np.arange(9,19,res)
    elif region == "India":
        xbin = np.arange(70, 91,res) 
        ybin = np.arange(5,31,res) 
    plotbins = [xbin,ybin]
    dist = stats.binned_statistic_2d(MCS_data["clon"],MCS_data["clat"],MCS_data["tcw"],"count",bins=plotbins)
    da=xr.DataArray(dist.statistic,coords=[xbin[:-1],ybin[:-1]],dims=["longitude","latitude"])
    if var!="count":
        if var[-4:] == "norm":
            dist_var = stats.binned_statistic_2d(MCS_data["clon"],MCS_data["clat"],MCS_data[var[:-5]],"mean",bins=plotbins)
            dist_var = dist_var.statistic/dist.statistic
        else:
            dist_var = stats.binned_statistic_2d(MCS_data["clon"],MCS_data["clat"],MCS_data[var],"mean",bins=plotbins).statistic
        da_var = xr.DataArray(dist_var,coords=[xbin[:-1],ybin[:-1]],dims=["longitude","latitude"])
        da = da_var.where(da>0)
    
    #counts=xr.DataArray(dist.statistic,coords=[xbin[:-1],ybin[:-1]],dims=["longitude","latitude"])
    return da

def MCS_dist_plotter(MCS_data,sim,region="Wafrica",var="count",cmap="viridis",vmax=0,vmin=0,dropna=False,save=False):
    if dropna==True and var!="count":
        MCS_data=MCS_data[["clat","clon","tcw",var]].dropna(axis=0)

    da=spatial_binning(MCS_data,var=var)
    if var=="count":
        vmax,cmap=25,"bone_r"
        if sim=="Obs":
            vmax=da.quantile(.99)
    if vmax==0 and vmin==0:
        vmax=da.quantile(.99)
        vmin=da.quantile(.01)
    
    fig,ax=plt.subplots(figsize=(9,6),subplot_kw=dict(projection=ccrs.PlateCarree()))
    da.T.plot(ax=ax,cbar_kwargs={"shrink":0.6,"label":var.capitalize().replace("_"," ")},vmax=vmax,vmin=vmin,cmap=cmap)
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_patch(mpatch.Rectangle((-12,9),24,10,fill=False,lw=2,edgecolor="maroon"))
    ax.set_facecolor("lightgrey")
    #ax.add_patch(mpatch.Rectangle((-12,15),24,7,fill=False,lw=2,edgecolor="maroon"))
    #ax.add_patch(mpatch.Rectangle((-12,8),24,7,fill=False,lw=2,edgecolor="maroon"))
    gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True,linestyle='--')
    gl.xlines = False
    gl.ylines = False
    gl.top_labels = False
    gl.right_labels = False
    #gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120, 180])
    #gl.ylocator = mticker.FixedLocator([-20, 0, 20])
    if region=="Wafrica":
        gl.ylocator = mticker.FixedLocator([5,10,15,20])
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    ax.set_title("{}, {}. MCS (n={}) spatial distribution".format(sim.replace("_"," "),region,len(MCS_data)))
    plt.show()

def fss(thr, scale, X_f, X_o, anom=False):
    """Accumulate forecast-observation pairs to an FSS object.
    Parameters
    -----------
    thr: threshold vakue for your input fields
    scale: scale of evaluation 
    X_f: array_like
        Array of shape (m, n) containing the forecast field.
    X_o: array_like
        Array of shape (m, n) containing the observation field.
    """
    if len(X_f.shape) != 2 or len(X_o.shape) != 2 or X_f.shape != X_o.shape:
        message = "X_f and X_o must be two-dimensional arrays"
        message += " having the same shape"
        raise ValueError(message)

    ###I added this bit! Probably not sensible...
    if anom==True:
        X_f=np.where(np.sign(X_f) == np.sign(X_o), np.abs(X_f), np.NaN)
        X_o=np.abs(X_o)
    ###
        
    # convert nan's to a value below the threshold
    I_f = X_f.astype(float)
    X_f[~np.isfinite(X_f)] = thr - 1
    X_o = X_o.copy()
    X_o[~np.isfinite(X_o)] = thr - 1
 
    # Convert to binary fields with the given intensity threshold
    I_f = (X_f >= thr).astype(float)
    I_o = (X_o >= thr).astype(float)
 
    # Compute fractions of pixels above the threshold within a square
    # neighboring area by applying a 2D moving average to the binary fields
    if scale > 1:
        S_f = uniform_filter(I_f, size=scale, mode="constant", cval=0.0)*(scale)
        S_o = uniform_filter(I_o, size=scale, mode="constant", cval=0.0)*(scale)
    else:
        S_f = I_f
        S_o = I_o
 
    sum_obs_sq = np.nansum(S_o**2)
    sum_fct_obs = np.nansum(S_f * S_o)
    sum_fct_sq = np.nansum(S_f**2)
    numer = sum_fct_sq - 2.0 * sum_fct_obs + sum_obs_sq
    denom = sum_fct_sq + sum_obs_sq
 
    return 1.0 - numer / denom, np.count_nonzero(I_o)/I_o.size