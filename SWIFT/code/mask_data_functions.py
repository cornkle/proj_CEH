
import os
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.crs as ccrs
from cartopy.feature import BORDERS, COASTLINE
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def regular_plot(file_lon, file_lat, dat_in,dat_count, title_in,file_out=None):
#def regular_plot(file_lon, file_lat, file_in, file_out=None):

    lon_m, fid = read_var(file_lon, "LON")
    lat_m, fid = read_var(file_lat, "LAT")


    # Extract a subsection of the full SEVIRI disc.
    bbox = (-20., 20., 0., 20.)

    #var_ss, lon_ss, lat_ss = subset_var(lai_in, lon_m, lat_m, bbox)
    var_ss, lon_ss, lat_ss = subset_var(dat_in, lon_m, lat_m, bbox)
    var_ss=dat_in

    #bbox = (-1.5, 3.5,14., 18.)
    bbox = (-2, 2,14., 18.)
    var_ss, lon, lat = subset_var(dat_in, lon_ss, lat_ss, bbox)
    var_ss_count, lon_ss, lat_ss = subset_var(dat_count, lon_ss, lat_ss, bbox)

    # Plotting.
    #cmap = cm.get_cmap("RdYlBu_r")
    #cmap = cm.get_cmap("copper")
    cmap = cm.get_cmap("coolwarm")
    plot_proj = ccrs.PlateCarree()
    F = plt.figure(figsize=(8, 8))
    A = F.add_subplot(211, projection=plot_proj)
    P = A.pcolormesh(lon, lat, var_ss, cmap=cmap, vmin=-12, vmax=0)
    A.set_extent(bbox, crs=plot_proj)
    #add_boundaries(A)
    add_gridlines(A)
    plt.colorbar(P, ax=A,extend="both")

    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
#    plt.title("LST @ 201906040915")
    plt.title(title_in)


    A = F.add_subplot(212, projection=plot_proj)
    P = A.pcolormesh(lon_ss, lat_ss, var_ss_count, cmap="cool", vmin=1, vmax=40)
    A.set_extent(bbox, crs=plot_proj)
    #add_boundaries(A)
    add_gridlines(A)
    plt.colorbar(P, ax=A,extend="max")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title("Number of contributing images")


    if file_out==None:
        plt.show()
    else:
        plt.savefig(file_out, dpi=100, bbox_inches="tight")
        plt.close(F)

    return



def subset_var(var, lon, lat, bbox):
    """Returns subsections of the input arrays.  When the input arrays are on
    curvilinear grids (e.g., SEVIRI) it's not possible to return a rectangular
    array subsection that follows lon/lat bounds exactly.  Instead this returns
    the var[col0:col1, row0:row1] subsection that contains all pixels that are
    within the requested region, which is likely to also include pixels outside
    of the requested lon/lat region.

    Subsection given by bbox = (lon0, lon1, lat0, lat1).
    """

    mask = ((bbox[0] <= lon) & (lon < bbox[1]) &
            (bbox[2] <= lat) & (lat < bbox[3]))
    c, r = np.where(mask.filled(False))
    cs, rs = slice(min(c), max(c)+1), slice(min(r), max(r))

    return var[cs, rs], lon[cs, rs], lat[cs, rs]


def add_gridlines(ax):
    ggl = ax.gridlines(color="0.5", linestyle="--", draw_labels=True)
    ggl.xlabels_top = False
    ggl.ylabels_right = False
    ggl.xformatter = LONGITUDE_FORMATTER
    ggl.yformatter = LATITUDE_FORMATTER
    return


def add_boundaries(ax):
    for feat in [COASTLINE, BORDERS]:
        #feat.scale = "10m"
        ax.add_feature(feat, alpha=0.8, zorder=100)
    return

def read_var(filename, varname, print_attrs=False):
    fid = h5py.File(filename, "r")
    var = fid[varname]

    mdi = var.attrs["MISSING_VALUE"]
    scale, offset = var.attrs["SCALING_FACTOR"], var.attrs["OFFSET"]

    var_m = np.ma.masked_equal(var, mdi)
    var_m = (var_m - offset) / scale

    return var_m, fid

def plot_stats_page(lon_ss,lat_ss,bbox,data_mean,data_stdev,data_count,file_out):
    
    cmap = cm.get_cmap("RdYlBu_r")
    plot_proj = ccrs.PlateCarree()
    #Plot of mean
    F = plt.figure(figsize=(2, 6))
    A = F.add_subplot(311, projection=plot_proj)
    P = A.pcolormesh(lon_ss, lat_ss, data_mean, cmap=cmap, vmin=0, vmax=20)
    A.set_extent(bbox, crs=plot_proj)
    #add_boundaries(A)
    add_gridlines(A)
    plt.colorbar(P, ax=A)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title("Mean")

    #Plot of stdev

    A = F.add_subplot(312, projection=plot_proj)
    P = A.pcolormesh(lon_ss, lat_ss, data_stdev, cmap=cmap, vmin=0, vmax=10)
    A.set_extent(bbox, crs=plot_proj)
    #add_boundaries(A)
    add_gridlines(A)
    plt.colorbar(P, ax=A)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title("Standard deviation")

    #plot counts       
    A = F.add_subplot(313, projection=plot_proj)
    P = A.pcolormesh(lon_ss, lat_ss,data_count, cmap=cmap, vmin=0)
    A.set_extent(bbox, crs=plot_proj)
    #add_boundaries(A)
    add_gridlines(A)
    plt.colorbar(P, ax=A)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title("Number of contributing points")

    plt.tight_layout()

    if file_out is None:
        plt.show()
    else:
        plt.savefig(file_out, dpi=100, bbox_inches="tight")
        plt.close(F)

    return

def plot_simple_stats_page(data_mean,data_stdev,data_count,use_title):    
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams.update({'font.size': 8})

    F = plt.figure(figsize=(7, 4))
    cmap = cm.get_cmap("viridis")
    plt.subplot(2,2,1)
    plt.imshow(data_mean/100,vmin=np.nanpercentile(data_mean/100,5),vmax=np.nanpercentile(data_mean/100,95),cmap=cmap)
    plt.title("Mean",pad=2)
    plt.colorbar(shrink=0.8,extend="both")
    plt.subplot(2,2,2)
    plt.imshow(data_stdev/100,vmin=np.nanpercentile(data_stdev/100,5),vmax=np.nanpercentile(data_stdev/100,95),cmap=cmap)
    plt.colorbar(shrink=0.8,extend="both")
    plt.title("Standard deviation",pad=2)
    plt.subplot(2,2,3)
    data_count=data_count.astype(float)
    data_count[np.where(np.isnan(data_mean))]=np.nan
    plt.imshow(data_count,vmin=0,cmap=cmap)
    plt.colorbar(shrink=0.8,extend="both")
    plt.title("Number of contributing points",pad=2)
    plt.subplot(2,2,4)
    cmap = cm.get_cmap("PuOr")
    plt.imshow(data_mean/data_stdev,vmax=3,vmin=-3,cmap=cmap)
    plt.colorbar(shrink=0.8,extend="both")
    plt.title("Mean/Standard deviation",pad=2)

    plt.suptitle(use_title,x=0.5,y=0.995) 
    plt.tight_layout()

def output_mean_stdev(data_all,plotdir,outdir,use_dates,dates_all,time,time_m1,nday,nrow,ncol,source_name,outstart,years,diff=True):

    from matplotlib.backends.backend_pdf import PdfPages
    if diff==True:
        pp = PdfPages(plotdir+"MeanStdev_"+outstart+"_"+time+"-"+time_m1+'.pdf')
    else:
        pp = PdfPages(plotdir+"MeanStdev_"+outstart+"_"+time+'.pdf')

    for i_mask in range(len(dates_all)):
        if i_mask<(nday)+1:
            i_min=0
            i_max=i_mask+nday+1
        elif (len(dates_all)-i_mask-1)<nday:
            i_min=i_mask-nday
            i_max=len(dates_all)
        else:  
            i_min=i_mask-nday
            i_max=i_mask+nday+1

        print(i_min,i_max)
        data_mean=np.nanmean(data_all[i_min:i_max,...],axis=(0,1))
        data_stdev=np.sqrt(np.nanvar(data_all[i_min:i_max,...],axis=(0,1)))
        data_count=np.sum(~np.isnan(data_all[i_min:i_max,...]),axis=(0,1))
        if not (os.path.exists(outdir+use_dates[i_mask][:2]+"/")):
                os.mkdir(outdir+use_dates[i_mask][:2]+"/")

        #write the output file
        if diff==True:
            fout=h5py.File(outdir+use_dates[i_mask][:2]+"/"+outstart+use_dates[i_mask]+"_"+time+"-"+time_m1,"w")  #create output file
        else:
            fout=h5py.File(outdir+use_dates[i_mask][:2]+"/"+outstart+use_dates[i_mask]+"_"+time,"w")  #create output file

        fout.create_dataset("Mean",data=data_mean)
        fout.create_dataset("Stdev",data=data_stdev)
        fout.create_dataset("count",data=data_count)

        for key in fout.keys():
            fout[key].attrs.create("START_YEAR",years[0])
            fout[key].attrs.create("END_YEAR",years[-1])
            fout[key].attrs.create("TIME",time)
            fout[key].attrs.create("N_DATES",2*nday+1)
            fout[key].attrs.create("DATA_SOURCE",source_name)
            fout[key].attrs.create("N_COLS",ncol)
            fout[key].attrs.create("N_LINES",nrow)
            fout[key].attrs.create("SCALING_FACTOR",100.0)
            fout[key].attrs.create("OFFSET",0.0)
            fout[key].attrs.create("UNITS","Degrees Celsium")

        fout.close()

        #make the plot
      #  plot_stats_page(lon_ss, lat_ss,bbox,data_mean,data_stdev,data_count,None)
        if diff==True:
            plot_simple_stats_page(data_mean,data_stdev,data_count,use_dates[i_mask][2:]+"/"+use_dates[i_mask][0:2]+" "+time+"-"+time_m1)
        else:
            plot_simple_stats_page(data_mean,data_stdev,data_count,use_dates[i_mask][2:]+"/"+use_dates[i_mask][0:2]+" "+time)
        pp.savefig()

    pp.close()
    return


def quick_test_plot(testfile,count_threshs):
    fid=h5py.File(testfile)
    data_mean=np.ma.masked_where(fid["count"][...]<count_thresh,fid["Mean"][...])
    data_stdev=np.ma.masked_where(fid["count"][...]<count_thresh,fid["Stdev"][...])
    plt.imshow(data_mean/data_stdev,vmin=-3,vmax=3,cmap="PuOr")
    plt.title(testfile+" \n CountThresh:"+str(count_thresh))
    plt.colorbar()

def quick_test_mean(testfile,count_thresh):    
    cmap = cm.get_cmap("viridis")
    fid=h5py.File(testfile)
    data_mean=np.ma.masked_where(fid["count"][...]<count_thresh,fid["Mean"][...])
    plt.imshow(data_mean/100,vmin=np.nanpercentile(fid["Mean"][...]/100,5),vmax=np.nanpercentile(fid["Mean"][...]/100,95),cmap=cmap)
    plt.title("Mean",pad=2)
    plt.colorbar(shrink=0.8,extend="both")
    plt.title(testfile+" \n CountThresh:"+str(count_thresh))



