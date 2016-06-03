#!/bin/env python

def read_deforest_data():
    
    import numpy as np
    
    # Read tree cover and loss
    varShape = (222,1384)
    varMDI = np.float32(13.5)
    coverFile = "/prj/vera/EO/MSG/hansen_treecover2000_MSGgrid.dat"
    # Read tree cover in 2000
    cover = np.fromfile(coverFile,dtype=varMDI.dtype)
    cover.shape = varShape
    cover = np.flipud(cover)

    # Read loss for each year between 2001-2014
    filePattern = "/prj/vera/EO/MSG/hansen_lossyear_year{0:d}_MSGgrid.dat"
    loss = []
    for yr in range(1,14):    
        lossFile = filePattern.format(yr)
        losstmp = np.fromfile(lossFile,dtype=varMDI.dtype)
        losstmp.shape = varShape
        losstmp = np.flipud(losstmp)
        loss.append(losstmp)
    return cover, loss


if __name__ == "__main__":

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    # Read lat/lon
    ll = np.load('/users/global/danbel/msg/tropWA/latlon.npz')
    lon = ll['lon']
    lat = ll['lat']

    # Load deforestation data
    cover, loss = read_deforest_data()
    
    # Create map projection
    LON_RANGE=[-10,-4]
    LAT_RANGE=[3.9,8.3]
    M = Basemap(projection='geos',resolution='i', lon_0=0, \
                llcrnrlat =LAT_RANGE[0],urcrnrlat=LAT_RANGE[1], \
                llcrnrlon =LON_RANGE[0],urcrnrlon=LON_RANGE[1] )
    xx,yy = M(lon,lat)

    # Plot tree cover in 2000
    F = plt.figure(figsize=(16,10.3))
    A = plt.subplot(111)
    title = "Tree cover in 2000"
    A.set_title(title,fontsize=15)
    # Plot
    IMAGE = M.contourf(xx,yy,cover,20)
    # draw coastlines
    M.drawcountries(linewidth=0.4)
    M.drawcoastlines(linewidth=0.4)
    # Draw a line around the map region.
    M.drawmapboundary()
    # draw parallels and meridians.
    parallels = np.arange(-90.,90.,1)
    M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=14)
    meridians = np.arange(-180.,180.,1)
    M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=14)
    # add colorbar
    plt.colorbar(IMAGE)
    #
    plt.tight_layout(pad=4.0)
    
    savefile = "/users/global/danbel/msg/tropWA/archive_tropWA/deforest/figures/cover2000.png"
    plt.savefig(savefile)
    #plt.show()
    plt.close(F)

    # Plot tree cover loss
    subplNo = 0
    for yr in range(1,14):
        subplNo = subplNo + 1
        F = plt.figure(figsize=(16,10.3))
        A = plt.subplot(111)
        ttl = "Tree cover loss in 20{0:02d}"
        title = ttl.format(yr)
        A.set_title(title,fontsize=15)
        # Plot
        clevs = np.arange(0,35,2)
        IMAGE = M.contourf(xx,yy,loss[subplNo-1],clevs)
        # draw coastlines
        M.drawcountries(linewidth=0.4)
        M.drawcoastlines(linewidth=0.4)
        # Draw a line around the map region.
        M.drawmapboundary()
        # draw parallels and meridians.
        parallels = np.arange(-90.,90.,1)
        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=14)
        meridians = np.arange(-180.,180.,1)
        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=14)
        # add colorbar
        plt.colorbar(IMAGE)
        #
        plt.tight_layout(pad=4.0)

        savePattern = "/users/global/danbel/msg/tropWA/archive_tropWA/deforest/figures/loss20{0:02d}.png"
        savefile = savePattern.format(yr)
        plt.savefig(savefile)
        #plt.show()
        plt.close(F)

