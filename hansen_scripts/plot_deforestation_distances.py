import matplotlib
#matplotlib.use('PS')
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from netCDF4 import Dataset
import numpy as np
from matplotlib.patches import Rectangle
from plot_utils import binned_cmap


def plot_distance_from_ocean(percent_threshold=0.):
    lats = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    lons = np.arange(-180, 180, 0.05) + 0.5 * 0.05
    lon_bounds = np.hstack((lons - 0.5*0.05, np.array([lons[-1]+0.5*0.05])))
    lat_bounds = np.hstack((lats - 0.5*0.05, np.array([lats[-1]+0.5*0.05])))
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/distance_from_ocean_0pt05deg.nc', 'r') as data:
        distance = data.variables['distance_from_ocean'][:]
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc', 'r') as data:
        defn_pc = data.variables['forest_cover_loss'][:]
    distance[defn_pc <= percent_threshold] = np.nan
    projection = ccrs.PlateCarree()
    fig = plt.figure(figsize=(15, 5)) 
    ax = plt.axes(projection=projection)
    ax.set_extent((-180, 180, -90, 90), crs=projection)
    contours = [0, 1, 2, 3, 4, 5, 10, 20, 50, 100]
    p = plt.pcolormesh(lon_bounds, lat_bounds, distance, cmap=cm.viridis)
    ax.coastlines(color='black', linewidth=1)
    cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.12, ax.get_position().width, 0.03])
    cbar = fig.colorbar(p, orientation='horizontal', cax=cax, aspect=40, pad=0.12)
    cbar.ax.set_xlabel('distance from coast (km)', fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    ax.set_xticks(np.arange(-180, 181, 90), crs=projection)
    ax.set_yticks(np.arange(-90, 91, 90), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=16)
    ax.tick_params(axis='x', pad=5)
    plt.savefig(f'../figures/distance_from_coast_pc{int(percent_threshold)}', dpi=400, bbox_inches='tight', pad_inches=0.1)
    plt.show()


def distance_histogram(lon_west=-180, lon_east=180, lat_south=-90, lat_north=90, percent_threshold=0.,
                       title='', title_color='k'):
    lats = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    lons = np.arange(-180, 180, 0.05) + 0.5 * 0.05
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/distance_from_ocean_0pt05deg.nc', 'r') as data:
        distance = data.variables['distance_from_ocean'][:]
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc', 'r') as data:
        defn_pc = data.variables['forest_cover_loss'][:].astype(float)
    defn_pc[defn_pc <= percent_threshold] = np.nan
    distance[defn_pc <= percent_threshold] = np.nan
    lon_west_idx = np.where(lons > lon_west)[0][0]
    lon_east_idx = np.where(lons < lon_east)[0][-1] + 1
    lat_south_idx = np.where(lats > lat_south)[0][0]
    lat_north_idx = np.where(lats < lat_north)[0][-1] + 1
    distance_box = distance[lat_south_idx:lat_north_idx, lon_west_idx:lon_east_idx]
    defn_box = defn_pc[lat_south_idx:lat_north_idx, lon_west_idx:lon_east_idx]
    plt.hist(distance_box.ravel(), bins=100)
    plt.title(title, color=title_color)
    plt.show()
    all_distances = distance_box.ravel()
    return np.histogram(all_distances[~np.isnan(all_distances)], bins=100)


def continent_boxes():
    boxes = {
    'South America': (-82, -30, -60, 12),
    'North America': (-170, -50, 10, 75),
    'Africa': (-20, 55, -40, 37),
    'Europe': (-20, 175, 37, 80),
    'Asia': (60, 175, -10, 37),
    'Oceania': (110, 180, -50, -10)
    }
    return boxes


def continent_colors():
    colors = {
    'South America': 'green',
    'North America': 'blue',
    'Africa': 'orange',
    'Europe': 'cyan',
    'Asia': 'black',
    'Oceania': 'red'
    }
    return colors


def plot_deforestation_maps():
    lats = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    lons = np.arange(-180, 180, 0.05) + 0.5 * 0.05

    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc', 'r') as data:
        defn_pc = data.variables['forest_cover_loss'][:]

    continents = continent_boxes()
    colors = continent_colors()

    projection = ccrs.PlateCarree()
    fig = plt.figure(figsize=(15, 5)) 
    ax = plt.axes(projection=projection)
    ax.set_extent((-180, 180, -90, 90), crs=projection)
    p = plt.pcolormesh(lons, lats, defn_pc, vmin=0., cmap=cm.gist_heat_r)
    ax.coastlines(color='black', linewidth=1)
    cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.12, ax.get_position().width, 0.03])
    cbar = fig.colorbar(p, orientation='horizontal', extend='min', cax=cax, aspect=40, pad=0.12)
    cbar.ax.set_xlabel('forest cover loss (%)', fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    ax.set_xticks(np.arange(-180, 181, 90), crs=projection)
    ax.set_yticks(np.arange(-90, 91, 90), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=16)
    ax.tick_params(axis='x', pad=5)
    plt.savefig('../figures/forest_cover_loss_percent', dpi=600, bbox_inches='tight', pad_inches=0.1)

    for continent, box in continents.items():
        box_west = box[0]
        box_east = box[1]
        box_south = box[2]
        box_north = box[3]
        ax.add_patch(Rectangle(xy=(box_west, box_south), width=box_east-box_west,
                               height=box_north-box_south, transform=projection,
                               facecolor='none', edgecolor=colors[continent], linewidth=2))

    plt.savefig('../figures/forest_cover_loss_percent_continents', dpi=600, bbox_inches='tight', pad_inches=0.1)
    plt.show()


def continent_distance_histograms(percent_threshold=0.):
    continents = continent_boxes()
    colors = continent_colors()
    for continent, box in continents.items():
        distance_histogram(lon_west=box[0], lon_east=box[1], lat_south=box[2], lat_north=[3],
                           percent_threshold=percent_threshold, title=f'{continent} deforestation > {percent_threshold}%',
                           title_color=colors[continent])


def region_distances(lon_west=-180, lon_east=180, lat_south=-90, lat_north=90, percent_threshold=0.,
                     forest_threshold=0.):
    lats = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    lons = np.arange(-180, 180, 0.05) + 0.5 * 0.05
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/distance_from_ocean_0pt05deg.nc', 'r') as data:
        distance = data.variables['distance_from_ocean'][:].data
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc', 'r') as data:
        defn_pc = data.variables['forest_cover_loss'][:].data.astype(float)
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/hansen_0pt05/global_percent_cover_2000_0pt05deg.nc', 'r') as data:
        cover2000 = data.variables['tree_cover_2000'][:].data.astype(float)
    distance[defn_pc <= percent_threshold] = np.nan
    distance[cover2000 <= forest_threshold] = np.nan
    lon_west_idx = np.where(lons > lon_west)[0][0]
    lon_east_idx = np.where(lons < lon_east)[0][-1] + 1
    lat_south_idx = np.where(lats > lat_south)[0][0]
    lat_north_idx = np.where(lats < lat_north)[0][-1] + 1
    distance_box = distance[lat_south_idx:lat_north_idx, lon_west_idx:lon_east_idx]
    all_distances = distance_box.ravel()
    valid_distances = all_distances[~np.isnan(all_distances)]
    return valid_distances


def global_distance_histogram(percent_threshold=0., cumulative=False):

    global_distances = region_distances(lon_west=-180, lon_east=180,
                                        lat_south=-90, lat_north=90,
                                        percent_threshold=percent_threshold)
    max_dist = np.max([global_distances.max(), global_distances.max()])
    bins = np.arange(0, (max_dist//20 + 2) * 20, 20)

    global_histogram, bin_edges = np.histogram(global_distances, bins=bins, density=True)
    if cumulative:
        global_histogram = np.cumsum(global_histogram)*20

    bin_centres = bin_edges[:-1] + (bin_edges[1] - bin_edges[0])/2.
    plt.figure(figsize=(9, 4.5))
    ax = plt.gca()
    plt.plot(bin_centres, global_histogram, '-o', color='k', linewidth=2, ms=4)
    ax.tick_params(labelsize=14)
    plt.xlabel('distance from coast (km)', fontsize=16)
    ax.set_xlim([0, bin_edges[-1]])
    plt.title(f'Forest cover loss > {int(percent_threshold)}%', fontsize=14)
    if cumulative:
        plt.ylabel('fraction of deforestation \n within distance', fontsize=16)
        plt.gca().set_ylim(bottom=0.)
        plt.tight_layout()
        plt.savefig(f'../figures/global_histogram_pc{int(percent_threshold)}_cumulative', dpi=400)
    else:
        plt.ylabel('pdf', fontsize=16)
        plt.tight_layout()
        plt.savefig(f'../figures/global_histogram_pc{int(percent_threshold)}')
    plt.show()


def west_africa_histogram(percent_threshold=0., cumulative=False):

    west_africa_distances = region_distances(lon_west=-15, lon_east=10,
                                        lat_south=4, lat_north=12,
                                        percent_threshold=percent_threshold)
    max_dist = west_africa_distances.max()
    bins = np.arange(0, (max_dist//20 + 2) * 20, 20)

    west_africa_histogram, bin_edges = np.histogram(west_africa_distances, bins=bins, density=True)
    if cumulative:
        west_africa_histogram = np.cumsum(west_africa_histogram)*20

    bin_centres = bin_edges[:-1] + (bin_edges[1] - bin_edges[0])/2.
    plt.figure(figsize=(9, 4.5))
    ax = plt.gca()
    plt.plot(bin_centres, west_africa_histogram, '-o', color='k', linewidth=2, ms=4)
    ax.tick_params(labelsize=14)
    plt.xlabel('distance from coast (km)', fontsize=16)
    ax.set_xlim([0, bin_edges[-1]])
    plt.title(f'Forest cover loss > {int(percent_threshold)}% - West Africa', fontsize=14)
    if cumulative:
        plt.ylabel('fraction of deforestation \n within distance', fontsize=16)
        plt.gca().set_ylim(bottom=0.)
        plt.tight_layout()
        plt.savefig(f'../figures/west_africa_histogram_pc{int(percent_threshold)}_cumulative', dpi=400)
    else:
        plt.ylabel('pdf', fontsize=16)
        plt.tight_layout()
        plt.savefig(f'../figures/west_africa_histogram_pc{int(percent_threshold)}')
    plt.show()


def compare_amazonia_west_africa(percent_threshold=0.):
    amazonia_box = (-80, -43, -18, 11)
    west_africa_box = (-15, 10, 4, 12)

    lats = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    lons = np.arange(-180, 180, 0.05) + 0.5 * 0.05
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc', 'r') as data:
        defn_pc = data.variables['forest_cover_loss'][:].astype(float)
    defn_pc[defn_pc <= percent_threshold] = np.nan

    projection = ccrs.PlateCarree()
    fig = plt.figure(figsize=(15, 5)) 
    ax = plt.axes(projection=projection)
    ax.set_extent((-180, 180, -90, 90), crs=projection)
    p = plt.pcolormesh(lons, lats, defn_pc, cmap=cm.gist_heat_r)
    ax.coastlines(color='black', linewidth=1)
    cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.12, ax.get_position().width, 0.03])
    cbar = fig.colorbar(p, orientation='horizontal', cax=cax, aspect=40, pad=0.12)
    cbar.ax.set_xlabel('forest cover loss (%)', fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    ax.set_xticks(np.arange(-180, 181, 90), crs=projection)
    ax.set_yticks(np.arange(-90, 91, 90), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=16)
    ax.tick_params(axis='x', pad=5)

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][0:2]
    for i, box in enumerate([west_africa_box, amazonia_box]):
        box_west = box[0]
        box_east = box[1]
        box_south = box[2]
        box_north = box[3]
        ax.add_patch(Rectangle(xy=(box_west, box_south), width=box_east-box_west,
                               height=box_north-box_south, transform=projection,
                               facecolor='none', edgecolor=colors[i], linewidth=2))

    plt.savefig('../figures/forest_cover_loss_west_africa_amazonia_boxes', dpi=600, bbox_inches='tight', pad_inches=0.1)


    west_africa_distances= region_distances(lon_west=-15, lon_east=10,
                                            lat_south=4, lat_north=12,
                                            percent_threshold=percent_threshold)
    amazonia_distances = region_distances(lon_west=-80, lon_east=-43,
                                          lat_south=-18, lat_north=11,
                                          percent_threshold=percent_threshold)
    max_dist = np.max([west_africa_distances.max(), amazonia_distances.max()])
    bins = np.arange(0, (max_dist//20 + 2) * 20, 20)

    west_africa_histogram, bin_edges = np.histogram(west_africa_distances, bins=bins, density=True)
    amazonia_histogram, bin_edges = np.histogram(amazonia_distances, bins=bins, density=True)

    bin_centres = bin_edges[:-1] + (bin_edges[1] - bin_edges[0])/2.
    plt.figure(figsize=(9, 4.5))
    ax = plt.gca()
    plt.plot(bin_centres, west_africa_histogram, '-o', label='West Africa', linewidth=2, ms=4)
    plt.plot(bin_centres, amazonia_histogram, '-v', label='Amazonia', linewidth=2, ms=4)
    plt.legend(loc='best', fontsize=16)
    ax.tick_params(labelsize=14)
    plt.xlabel('distance from coast (km)', fontsize=16)
    ax.set_xlim([0, bin_edges[-1]])
    plt.ylabel('pdf', fontsize=16)
    plt.title(f'Forest cover loss > {int(percent_threshold)}%', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'../figures/west_africa_amazonia_histograms_pc{int(percent_threshold)}')
    plt.show()


def compare_africa_south_america(percent_threshold=0.):
    africa_box = (-20, 55, -40, 37)
    south_america_box = (-82, -30, -60, 13)

    lats = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    lons = np.arange(-180, 180, 0.05) + 0.5 * 0.05
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc', 'r') as data:
        defn_pc = data.variables['forest_cover_loss'][:].astype(float)
    defn_pc[defn_pc <= percent_threshold] = np.nan

    if percent_threshold == 0.:
        projection = ccrs.PlateCarree()
        fig = plt.figure(figsize=(15, 5)) 
        ax = plt.axes(projection=projection)
        ax.set_extent((-180, 180, -90, 90), crs=projection)
        p = plt.pcolormesh(lons, lats, defn_pc, cmap=cm.gist_heat_r)
        ax.coastlines(color='black', linewidth=1)
        cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.12, ax.get_position().width, 0.03])
        cbar = fig.colorbar(p, orientation='horizontal', cax=cax, aspect=40, pad=0.12)
        cbar.ax.set_xlabel('forest cover loss (%)', fontsize=14)
        cbar.ax.tick_params(labelsize=14)
        ax.set_xticks(np.arange(-180, 181, 90), crs=projection)
        ax.set_yticks(np.arange(-90, 91, 90), crs=projection)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.tick_params(labelsize=16)
        ax.tick_params(axis='x', pad=5)

        colors = ['#6fe377', 'purple']
        for i, box in enumerate([africa_box, south_america_box]):
            box_west = box[0]
            box_east = box[1]
            box_south = box[2]
            box_north = box[3]
            ax.add_patch(Rectangle(xy=(box_west, box_south), width=box_east-box_west,
                                   height=box_north-box_south, transform=projection,
                                   facecolor='none', edgecolor=colors[i], linewidth=2))

        plt.savefig('../figures/forest_cover_loss_africa_south_america_boxes', dpi=600, bbox_inches='tight', pad_inches=0.1)

    africa_distances= region_distances(lon_west=-20, lon_east=55,
                                       lat_south=-40, lat_north=37,
                                       percent_threshold=percent_threshold)
    south_america_distances = region_distances(lon_west=-82, lon_east=-30,
                                               lat_south=-60, lat_north=12,
                                               percent_threshold=percent_threshold)
    max_dist = np.max([africa_distances.max(), south_america_distances.max()])
    bins = np.arange(0, (max_dist//20 + 2) * 20, 20)

    africa_histogram, bin_edges = np.histogram(africa_distances, bins=bins, density=True)
    south_america_histogram, bin_edges = np.histogram(south_america_distances, bins=bins, density=True)

    bin_centres = bin_edges[:-1] + (bin_edges[1] - bin_edges[0])/2.
    plt.figure(figsize=(9, 4.5))
    ax = plt.gca()
    plt.plot(bin_centres, africa_histogram, '-o', label='Africa', linewidth=2, ms=4, color='#6fe377')
    plt.plot(bin_centres, south_america_histogram, '-v', label='South America', linewidth=2, ms=4, color='purple')
    plt.legend(loc='best', fontsize=16)
    ax.tick_params(labelsize=14)
    plt.xlabel('distance from coast (km)', fontsize=16)
    ax.set_xlim([0, bin_edges[-1]])
    plt.ylabel('pdf', fontsize=16)
    plt.title(f'Forest cover loss > {int(percent_threshold)}%', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'../figures/africa_south_america_histograms_pc{int(percent_threshold)}')
    plt.show()


def plot_west_africa_deforestation():
    lats = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    lons = np.arange(-180, 180, 0.05) + 0.5 * 0.05
    lon_bounds = np.hstack((lons - 0.5*0.05, np.array([lons[-1]+0.5*0.05])))
    lat_bounds = np.hstack((lats - 0.5*0.05, np.array([lats[-1]+0.5*0.05])))
    with Dataset("/users/global/bethar/python/deforestation-coast-distance/data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc", 'r') as data:
        defn_pc = data.variables['forest_cover_loss'][:]
    west_africa_box = (-15, 10, 4, 12)
    contours = [1, 10, 20]
    cmap, norm = binned_cmap(contours, 'Reds', extend='both', fix_colours=[(-1, 'white')])
    cmap.set_bad('gray', alpha=0)
    projection = ccrs.PlateCarree()
    fig = plt.figure(figsize=(15, 5)) 
    ax = plt.axes(projection=projection)
    ax.set_extent((-15, 10, 4, 12), crs=projection)
    p = plt.pcolormesh(lon_bounds, lat_bounds, defn_pc, cmap=cmap, norm=norm, transform=projection)
    ax.coastlines(color='black', linewidth=1)
    cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.12, ax.get_position().width, 0.03])
    cbar = fig.colorbar(p, orientation='horizontal', cax=cax, aspect=40, pad=0.12, extend='both')
    cbar.set_ticks(contours)
    cbar.ax.set_xlabel('forest cover loss (%)', fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    ax.set_xticks(np.arange(-15, 11, 5), crs=projection)
    ax.set_yticks(np.arange(4, 13, 4), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=16)
    ax.tick_params(axis='x', pad=5)
    plt.savefig('../figures/forest_cover_loss_percent_west_africa', dpi=400, bbox_inches='tight', pad_inches=0.1)
    plt.show()


def compare_tropical_continents(percent_threshold=0., cumulative=True):
    africa_box = (-25, 53, -30, 30)
    america_box = (-120, -25, -30, 30)
    asia_box = (53, 180, -30, 30)
    all_tropics_box = (-180, 180, -30, 30)

    lats = np.arange(-90, 90, 0.05) + 0.5 * 0.05
    lons = np.arange(-180, 180, 0.05) + 0.5 * 0.05
    with Dataset('/users/global/bethar/python/deforestation-coast-distance/data/hansen_0pt05/global_percent_forest_loss_0pt05deg.nc', 'r') as data:
        defn_pc = data.variables['forest_cover_loss'][:].astype(float)
    lon_bounds = np.hstack((lons - 0.5*0.25, np.array([lons[-1]+0.5*0.25])))
    lat_bounds = np.hstack((lats - 0.5*0.25, np.array([lats[-1]+0.5*0.25])))
    defn_pc[defn_pc <= percent_threshold] = np.nan

    tropical_continent_colours = ['#1f78b4', '#ff7f00', '#33a02c']

    if percent_threshold == 0.:
        projection = ccrs.PlateCarree()
        fig = plt.figure(figsize=(10, 5)) 
        ax = plt.axes(projection=projection)
        ax.set_extent((-180, 180, -30, 30), crs=projection)
        p = plt.pcolormesh(lon_bounds, lat_bounds, defn_pc, cmap=cm.gist_heat_r, vmin=0., vmax=100., rasterized=True)
        ax.coastlines(color='black', linewidth=1)
        cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.12, ax.get_position().width, 0.03])
        cbar = fig.colorbar(p, orientation='horizontal', cax=cax, aspect=40, pad=0.12)
        cbar.ax.set_xlabel('forest cover loss (%)', fontsize=14)
        cbar.set_ticks(np.arange(0, 101, 20))
        cbar.ax.tick_params(labelsize=14)
        ax.set_xticks(np.arange(-180, 181, 90), crs=projection)
        ax.set_yticks(np.arange(-30, 31, 30), crs=projection)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.tick_params(labelsize=16)
        ax.tick_params(axis='x', pad=5)

        for i, box in enumerate([america_box, africa_box, asia_box]):
            box_west = box[0]
            box_east = box[1]
            box_south = box[2]
            box_north = box[3]
            ax.add_patch(Rectangle(xy=(box_west, box_south), width=box_east-box_west,
                                   height=box_north-box_south, transform=projection,
                                   facecolor='none', edgecolor=tropical_continent_colours[i], linewidth=2, zorder=5, clip_on=False))

        ax.set_zorder(0)
        plt.gcf().text(0.06, 0.68, '(a)', fontsize=16)
        plt.savefig('../figures/forest_cover_loss_tropical_continent_boxes.eps', dpi=200, bbox_inches='tight', pad_inches=0.1)
        plt.savefig('../figures/forest_cover_loss_tropical_continent_boxes.png', dpi=200, bbox_inches='tight', pad_inches=0.1)

    all_tropics_distances = region_distances(lon_west=all_tropics_box[0], lon_east=all_tropics_box[1],
                                             lat_south=all_tropics_box[2], lat_north=all_tropics_box[3],
                                             percent_threshold=percent_threshold)
    america_distances = region_distances(lon_west=america_box[0], lon_east=america_box[1],
                                         lat_south=america_box[2], lat_north=america_box[3],
                                         percent_threshold=percent_threshold)
    africa_distances= region_distances(lon_west=africa_box[0], lon_east=africa_box[1],
                                       lat_south=africa_box[2], lat_north=africa_box[3],
                                       percent_threshold=percent_threshold)
    asia_distances= region_distances(lon_west=asia_box[0], lon_east=asia_box[1],
                                     lat_south=asia_box[2], lat_north=asia_box[3],
                                     percent_threshold=percent_threshold)
    
    max_dist = np.max([all_tropics_distances.max(), africa_distances.max(),
                       america_distances.max(), asia_distances.max()])
    bin_size = 25
    bins = np.arange(0, (max_dist//bin_size + 2) * bin_size, bin_size)

    all_tropics_histogram, bin_edges = np.histogram(all_tropics_distances, bins=bins, density=True)
    africa_histogram, bin_edges = np.histogram(africa_distances, bins=bins, density=True)
    america_histogram, bin_edges = np.histogram(america_distances, bins=bins, density=True)
    asia_histogram, bin_edges = np.histogram(asia_distances, bins=bins, density=True)

    if cumulative:
        all_tropics_histogram = np.cumsum(all_tropics_histogram)*bin_size
        africa_histogram = np.cumsum(africa_histogram)*bin_size
        america_histogram = np.cumsum(america_histogram)*bin_size
        asia_histogram = np.cumsum(asia_histogram)*bin_size

    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    plt.plot(bin_edges[1:], all_tropics_histogram, 'k-o', label='All tropics', linewidth=1.5, ms=3.5, zorder=5)
    plt.plot(bin_edges[1:], america_histogram, '-v', color=tropical_continent_colours[0], label='Tropical America', linewidth=1.5, ms=4)
    plt.plot(bin_edges[1:], africa_histogram, '-d', color=tropical_continent_colours[1], label='Tropical Africa', linewidth=1.5, ms=4)
    plt.plot(bin_edges[1:], asia_histogram, '-s', color=tropical_continent_colours[2], label='Tropical Asia/Australasia', linewidth=1.5, ms=3)
    plt.legend(loc='best', fontsize=16)
    ax.tick_params(labelsize=14)
    plt.xlabel('distance from coast (km)', fontsize=16)
    ax.set_xlim([0, bin_edges[-1]])
    plt.title(f'Forest loss > {int(percent_threshold)}%', fontsize=14)
    if cumulative:
        plt.ylabel('fraction of deforested pixels \n within distance', fontsize=16)
        plt.gca().set_ylim(bottom=0.)
        plt.gca().axvline(50, color='gray', alpha=0.5)
        plt.gca().axvline(300, color='gray', alpha=0.5)
        plt.tight_layout()
        plt.gcf().text(0.04, 0.96, '(b)', fontsize=16)
        plt.savefig(f'../figures/tropical_continents_histograms_pc{int(percent_threshold)}_cumulative.eps')
        plt.savefig(f'../figures/tropical_continents_histograms_pc{int(percent_threshold)}_cumulative.png', dpi=300)
    else:
        plt.ylabel('pdf', fontsize=16)
        plt.tight_layout()
        plt.savefig(f'../figures/tropical_continents_histograms_pc{int(percent_threshold)}.eps')
    plt.show()


if __name__ == '__main__':
    compare_tropical_continents(percent_threshold=0.)
    compare_tropical_continents(percent_threshold=20., cumulative=True)
