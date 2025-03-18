"""
Download data from the Hansen et al. (2013) Global Forest Change dataset
v1.7 (2000-2019).

See https://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.7.html
and the original journal article doi.org/10.1126/science.1244693
for details of the dataset.

This script can be run as
>> python download_hansen_tiles.py /my/example/directory/ dataset_name
where /my/example/directory/ is the directory to which the files will be
downloaded and dataset_name is the data layer (e.g. gain).
This will download all the 10x10 degree tiles of data globally.

Alternatively, this script can be imported as a module and the
download_tile() function used to download individual 10x10 degree tiles,
or the download_all_tiles() function used to download a complete global layer.

Bethan L. Harris, UK Centre for Ecology & Hydrology, 15th April 2021.
"""


import sys, os
from urllib.request import urlopen
from shutil import copyfileobj
from tqdm import tqdm


def download_tile(download_directory, filename):
    """
    Download a 10x10 degree tile of data from the Hansen et al. v1.7 dataset.

    The download will be skipped if a file with that name already exists
    in the download directory. The filename will be printed out if the download fails.
    Parameters:
    download_directory (str): directory to which files will be saved.
    filename (str): filename of data to be downloaded, e.g. 'Hansen_GFC-2019-v1.7_treecover2000_40N_080W.tif'.
    Returns:
    None
    """

    save_filename = f'{download_directory}/{filename}'
    if not os.path.exists(save_filename):
        url = f'http://storage.googleapis.com/earthenginepartners-hansen/GFC-2019-v1.7/{filename}'
        try:
            with urlopen(url) as tile, open(save_filename, 'wb') as save_file:
                copyfileobj(tile, save_file)
        except:
            print(f'{filename} download failed')


def download_all_tiles(download_directory, dataset_name):
    """
    Download all 10x10 degree tiles globally for a given data field from the Hansen et al. dataset.

    An error will be returned if the name of the data field is not one found in the data catalogue.
    Parameters:
    download_directory (str): directory to which files will be saved.
    dataset_name (str): name of data layer (e.g. 'treecover2000').
    Returns:
    None
    """

    valid_datasets = ['gain', 'lossyear', 'datamask', 'treecover2000', 'first', 'last']
    if not dataset_name in valid_datasets:
        raise KeyError(f'{dataset_name} is not a valid data field. Valid fields are {", ".join(valid_datasets)}.')
    if not os.path.isdir(download_directory):
        os.makedirs(download_directory)
    lons = [str(western_lon).zfill(3) + "W" for western_lon in range(180, 0, -10)]
    lons.extend([str(eastern_lon).zfill(3) + "E" for eastern_lon in range(0, 180, 10)])
    lats = [str(southern_lat).zfill(2) + "S" for southern_lat in range(50, 0, -10)]
    lats.extend([str(northern_lat).zfill(2) + "N" for northern_lat in range(0, 90, 10)])
    all_filenames = [f'Hansen_GFC-2019-v1.7_{dataset_name}_{lat}_{lon}.tif' for lon in lons for lat in lats]
    progress_format = '{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}'
    for filename in tqdm(all_filenames, desc=dataset_name, bar_format=progress_format):
        download_tile(download_directory, filename)

 
if __name__ == "__main__":
    # Read the download directory and data layer name from command line arguments
    # and download all available 10x10 degree tiles for that layer. 
    download_directory, dataset_name = sys.argv[1:]
    download_all_tiles(download_directory, dataset_name)


