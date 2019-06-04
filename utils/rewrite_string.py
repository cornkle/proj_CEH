import glob
import os
import ipdb

path = '/users/global/cornkle/shared/data/OBS/AMSRE/aqua/raw_night_AMSR2'


files = glob.glob(path + '/*_2012*.nc4')

for f in files:

    ssplit = f.split('.')
    new = ssplit[0][0:-6] + '.' +  ssplit[1]

    #ipdb.set_trace()
    os.rename(f,new)