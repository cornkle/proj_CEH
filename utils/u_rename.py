import os
import glob

str_replace = 'cores_-40'
str_new = 'cores_MSG_-40'
dirpath = '/prj/vera/cores/'

files = glob.glob(dirpath+str_replace+'*.nc')

for f in files:
    if str_replace in f:
        out = f.replace(str_replace, str_new)
        os.rename(f, out)
    
