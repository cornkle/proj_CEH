import pyproj
import matplotlib
import matplotlib.pyplot as plt

print(pyproj.__version__)
print(matplotlib.__version__)

wgs84 = pyproj.Proj(proj='latlong', datum='WGS84')
fig = plt.figure()

srs = '+units=m +proj=lcc +lat_1=29.0 +lat_2=29.0 +lat_0=29.0 +lon_0=89.8'

proj_out = pyproj.Proj("+init=EPSG:4326", preserve_units=True)
proj_in = pyproj.Proj(srs, preserve_units=True)

lon, lat = pyproj.transform(proj_in, proj_out, -2235000, -2235000)

print(lon)