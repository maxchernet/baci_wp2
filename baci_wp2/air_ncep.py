# Upscaling air temperature to MODIS
# Air temperature data are from NCEP reanalysis

import numpy as np
import netCDF4 as nc
import optparse
from scipy import interpolate
import os

parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(),
                                   usage=globals()['__doc__'])

parser.add_option('--tile', action="store", dest="tile",
                      type=str, help="tile")

parser.add_option('--year', action="store", dest="year",
                      type=str, help="year")

parser.add_option('--day_from', action="store", dest="day_from",
                      type=str, help="day from")

parser.add_option('--day_to', action="store", dest="day_to",
                      type=str, help="day to")

(options, args) = parser.parse_args()

year = int(options.year)
# day = int(options.day)
days = np.arange(int(options.day_from), int(options.day_to))

ds = nc.Dataset('/group_workspaces/jasmin2/baci/sigil/air_temperature/2.5_deg/air.%d.nc' % year)
ds1 = np.load('latlon/latlon_%s.npz' % options.tile)

x = np.arange(0, 1200)
y = np.arange(0, 1200)
f = interpolate.interp2d(x, y, ds1['lat'][:], kind='linear')
xnew = np.arange(0, 1200, 50)
ynew = np.arange(0, 1200, 50)
lat_low = f(xnew, ynew)

f = interpolate.interp2d(x, y, ds1['lon'][:], kind='linear')
lon_low = f(xnew, ynew)

min_lat = np.zeros((ds.variables['lat'].shape[0], 1200, 1200))
min_lat_low = np.zeros((ds.variables['lat'].shape[0], lat_low.shape[0], lat_low.shape[1]))
for i in xrange(ds.variables['lat'].shape[0]):
    min_lat[i, :, :] = (ds1['lat'] - ds.variables['lat'][i]) ** 2
    min_lat_low[i, :, :] = (lat_low - ds.variables['lat'][i]) ** 2

min_lon = np.zeros((ds.variables['lon'].shape[0], 1200, 1200))
min_lon_low = np.zeros((ds.variables['lon'].shape[0], lon_low.shape[0], lon_low.shape[1]))
for i in xrange(ds.variables['lon'].shape[0]):
    min_lon[i, :, :] = (ds1['lon'] - ds.variables['lon'][i]) ** 2
    min_lon_low[i, :, :] = (lon_low - ds.variables['lon'][i]) ** 2

for day in days:
    air_new = np.zeros((1200, 1200))
    air_new_low = np.zeros((xnew.shape[0], ynew.shape[0]))

    ind_lat = []
    ind_lon = []
    print 'ds.variables[air].shape[0]:', ds.variables['air'].shape[0]
    for i in xrange(0, 1200):
        for j in xrange(0, 1200):
            ind_lat0 = np.argmin(min_lat[:, i, j])
            ind_lon0 = np.argmin(min_lon[:, i, j])
            air_new[i, j] = ds.variables['air'][day, 0, ind_lat0, ind_lon0]


    ind_lat = []
    ind_lon = []
    print 'ds.variables[air].shape[0]:', ds.variables['air'].shape[0]
    for i in xrange(0, xnew.shape[0]):
        for j in xrange(0, ynew.shape[0]):
            ind_lat0 = np.argmin(min_lat_low[:, i, j])
            ind_lon0 = np.argmin(min_lon_low[:, i, j])
            air_new_low[i, j] = ds.variables['air'][day, 0, ind_lat0, ind_lon0]

    f = interpolate.interp2d(xnew, ynew, air_new_low, kind='cubic')
    air_new_interp = f(x, y)

    out_dir = '/group_workspaces/jasmin2/baci/sigil/air_temperature/1_km/%s/' % options.tile
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    np.savez(out_dir + 'air_new_%s_%d_d%d' % (options.tile, year, day),\
             air_new=air_new, air_new_interp=air_new_interp, ind_lat=ind_lat, ind_lon=ind_lon)