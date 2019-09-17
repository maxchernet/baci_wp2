import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import gdal
import os
import optparse

dir_main = '/group_workspaces/jasmin2/baci/sigil/baci_wp2_files/13_europe/'

parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(),
                                   usage=globals()['__doc__'])
parser.add_option('--file_name', action="store", dest="file_name",
                      type=str, help="File name")
parser.add_option('--var', action="store", dest="var",
                      type=str, help="Variable")
parser.add_option('--band', action="store", dest="band",
                      type=str, help="Band")

(options, args) = parser.parse_args()

file_name = options.file_name
var = options.var
band = options.band

print 'File name:'
print file_name

#ds = nc.Dataset(dir_main + 'lst/lst_h17v04_2015_7day.nc')
ds = nc.Dataset(file_name)
print ds

x_size = 20
y_size = 20
size = 1200

mean_arr = np.zeros((ds[var].shape[0], size/x_size, size/y_size))
std_arr =  np.zeros((ds[var].shape[0], size/x_size, size/y_size))

drv_tif = gdal.GetDriverByName('GTiff')

out_name_mean = dir_main + 'tiff/' + file_name.split('/')[-1].split('.')[0] + '_' + band + '.tif'
gdal_out_mean = drv_tif.Create(out_name_mean, mean_arr.shape[1], mean_arr.shape[2], mean_arr.shape[0], gdal.GDT_Float32)

out_name_sd = dir_main + 'tiff/' + file_name.split('/')[-1].split('.')[0] + '_' + band + '_sd.tif'
gdal_out_sd = drv_tif.Create(out_name_sd, mean_arr.shape[1], mean_arr.shape[2], mean_arr.shape[0], gdal.GDT_Float32)

try:
        gdal_out_mean.SetGeoTransform(ds['crs'].GeoTransform)
        gdal_out_mean.SetProjection(str(ds['crs'].spatial_ref))
        gdal_out_sd.SetGeoTransform(ds['crs'].GeoTransform)
        gdal_out_sd.SetProjection(str(ds['crs'].spatial_ref))

except:
        print 'crs does not exist'

for i0, i in enumerate(xrange(0, size, x_size)):
        for j0, j in enumerate(xrange(0, size, y_size)):
                mean_arr[:, i0, j0] = np.mean(ds[var][:, i:i+x_size, j:j+y_size], axis=(1,2))
		std_arr[:, i0, j0] = np.std(ds[var][:, i:i+x_size, j:j+y_size], axis=(1,2))

#print mean_arr.shape, std_arr.shape
mean_arr = np.ma.MaskedArray(mean_arr, mask=mean_arr==0)
std_arr = np.ma.MaskedArray(std_arr, mask=std_arr==0)

#print gdal_out.GetGeoTransform()
#print gdal_out.GetProjection()

for i in xrange(mean_arr.shape[0]):
        gdal_out_mean.GetRasterBand(i+1).WriteArray(mean_arr[i, :, :])
        gdal_out_sd.GetRasterBand(i+1).WriteArray(std_arr[i, :, :])

gdal_out_mean = None
gdal_out_sd = None

#os.system('gdalinfo %s' % out_name)

#plt.subplot(121)
#plt.imshow(mean_arr[0,:,:])

#plt.subplot(122)
#plt.imshow(std_arr[0,:,:])

#plt.show()
