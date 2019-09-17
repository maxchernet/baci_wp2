import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import h5py
import julday
import os
import optparse
import gdal, osr





# ***********************************************************************************




def get_latlon(mx, my):
    """
    Transfrom coords from meters of MODIS sinusoidal to lat/lon WGS84.

    Parameters
    ----------
    mx: float
    my: float

    Returns
    -------
    lat: float
    lon: float
    """

    # Get spatial reference
    wgs84 = osr.SpatialReference()
    # and initialize it by EPSG
    wgs84.ImportFromEPSG(4326)
    modis_sinu = osr.SpatialReference()
    # initialize by a string
    modis_sinu.ImportFromProj4("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
    # define coordinate transformation
    # in this case from WGS84 to MODIS sinusoidal
    tx = osr.CoordinateTransformation(modis_sinu, wgs84)
    # transform couple of coordinates to MODIS
    # print 'mx %.6f, my %.6f:' % (mx, my)
    # print type(mx), type(my)
    # print gdal.VersionInfo()
    lon, lat, mz = tx.TransformPoint(float(mx), float(my))

    return lat, lon



# **********************************************************************************************************************



def do_average(file_in, file_out, year, step=30, ref_file='', meter_file=''):
    """
    get an array from a geo image.

    Parameters
    -----------
    f_name : str
        file name
    px : int
        x pixel
    py : int
        y pixel
    cols : int
        number of columns
    rows : int
        number of rows
    sub_ds_n : int
        number of a subdataset
    """

    print 'Read input file:'
    ds_in = h5py.File(file_in, 'r')

    # ds_out = nc.Dataset(path_out + 'lst_average_%s_%d_%dday.nc' % (tile, year, step), 'w')
    ds_out = nc.Dataset(file_out, 'w')
    ds_out.createDimension('x', 1200)
    ds_out.createDimension('y', 1200)
    # date_step = np.arange(1, ds_in['doy'].shape[0] + 1, step)
    date_step = np.arange(1, 366, step)
    ds_out.createDimension('time', date_step.shape[0])

    print 'Create variables:'
    lstVar = ds_out.createVariable('lst', 'f8', ('time', 'x', 'y'), zlib=True)
    lstVar.standard_name = 'LST'
    lstVar.long_name = 'Land Surface Temperature'
    lstVar.units = 'K'

    ds_out.createVariable('lst_sd', 'f8', ('time', 'x', 'y'), zlib=True)
    lstVar.standard_name = 'LST SD'
    lstVar.long_name = 'Land Surface Temperature Uncertainty'
    lstVar.units = 'K'

    timeVar = ds_out.createVariable('time', 'f8', 'time')
    timeVar.units = 'days since -4712-01-01 00:00:00'


    xVar = ds_out.createVariable('x', 'f8', ('x'), zlib=True)
    xVar.standard_name = 'projection_x_coordinate'
    xVar.long_name = 'x distance on the projection plane from the origin'
    xVar.units = 'm'
    yVar = ds_out.createVariable('y', 'f8', ('y'), zlib=True)
    yVar.standard_name = 'projection_y_coordinate'
    yVar.long_name = 'y distance on the projection plane from the origin'
    yVar.units = 'm'

    meters = np.load(meter_file)
    xVar[:] = meters['mx']
    yVar[:] = meters['my']

    lonVar = ds_out.createVariable('lon', 'f8', ('y', 'x'), zlib=True)
    latVar = ds_out.createVariable('lat', 'f8', ('y', 'x'), zlib=True)

    lonVar.long_name = "longitude"
    lonVar.units = 'degrees_east'

    latVar.long_name = "latitude"
    latVar.units = 'degrees_north'


    # lat = np.zeros((yVar.shape[0], xVar.shape[0]))
    # lon = np.zeros((yVar.shape[0], xVar.shape[0]))
    # get lat-lon by projected (map) coordinates
    print 'Get Lat/Lon:'
    for i in xrange(xVar.shape[0]):
        for j in xrange(yVar.shape[0]):
            # print 'xVar[i], yVar[j]', xVar[i], yVar[j]
            latVar[j, i], lonVar[j, i] = get_latlon(xVar[i], yVar[j])

    # lonVar[:] = dsin['lon'][:]
    # latVar[:] = dsin['lat'][:]
    print 'Create crs:'
    crs = ds_out.createVariable('crs', 'S1')
    crs.grid_mapping_name = 'sinusoidal'
    crs.longitude_of_central_meridian = '0.0'
    crs.longitude_of_prime_meridian = '0.0'
    crs.semi_major_axis = '6371007.181'
    crs.false_easting = '0.0'
    crs.false_northing = '0.0'
    crs.earth_radius = '6371.229'
    gds = gdal.Open(ref_file)
    geo_ref = gds.GetGeoTransform()
    crs.spatial_ref = gds.GetProjection()
    crs.GeoTransform = gds.GetGeoTransform()

    print 'Do average:'
    for i, d in enumerate(date_step):
        ind = np.logical_and(ds_in['doy'][:] >= d, ds_in['doy'][:] < (d + step))
        if np.max(ind) == True:
            img = ds_in['lst'][ind, :, :] * 0.02
            img[img == 0] = np.NaN
            ds_out['lst'][i, :, :] = np.nanmean(img, axis=0)
            ds_out['lst_sd'][i, :, :] = np.nanstd(img, axis=0)
        else:
            ds_out['lst'][i, :, :] = np.NaN
            ds_out['lst_sd'][i, :, :] = np.NaN

        ds_out['lst'][i,:,:][np.isnan(ds_out['lst'][i,:,:])] = ma.masked

        ds_out['time'][i] = julday.date2jul(julday.doy2date(year, d))

    return 0




# *******************************************************************************************




if __name__ == "__main__":

    parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(),
                                   usage=globals()['__doc__'])

    # parser.add_option('--path_in', action="store", dest="path_in",
    #                   type=str, help="Input path")

    # parser.add_option('--path_out', action="store", dest="path_out",
    #                   type=str, help="Output paht")

    parser.add_option('--file_in', action="store", dest="file_in",
                      type=str, help="Input file")

    parser.add_option('--file_out', action="store", dest="file_out",
                      type=str, help="Output file")

    parser.add_option('--year', action="store", dest="year",
                      type=str, help="year")

    parser.add_option('--tile', action="store", dest="tile",
                      type=str, help="MODIS tile")

    parser.add_option('--step', action="store", dest="step",
                      type=str, help="Step in number of days")

    parser.add_option('--ref_file', action="store", dest="ref_file",
                      type=str, help="MODIS reference hdf file")

    parser.add_option('--meter_file', action="store", dest="meter_file",
                      type=str, help="Meter coordinates")

    (options, args) = parser.parse_args()

    # path_out = options.path_out,\
    # path_in = options.path_in,\
    # tile = options.tile,\

    do_average(file_in = options.file_in, \
               file_out = options.file_out, \
               year = int(options.year), \
               step = int(options.step),\
               ref_file = options.ref_file,\
               meter_file = options.meter_file)
