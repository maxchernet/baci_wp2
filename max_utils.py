"""
Some useful utils which are used by different codes
"""

import netCDF4 as nc
import osr
import gdal
import os
import numpy as np
import glob
#import modis_data


def make_vrt(tile, year, loc, dir_vrt):
    # Create virtual datasets based of MODIS images

    m = modis_data.MODISFinder(tile=tile, year=year, head_loc=loc)
    m.create_virtual_data(dir_vrt)


# *************************************************************************

def locate(pattern, root=os.getcwd()):
    """
    Recursively locate a file matching a pattern (Jose).
    """

    import fnmatch
    for path, dirs, files in os.walk(root):
        for filename in [os.path.abspath(os.path.join(path, filename)) \
                         for filename in files if fnmatch.fnmatch(filename, pattern)]:
            # Check file size. Max
            if os.path.getsize(filename) > 100000:
                # yield filename, int(filename.split("/")[-1].split(".")[1][1:5]), \
                #       int(filename.split("/")[-1].split(".")[1][5:])

                yield filename



def get_n_days(year):
    n_days = (datetime.date(year, 12, 31) - datetime.date(year, 1, 1)).days + 1
    date_range = np.array([(datetime.date(year, 1, 1) + datetime.timedelta(days=i)).strftime('%Y.%m.%d') \
                           for i in range(0, n_days)])
    return n_days, date_range



def get_modis_coord(lat, lon):
    """
    Transfrom coords from Sentinel to MODIS
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
    tx = osr.CoordinateTransformation(wgs84, modis_sinu)
    # transform couple of coordinates to MODIS
    mx, my, mz = tx.TransformPoint(lon, lat)

    return mx, my




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
    lon, lat, mz = tx.TransformPoint(mx, my)

    return lat, lon





def create_netcdf(netcdf_file, cols, rows, n_days, do_zlib=True):
    """
    Create an netCDF file where we store the results

    Parameters
    ----------
    netcdf_file: str
        file to write

    Returns
    --------
    rootgrp: netCDF4
        a pointer to a netcdf file
    """

    # least_significant_digit for createVariable(...)
    lsd = None

    # create a new netCDF file
    rootgrp = nc.Dataset(netcdf_file, "w")
    rootgrp.description = "NBAR reflectance and BB albedo around (9x9) the Hainich fluxnet site"

    # Create groups for NBAR reflectance and BB albedo
    refl_grp = rootgrp.createGroup("reflectance")
    refl_grp.description = "Normalized at nadir (NBAR) reflectance"
    albedo_grp = rootgrp.createGroup("albedo")
    albedo_grp.description = "Broad band (BB) albedo"

    # Create dimensions for reflectance data. I.e. time series of reflectance
    # are determined by x, y and a day since 2000
    x = rootgrp.createDimension("x", cols)
    y = rootgrp.createDimension("y", rows)
    day = rootgrp.createDimension("day", n_days)
    band = rootgrp.createDimension("band", 7)
    str_dim = rootgrp.createDimension("str_dim", 10)

    # We can set zlib=True for compression as an argument
    # of the createVariable function
    date_str = rootgrp.createVariable("date_str", "S1", ("day", "str_dim"), zlib=do_zlib, least_significant_digit=lsd)
    date_str.units = "string representation of date: yyyy.mm.dd"
    # date = rootgrp.createVariable("julday", "i")
    date = rootgrp.createVariable("julday", "i", "day", zlib=do_zlib, least_significant_digit=lsd)
    date.units = "Julian day"

    # Create variables of NBAR reflectance: 7 MODIS bands
    for i in range(1, 8):
        refl = refl_grp.createVariable("refl_b%d" % i, "f4", ("x", "y", "day"), zlib=do_zlib, \
                                       least_significant_digit=lsd)
        refl.units = "surface bidirectional reflectance, band %d" % i

        # reflectance uncertainty
        refl_sd = refl_grp.createVariable("refl_b%d_sd" % i, "f4", ("x", "y", "day"), zlib=do_zlib, \
                                          least_significant_digit=lsd)
        refl_sd.units = "uncertainty of surface bidirectional reflectance, band %d" % i

    y_fwd = refl_grp.createVariable("y_fwd", "f4", ("x", "y", "day", "band"), zlib=do_zlib, \
                                    least_significant_digit=lsd)

    albedo_vis = albedo_grp.createVariable("albedo_vis", "f4", ("x", "y", "day"), zlib=do_zlib, \
                                           least_significant_digit=lsd)
    albedo_vis.units = "broad band albedo"
    albedo_nir = albedo_grp.createVariable("albedo_nir", "f4", ("x", "y", "day"), zlib=do_zlib, \
                                           least_significant_digit=lsd)
    albedo_nir.units = "broad band albedo"
    albedo_swir = albedo_grp.createVariable("albedo_swir", "f4", ("x", "y", "day"), zlib=do_zlib, \
                                            least_significant_digit=lsd)
    albedo_swir.units = "broad band albedo"

    # albedo uncertainty
    albedo_vis_sd = albedo_grp.createVariable("albedo_vis_sd", "f4", ("x", "y", "day"), zlib=do_zlib, \
                                              least_significant_digit=lsd)
    albedo_vis_sd.units = "albedo standard deviation"
    albedo_nir_sd = albedo_grp.createVariable("albedo_nir_sd", "f4", ("x", "y", "day"), zlib=do_zlib, \
                                              least_significant_digit=lsd)
    albedo_nir_sd.units = "albedo standard deviation"
    albedo_swir_sd = albedo_grp.createVariable("albedo_swir_sd", "f4", ("x", "y", "day"), zlib=do_zlib, \
                                               least_significant_digit=lsd)
    albedo_swir_sd.units = "albedo standard deviation"

    # latitude and longitude arrays
    lat = rootgrp.createVariable('lat', "f4", ("x", "y"), zlib=do_zlib)
    lat.units = "latitude"
    lon = rootgrp.createVariable('lon', "f4", ("x", "y"), zlib=do_zlib)
    lon.units = "longitude"

    return rootgrp
    
    
    
def get_corner(fs, cols, rows, out_dir, product_name, opt_py_file, pattern='', root=''):
    """
    Get upper left and lower right corners of image which was divided on blocks
    """
    
    if fs.site['ref_file'] == '':
        # if we don't have a reference for this site use as
        # a reference first file in list
        f_list = np.array([file for file in glob.glob(fs.site['loc_sentinel-1_modis'] + pattern)])
    else:
        f_list = np.array([file for file in glob.glob(fs.site['loc_sentinel-1_modis'] + fs.site['ref_file'])])
        
    if f_list != []:
        # Open a reference image
        ref_g = gdal.Open(f_list[0])
        geo = ref_g.GetGeoTransform()
        #proj = ref_g.GetProjection()
        ref_data = ref_g.ReadAsArray()
        print f_list[0].split('/')[-1], ':',  ref_data.shape

        # Get upper left and lower right corners in degrees or in map coord.
        ulx = 0 * geo[1] + geo[0]
        lrx = ref_data.shape[2] * geo[1] + geo[0]
        uly = 0 * geo[5] + geo[3]
        lry = ref_data.shape[1] * geo[5] + geo[3]
    else:
        print 'Reference file %s for this area does not exists' % sent1_file
        print 'use default size of area'
        ulx, uly = get_modis_coord(fs.site['lat'] + 0.6, fs.site['lon'] - 0.6)
        lrx, lry = get_modis_coord(fs.site['lat'] - 0.6, fs.site['lon'] + 0.6)
        
    return ulx, uly, lrx, lry, geo[1], geo[5]
    
    
def make_vrt(tile, year, loc, dir_vrt):

        # Create virtual datasets based of MODIS images

        m = modis_data.MODISFinder( tile=tile, year=year, head_loc=loc)
        m.create_virtual_data(dir_vrt)
