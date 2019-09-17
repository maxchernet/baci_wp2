"""
Combine MODIS normalized time series from a set of netCDF files. I.e. do
a mosaic - make one big image from many small blocks, where each block is a netCDF file.
We have MODIS optical, MODIS LST and Sentinel-1 backscater in VV and HH polarization.
"""

from baci_fs import *
import os
import sys
import netCDF4 as nc
import gdal
import osr
import numpy as np
import glob
import julday
import datetime as dt
import socket
import optparse
import pdb

import max_utils
# import baci_netcdf as nc_baci


def filter_window(img):
    """
    filter artefacts
    :param img: array
    :return: img
    """
    for i in xrange(1, img.shape[0]-1):
        for j in xrange(1, img.shape[0]-1):
            tmp = np.zeros((3, 3))
            tmp[:, :] = img[i - 1:i + 2, j - 1:j + 2]
            tmp[1, 1] = 0
            if np.logical_and(np.mean(tmp) == 0., img[i, j] != 0):
                # print tmp
                img[i, j] = 0.
    return img



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





def create_nc_sent1(ncd_file, M, N, ndays, ndays_orig):
    """
    Create a netCDF file for Sentinel-1 backscatter
    """
    
    nc_big = nc.Dataset(ncd_file, 'w')

    nc_big.description = 'Smoothed Sentinel-1 backscatter with uncertainties'

    nchar_proj = 1000
    nchar_geo = 1000
    # nc_big.createDimension('nchar_proj', nchar_proj)
    # nc_big.createDimension('nchar_geo', nchar_geo)

    nc_big.createDimension('x', N)
    nc_big.createDimension('y', M)
    nc_big.createDimension('date', ndays)
    
    nc_big.createVariable('bs', 'f4', ('date', 'y', 'x'), zlib=True)
    nc_big.createVariable('bs_sd', 'f4', ('date', 'y', 'x'), zlib=True)
    nc_big.createVariable('bs_orig', 'f4', ('date', 'y', 'x'), zlib=True)
    nc_big.createVariable('bs_orig_sd', 'f4', ('date', 'y', 'x'), zlib=True)
    nc_big.createVariable('julday', 'i', 'date')

    # nc_big.createVariable('proj', 'S1', ('nchar_proj'))
    # nc_big.createVariable('geo_transform', 'S1', 'nchar_geo')

    # latitude and longitude arrays
    # lat = nc_big.createVariable('lat', "f4", ("y", "x"), zlib=True)
    # lat.units = "latitude"
    # lon = nc_big.createVariable('lon', "f4", ("y", "x"), zlib=True)
    # lon.units = "longitude"

    # crs = nc_big.createVariable('crs', 'S1')
    # crs.grid_mapping_name = 'latitude_longitude'
    # crs.longitude_of_prime_meridian = '0.0'
    # crs.semi_major_axis = ''
    # crs.inverse_flattening = ''
    #
    # print nc_big.variables['crs']
    #
    # nc_big.close()
    
    return nc_big, 0



def copy_ds_mod09(nc_big, ds, x1, x2, y1, y2, year):

    print 'copy ds mod09'
    print 'nc_big:', nc_big
    print 'ds:', ds

    for band in xrange(7):

        # try:

        # nc_big[0].groups['reflectance'].variables['refl_b%d' % (band + 1)][x1:x2, y1:y2, :] = \
        #     ds.groups['reflectance'].variables['refl_b%d' % (band + 1)][:]
        try:
            img = ds.groups['reflectance'].variables['refl_b%d' % (band + 1)][:]
            img[img<0] = 0
            img[img>1] = 0
            print 'nc_big[0].groups[reflectance]:', \
                nc_big[0].groups['reflectance'].variables['refl_b%d' % (band + 1)][:, x1:x2, y1:y2].shape
            print 'np.ma.masked_values(img, 0):', np.ma.masked_values(img, 0).shape

            nc_big[0].groups['reflectance'].variables['refl_b%d' %\
                                                      (band + 1)][:, x1:x2, y1:y2] = np.ma.masked_values(img, 0)

            nc_big[0].groups['reflectance'].variables['refl_b%d_sd' % (band + 1)][:, x1:x2, y1:y2] = \
                ds.groups['reflectance'].variables['refl_b%d_sd' % (band + 1)][:]

        except:
            print 'reflectance does not work'
            return nc_big

    try:
        print 'nc_big[1].variables', nc_big[1].variables

        nc_big[1].variables['y_fwd'][ :, :, x1:x2, y1:y2] = \
            ds.groups['reflectance'].variables['y_fwd'][:]

        nc_big[1].variables['y_orig'][:, :, x1:x2, y1:y2] = \
            ds.groups['reflectance'].variables['y_orig'][:]

        img = ds.groups['albedo'].variables['albedo_vis'][:]
        img[img < 0] = 0
        img[img > 1] = 0
        nc_big[0].groups['albedo'].variables['albedo_vis'][:, x1:x2, y1:y2] = np.ma.masked_values(img, 0)
            # ds.groups['albedo'].variables['albedo_vis'][:]

        img = ds.groups['albedo'].variables['albedo_nir'][:]
        img[img < 0] = 0
        img[img > 1] = 0
        nc_big[0].groups['albedo'].variables['albedo_nir'][:, x1:x2, y1:y2] = np.ma.masked_values(img, 0)
        # nc_big[0].groups['albedo'].variables['albedo_nir'][x1:x2, y1:y2, :] =\
        #     ds.groups['albedo'].variables['albedo_nir'][:]

        img = ds.groups['albedo'].variables['albedo_swir'][:]
        img[img < 0] = 0
        img[img > 1] = 0
        nc_big[0].groups['albedo'].variables['albedo_swir'][:, x1:x2, y1:y2] = np.ma.masked_values(img, 0)
        # nc_big[0].groups['albedo'].variables['albedo_swir'][x1:x2, y1:y2, :] =\
        #     ds.groups['albedo'].variables['albedo_swir'][:]

        nc_big[0].groups['albedo'].variables['albedo_vis_sd'][:, x1:x2, y1:y2] =\
            ds.groups['albedo'].variables['albedo_vis_sd'][:]

        nc_big[0].groups['albedo'].variables['albedo_nir_sd'][:, x1:x2, y1:y2] =\
            ds.groups['albedo'].variables['albedo_nir_sd'][:]

        nc_big[0].groups['albedo'].variables['albedo_swir_sd'][:, x1:x2, y1:y2] =\
            ds.groups['albedo'].variables['albedo_swir_sd'][:]

    except:
        print 'albedo does not work'
        return nc_big

    # try:
        # nc_big[0].variables['lat'][x1:x2, y1:y2] = ds.variables['lat'][:].T
        # nc_big[0].variables['lon'][x1:x2, y1:y2] = ds.variables['lon'][:].T

    # for some blocks 'julday' is empty
    if np.ma.is_masked(ds.variables['julday'][:]) == False:
        nc_big[0].variables['julday'][:] = ds.variables['julday'][:]

        # nc_big[0].variables['date_str'][:] = ds.variables['date_str'][:]

    # except:
    #     print 'lat/lon, julday or date does not work'
    #     return nc_big

    return nc_big




def get_ndays(tile, hdf_dir, ncd_file_pattern, year, file_mosaic):

    mosaic = np.loadtxt(file_mosaic, skiprows=1)
    print 'mosaic.shape', mosaic.shape
    isfile = False

    for i in xrange(mosaic.shape[0]):
        N = mosaic[i, 0]
        M = mosaic[i, 1]
        ulx = mosaic[i, 2]
        uly = mosaic[i, 3]
        #cols = mosaic[0, 4]
        #rows = mosaic[0, 5]

        nc_file = hdf_dir + ncd_file_pattern % (tile, year, ulx, uly)
        isfile =os.path.isfile(nc_file)
        if isfile:
            ds = nc.Dataset(nc_file, 'r')
            break
        else:
            print 'file doesnt exist:', nc_file

    if isfile == False:
        print 'no files exist'
        return -1, -1, -1, -1, -1

    try:
        ndays = ds.variables['date_jul'].shape[0]
    except:
        try:
            ndays = ds.variables['julday'].shape[0]
        except:
            # ndays = 365
            ndays = (dt.date(year, 12, 31) - dt.date(year, 1, 1)).days + 1

    try:
        ndays_orig = ds['date_jul_obs'][:].shape[0]
    except:
        ndays_orig = (dt.date(year, 12, 31) - dt.date(year, 1, 1)).days + 1

    return int(ulx), int(uly), int(N), int(M), ndays, ndays_orig





def copy_ds_sent1(nc_big, ds, x1, x2, y1, y2, year):
    """
    copy data from small datasets to one big dataset backscatter
    
    Parameters
    __________
    
    nc_big: netCDF
      a netCDF object - one big dataset backscatter
    ds: netCDF
      small netCDF dataset. I.e. a block
    x1: int
      upper x limit
    x2: int
      lower x limit
    y1: int
      upper y limit
    y2: int
      lower y limit
    """
    
    x3 = nc_big[0].variables['bs'][:, y1:y2, x1:x2].shape[2]
    y3 = nc_big[0].variables['bs'][:, y1:y2, x1:x2].shape[1]

    data = ds.variables['bs'][:, 0:y3, 0:x3].copy()
    data[np.isnan(data)] = 0

    # pdb.set_trace()

    nc_big[0].variables['bs'][:, y1:y2, x1:x2] = np.ma.masked_values(data, 0)

    data = ds.variables['bs_sd'][:, 0:y3, 0:x3].copy()
    data[np.isnan(data)] = 0
    nc_big[0].variables['bs_sd'][:, y1:y2, x1:x2] = np.ma.masked_values(data, 0)

    if nc_big[1] != 0:
        data = ds.variables['bs_orig'][:, 0:y3, 0:x3].copy()
        data[np.isnan(data)] = 0
        nc_big[1].variables['bs_orig'][:, y1:y2, x1:x2] = np.ma.masked_values(data, 0)

        data = ds.variables['bs_orig_sd'][:, 0:y3, 0:x3].copy()
        data[np.isnan(data)] = 0
        nc_big[1].variables['bs_orig_sd'][:, y1:y2, x1:x2] = np.ma.masked_values(data, 0)

    nc_big[0].variables['julday'][:] = ds.variables['date_jul'][:]
    
    
    
# *********************************************************************************************************

    
def copy_ds_mod11(nc_big, ds, x1, x2, y1, y2, year):
    """
    copy data from small datasets to one big netCDF dataset with MOD11 LST
    """

    # ds.variables['lst'][:][ds.variables['lst'][:]==0] = np.ma.masked
    # ds.variables['lst_sd'][:][ds.variables['lst_sd'][:] == 0] = np.ma.masked

    x3 = nc_big[0].variables['lst'][:, y1:y2, x1:x2].shape[2]
    y3 = nc_big[0].variables['lst'][:, y1:y2, x1:x2].shape[1]

    print 'size block:', ds.variables['lst'][:, 0:y3, 0:x3].shape[0]
    print 'size big:', nc_big[0].variables['lst'][:, y1:y2, x1:x2].shape[0]
    # MaskedArrayFutureWarning: setting an item on a masked array which has a shared mask will not copy the mask and
    # also change the original mask array in the future.
    # Copy one of the blocks to output array
    nc_big[0].variables['lst'][:, y1:y2, x1:x2] = np.ma.masked_values(ds.variables['lst'][:, 0:y3, 0:x3], 0)
    # nc_big[0].variables['lst'][:][y1:y2, x1:x2, :].unshare_mask()
    # Mask zero-pixels (water, etc.)
    # ind = nc_big[0].variables['lst'][y1:y2, x1:x2, :] == 0
    # nc_big[0].variables['lst'][:][y1:y2, x1:x2, :][ind] = True #np.ma.masked

    # nc_big[0].variables['lst_sd'][:][y1:y2, x1:x2, :] = ds.variables['lst_sd'][0:y3, 0:x3, :]
    nc_big[0].variables['lst_sd'][:, y1:y2, x1:x2] = np.ma.masked_values(ds.variables['lst_sd'][:, 0:y3, 0:x3], 0)
    # nc_big[0].variables['lst_sd'][:][y1:y2, x1:x2, :].unshare_mask()
    # ind = nc_big[0].variables['lst_sd'][:][y1:y2, x1:x2, :] == 0
    # nc_big[0].variables['lst_sd'][:][y1:y2, x1:x2, :][ind] = np.ma.masked

    # if np.min(ds.variables['lst'][:]) == 0:
    #     pdb.set_trace()


    #print 'min/max:', np.min(nc_big.variables['lst'][0:y3, 0:x3, :]), np.max(ds.variables['lst'][0:y3, 0:x3, :])
    #print 'min/max:', np.min(nc_big.variables['lst'][y1:y2, x1:x2, :]), np.max(ds.variables['lst'][0:y3, 0:x3, :])

    #nc_big.variables['date_jul'][:] = ds.variables['date_jul'][:]
    jd1 = julday.date2jul(dt.date(year, 1, 1))
    jd = np.arange(jd1, jd1 + ds.variables['lst'].shape[0])
    nc_big[0].variables['julday'][:] = ds.variables['date_jul'][:] #jd

    #nc_big.variables['lat'][y1:y2, x1:x2] = ds.variables['lat'][:].T
    #nc_big.variables['lon'][y1:y2, x1:x2] = ds.variables['lon'][:].T

    #********
    # original dates are always different
    #julday.jul2doy()

    # pdb.set_trace()

    # mask = np.array([julday.jul2doy(d) for d in ds.variables['date_jul_obs'][:]]) - 1
    #print mask
    # print ds.variables['date_jul_obs'][:].shape
    #print ds.variables['lst_orig'][0:y3, 0:x3, :]
    #print mask.shape
    #print nc_big.variables['lst_orig'][y1:y2, x1:x2, :].shape

    orig_size =  ds.variables['lst_orig'][:, 0:y3, 0:x3].shape[0]

    # print 'nc_big[1].variables[lst_orig][:orig_size, y1:y2, x1:x2]: ', nc_big[1].variables['lst_orig'][:orig_size, y1:y2, x1:x2].shape
    # print 'ds.variables[lst_orig][:, 0:y3, 0:x3]: ', ds.variables['lst_orig'][:, 0:y3, 0:x3].shape

    nc_big[1].variables['lst_orig'][:orig_size, y1:y2, x1:x2] =\
        np.ma.masked_values(ds.variables['lst_orig'][:, 0:y3, 0:x3], 0)
    nc_big[1].variables['lst_orig_sd'][:orig_size, y1:y2, x1:x2] = \
        np.ma.masked_values(ds.variables['lst_orig_sd'][:, 0:y3, 0:x3], 0)

    nc_big[1].variables['julday_obs'][:orig_size] = ds.variables['date_jul_obs'][:]


    # str_len = ds.variables['proj'][:].shape[0]
    # nc_big.variables['proj'][:str_len] = ds.variables['proj'][:]
    # str_len = ds.variables['geo_transform'][:].shape[0]
    # nc_big.variables['geo_transform'][:str_len] = ds.variables['geo_transform'][:]

    




def create_nc_mod09(netcdf_file, cols, rows, n_days, ndays_orig, do_zlib=True):
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
    rootgrp.description = "NBAR reflectance and BB albedo"

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
        refl = refl_grp.createVariable("refl_b%d" % i, "f4", ("day", "x", "y"), zlib=do_zlib, \
                                       least_significant_digit=lsd)
        refl.units = "surface bidirectional reflectance, band %d" % i

        # reflectance uncertainty
        refl_sd = refl_grp.createVariable("refl_b%d_sd" % i, "f4", ("day", "x", "y"), zlib=do_zlib, \
                                          least_significant_digit=lsd)
        refl_sd.units = "uncertainty of surface bidirectional reflectance, band %d" % i



    albedo_vis = albedo_grp.createVariable("albedo_vis", "f4", ("day", "x", "y"), zlib=do_zlib, \
                                           least_significant_digit=lsd)
    albedo_vis.units = "broad band albedo"
    albedo_nir = albedo_grp.createVariable("albedo_nir", "f4", ("day", "x", "y"), zlib=do_zlib, \
                                           least_significant_digit=lsd)
    albedo_nir.units = "broad band albedo"
    albedo_swir = albedo_grp.createVariable("albedo_swir", "f4", ("day", "x", "y"), zlib=do_zlib, \
                                            least_significant_digit=lsd)
    albedo_swir.units = "broad band albedo"

    # albedo uncertainty
    albedo_vis_sd = albedo_grp.createVariable("albedo_vis_sd", "f4", ("day", "x", "y"), zlib=do_zlib, \
                                              least_significant_digit=lsd)
    albedo_vis_sd.units = "albedo standard deviation"
    albedo_nir_sd = albedo_grp.createVariable("albedo_nir_sd", "f4", ("day", "x", "y"), zlib=do_zlib, \
                                              least_significant_digit=lsd)
    albedo_nir_sd.units = "albedo standard deviation"
    albedo_swir_sd = albedo_grp.createVariable("albedo_swir_sd", "f4", ("day", "x", "y"), zlib=do_zlib, \
                                               least_significant_digit=lsd)
    albedo_swir_sd.units = "albedo standard deviation"

    # latitude and longitude arrays
    lat = rootgrp.createVariable('lat', "f4", ("x", "y"), zlib=do_zlib)
    lat.units = "latitude"
    lon = rootgrp.createVariable('lon', "f4", ("x", "y"), zlib=do_zlib)
    lon.units = "longitude"

    name_fwd = ''.join(netcdf_file.split('.')[:-1]) + '_fwd.nc'
    root_fwd = nc.Dataset(name_fwd, 'w')
    x = root_fwd.createDimension("x", cols)
    y = root_fwd.createDimension("y", rows)
    day = root_fwd.createDimension("day", n_days)
    band = root_fwd.createDimension("band", 7)
    y_fwd = root_fwd.createVariable("y_fwd", "f4", ("day", "band", "x", "y"), zlib=do_zlib, \
                                    least_significant_digit=lsd)
    y_orig = root_fwd.createVariable("y_orig", "f4", ("day", "band", "x", "y"), zlib=do_zlib, \
                                    least_significant_digit=lsd)

    return rootgrp, root_fwd



    
def create_nc_mod11(ncd_file, M, N, ndays, ndays_orig):
    """
    Create a netCDF file for MOD11 LST
    """

    nchar_proj = 1000
    nchar_geo = 1000
    
    root_inv = nc.Dataset(ncd_file, 'w')
    root_inv.description = 'Smoothed LST with uncertainties'
    
    root_inv.createDimension('x', N)
    root_inv.createDimension('y', M)
    root_inv.createDimension('time', ndays)

    # root_inv.createDimension('x_obs', N)
    # root_inv.createDimension('y_obs', M)
    # root_inv.createDimension('time_obs', ndays_orig)

    # root_inv.createDimension('x_obs', N)
    # root_inv.createDimension('y_obs', M)
    # root_inv.createDimension('date_obs', ndays)

    # root_inv.createDimension('nchar_proj', nchar_proj)
    # root_inv.createDimension('nchar_geo', nchar_geo)
    
    root_inv.createVariable('lst', 'f4', ('time', 'y', 'x'), zlib=True)
    root_inv.createVariable('lst_sd', 'f4', ('time', 'y', 'x'), zlib=True)

    # root_inv.createVariable('lst_orig', 'f4', ('time_obs', 'x_obs', 'y_obs'), zlib=True)
    # root_inv.createVariable('lst_orig_sd', 'f4', ('time_obs', 'x_obs', 'y_obs'), zlib=True)

    root_inv.createVariable('julday', 'i', 'time')
    # nc_out.createVariable('date_jul_obs', 'i', 'time_obs')
    # root_inv.createVariable('date_jul_obs', 'i', 'date_obs')

    # root_inv.createVariable('proj', 'S1', ('nchar_proj'))
    # root_inv.createVariable('geo_transform', 'S1', 'nchar_geo')

    # latitude and longitude arrays
    lat = root_inv.createVariable('lat', "f4", ("y", "x"), zlib=True)
    lat.units = "latitude"
    lon = root_inv.createVariable('lon', "f4", ("y", "x"), zlib=True)
    lon.units = "longitude"

    
    name_fwd = ''.join(ncd_file.split('.')[:-1]) + '_fwd.nc'
    print 'name_fwd: ', name_fwd
    root_fwd = nc.Dataset(name_fwd, 'w')

    root_fwd.description = 'Smoothed LST with uncertainties'

    # root_fwd.createDimension('x', N)
    # root_fwd.createDimension('y', M)
    # root_fwd.createDimension('date', ndays)

    root_fwd.createDimension('x_obs', N)
    root_fwd.createDimension('y_obs', M)
    root_fwd.createDimension('time_obs', ndays_orig)

    root_fwd.createVariable('julday_obs', 'i', 'time_obs')
    root_fwd.createVariable('lst_orig', 'f4', ('time_obs', 'y_obs', 'x_obs'), zlib=True)
    root_fwd.createVariable('lst_orig_sd', 'f4', ('time_obs', 'y_obs', 'x_obs'), zlib=True)
    
    return root_inv, root_fwd
    
    
    
def scale_mod11(fs, nc_1km, ndays, year):
    """
    Rescale MOD11 1km to MODIS 500m
    
    Parameters
    -----------
    nc_1km: netCDF object
      1km netCDF LST
    reference_500m: string
      name of a file of reference 500m size
    """
    
    # /group_workspaces/cems/baci/sigil/baci_wp2_files/sent1_viterbo_2015_vv_677x458.nc
    
    # We need the same size for all images in time series.
    # For this we use one reference image
    # All images must have the same projection and geo. units!
    if fs.site['ref_file'] == '':
        # if we don't have a reference for this site use as
        # a reference first file in list
        f_list = np.array([file for file in glob.glob(fs.site['loc_sentinel-1_modis'] + '*')])
    else:
        f_list = np.array([file for file in glob.glob(fs.site['loc_sentinel-1_modis'] + fs.site['ref_file'])])
        
    
    print nc_1km.variables['lst'].shape
    in_drv =  gdal.GetDriverByName('MEM')
    out_drv = gdal.GetDriverByName('MEM')
    
    # get size of 1km mod11 image
    xs1 = nc_1km.variables['lst'].shape[1]
    ys1 = nc_1km.variables['lst'].shape[0]
    
    # get size of a reference image
    # ds_ref = nc.Dataset(reference_500m)
    # xs2 = ds_ref.variables['bs'].shape[1]
    # ys2 = ds_ref.variables['bs'].shape[0]
    
    gds_500 = gdal.Open(f_list[0])
    xs2 = gds_500.RasterXSize
    ys2 = gds_500.RasterYSize
    
    in_wgs =  in_drv.Create('', xs1, ys1, 1, gdal.GDT_CFloat32)
    out_wgs = out_drv.Create('', xs2, ys2, 1, gdal.GDT_CFloat32)
    
    # geo_in = gds_500.GetGeoTransform()
    geo_out = gds_500.GetGeoTransform()
    proj_in = gds_500.GetProjection()
    proj_out = gds_500.GetProjection()
    # projections and geo transformations are the same except
    # pixel sizes
    # geo_in[1] = 926.6254331391667
    # geo_in[5] = -926.625433139167
    
    geo_in = (geo_out[0], 926.6254331391667, geo_out[2], geo_out[3], geo_out[4], -926.625433139167)
    
    nc_500 = create_nc_mod11('/group_workspaces/cems/baci/sigil/baci_wp2_files/mod11_500m_viterbo_%d.nc' % year, xs2, ys2, ndays)
    
    in_wgs.SetGeoTransform(geo_in)
    out_wgs.SetGeoTransform(geo_out)
    in_wgs.SetProjection(proj_in)
    out_wgs.SetProjection(proj_out)
    
    nc_500.variables['date_jul'] = nc_1km.variables['date_jul']
    dic = ['lst', 'lst_sd']
    for i in xrange(nc_1km.variables['lst'].shape[2]):
        for s in dic:
            in_wgs.GetRasterBand(1).WriteArray(nc_1km.variables[s][:,:,i])
            wgs84 = osr.SpatialReference()
            wgs84.ImportFromEPSG(4326)
            
            res = gdal.ReprojectImage(in_wgs, out_wgs,\
                                      wgs84.ExportToWkt(), wgs84.ExportToWkt(),
                                      gdal.GRA_Average)
                                      
            nc_500.variables[s][:,:,i] = out_wgs.GetRasterBand(1).ReadAsArray()
            # nc_500.variables['lst_sd'][:,:,i] = out_wgs.GetRasterBand(1).ReadAsArray()
    
    return nc_500



def scale_mod09(fs, nc_big, ndays):
    return 0



def scale_sent1(fs, nc_big, ndays):
    return 0


# where LST files are stored:
# /group_workspaces/cems/baci/sigil/mod11a1/h18v04/

def run_job(tile, year, file_mosaic, ncd_file_pattern, hdf_dir, out_file, latlon_file, create_func,\
            copy_func, scale_func, ref_file):

    os.system('head -n 1 %s' % file_mosaic)

    if os.path.isfile(file_mosaic):
        mosaic = np.loadtxt(file_mosaic, skiprows=1)
        print 'mosaic.shape', mosaic.shape
    else:
        return -1

    # N = mosaic[0,0]
    # M = mosaic[0,1]
    # ulx1 = mosaic[0,2]
    # uly1 = mosaic[0,3]
    # cols = mosaic[0,4]
    # rows = mosaic[0,5]

    ulx0, uly0, N, M, ndays, n_days_orig = get_ndays(tile, hdf_dir, ncd_file_pattern, year, file_mosaic)

    print 'number of days in %d = %d' % (year, ndays)
    if ndays == -1:
        return -1

    # latlon_file = latlon_file % (fs.site['name'], year, N, M)

    print 'N, M', N, M

    backs = np.zeros((M, N, ndays))
    # rho_fwd = np.zeros((M, N, 365, 7))
    lat = np.zeros((M, N))
    lon = np.zeros((M, N))

    if os.path.isfile(out_file):
        os.system('rm ' + out_file)

    func_str = {'create_nc_mod09': create_nc_mod09,\
                'create_nc_mod11': create_nc_mod11,\
                'create_nc_sent1': create_nc_sent1,\
                'copy_ds_mod09': copy_ds_mod09,\
                'copy_ds_mod11': copy_ds_mod11,\
                'copy_ds_sent1': copy_ds_sent1}

    # nc_big = create_func(*(out_file, M, N, ndays))
    nc_big = func_str[create_func](*(out_file, M, N, ndays, n_days_orig))

    # nc_big = create_nc_mod11(out_file, M, N, ndays)
    # nc_big = max_utils.create_netcdf(out_file, M, N, ndays)
    # print 'nc_big: ', nc_big

    img_rgb = np.zeros((M, N, 3))

    print 'img_rgb.shape', img_rgb.shape
    print 'mosaic.shape[0]:', mosaic.shape[0]

    f_err = open(hdf_dir + 'missed_files.txt', 'w')
    for i in xrange(mosaic.shape[0]):

        N = int(mosaic[i, 0])
        M = int(mosaic[i, 1])
        ulx = mosaic[i, 2]
        uly = mosaic[i, 3]
        cols = mosaic[i, 4]
        rows = mosaic[i, 5]

        # nc_file = hdf_dir + ncd_file_pattern % (fs.site['name'], year, pol, N, M, ulx, uly)
        nc_file = hdf_dir + ncd_file_pattern % (tile, year, ulx, uly)
        if os.path.isfile(nc_file):
            print 'opening file N %d: %s' % (i, nc_file)
            try:
                ds = nc.Dataset(nc_file, 'r')
            except:
                print 'the file is corrupted'

            # y1 = int(uly) - uly0
            # y2 = int(uly + cols) - uly0
            # x1 = int(ulx) - ulx0
            # x2 = int(ulx + rows) - ulx0

            y1 = int(uly)
            y2 = int(uly + cols)
            x1 = int(ulx)
            x2 = int(ulx + rows)


            print 'x1=%d, y1=%d, x2=%d, y2=%d, N=%d, M=%d, col=%d, rows=%d' % (x1, y1, x2, y2, N, M, cols, rows)
            if x2 < (N + cols) and y2 < (M + rows):
                # try:
                print 'copying data'
                func_str[copy_func](*(nc_big, ds, x1, x2, y1, y2, year))
                # pdb.set_trace()

                # var = nc_big[0].groups['reflectance'].variables['refl_b1'][:]
                # print 'nc_big all min max:', np.min(var), np.max(var)

                # except:
                #    print 'something wrong with the file...'
                # copy_ds_mod11(nc_big, ds, x1, x2, y1, y2)
                # copy_ds_mod09(nc_big, ds, x1, x2, y1, y2)
            elif np.logical_and(N < 1200, M < 1200):
                print 'copying data'
                func_str[copy_func](*(nc_big, ds, x1, x2, y1, y2, year))
            else:
                print 'data not copied: x2=%d, y2=%d, N=%d, M=%d' % (x2, y2, N, M)

        else:
            print 'file %s does not exist' % nc_file
            f_err.write('%s\n' % nc_file)
    f_err.close()

    add_geo = True
    if add_geo:
        xVar = nc_big[0].createVariable('x', 'f8', ('x'), zlib=True)
        # xVar.long_name = "x-coordinate in Cartesian system"
        xVar.standard_name = 'projection_x_coordinate'
        xVar.long_name = 'x distance on the projection plane from the origin'
        xVar.units = 'm'
        yVar = nc_big[0].createVariable('y', 'f8', ('y'), zlib=True)
        # yVar.long_name = "y-coordinate in Cartesian system"
        yVar.standard_name = 'projection_y_coordinate'
        yVar.long_name = 'y distance on the projection plane from the origin'
        yVar.units = 'm'

        meters = np.load('latlon/%s_2015_meters.npz' % tile)
        xVar[:] = meters['mx']
        yVar[:] = meters['my']

        lonVar = nc_big[0].createVariable('lon', 'f8', ('y', 'x'), zlib=True)
        latVar = nc_big[0].createVariable('lat', 'f8', ('y', 'x'), zlib=True)

        lonVar.long_name = "longitude"
        lonVar.units = 'degrees_east'

        latVar.long_name = "latitude"
        latVar.units = 'degrees_north'

        # get lat-lon by projected (map) coordinates
        for i in xrange(N):
            for j in xrange(M):
                latVar[j, i], lonVar[j, i] = get_latlon(xVar[i], yVar[j])

        # nc_big[0].variables['lat'][:] = lat
        # nc_big[0].variables['lon'][:] = lon



        # try:

        crs = nc_big[0].createVariable('crs', 'S1')
        crs.grid_mapping_name = 'sinusoidal'
        crs.longitude_of_central_meridian = '0.0'
        crs.longitude_of_prime_meridian = '0.0'
        crs.semi_major_axis = '6371007.181'
        crs.false_easting = '0.0'
        crs.false_northing = '0.0'
        crs.earth_radius = '6371.229'
        # crs.inverse_flattening = ''
        # crs.spatial_ref = proj
        # crs.GeoTransform = geo

        # nc_big = nc.Dataset(out_file)
        # crs = nc_big.variables['crs']
        # crs.grid_mapping_name = 'latitude_longitude'
        # crs.longitude_of_prime_meridian = '0.0'
        # crs.semi_major_axis = ''
        # crs.inverse_flattening = ''

        print 'lat/lon reference file: ' + ref_file
        gds = gdal.Open(ref_file)
        geo_ref = gds.GetGeoTransform()
        crs.spatial_ref = gds.GetProjection()
        crs.GeoTransform = gds.GetGeoTransform()

        # define arrays for lat-lon information

        # lat = np.zeros((M, N))
        # lon = np.zeros((M, N))
        # mx = np.zeros(N)
        # my = np.zeros(M)

        # make vectors with x and y coordinates of image
        # mx = np.arange(0, N) * geo_ref[1] + geo_ref[0]
        # my = np.arange(0, M) * geo_ref[5] + geo_ref[3]




        # except:
        #     print 'error in coping geo/proj strings'


    
    


if __name__ == '__main__':

    parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(),
                                   usage=globals()['__doc__'])

    # file_mosaic, ncd_file_pattern, hdf_dir, out_file, latlon_file, create_func, copy_func, scale_func
    parser.add_option('--tile', action="store", dest="tile",
                      type=str, help="tile")

    parser.add_option('--year', action="store", dest="year",
                      type=str, help="year")

    parser.add_option('--file_mosaic', action="store", dest="file_mosaic",
                      type=str, help="file_mosaic")

    parser.add_option('--ncd_file_pattern', action="store", dest="ncd_file_pattern",
                      type=str, help="ncd_file_pattern")

    parser.add_option('--hdf_dir', action="store", dest="hdf_dir",
                      type=str, help="hdf_dir")

    parser.add_option('--out_file', action="store", dest="out_file",
                      type=str, help="out_file")

    parser.add_option('--latlon_file', action="store", dest="latlon_file",
                      type=str, help="latlon_file")

    parser.add_option('--create_func', action="store", dest="create_func",
                      type=str, help="create_func")

    parser.add_option('--copy_func', action="store", dest="copy_func",
                      type=str, help="copy_func")

    parser.add_option('--scale_func', action="store", dest="scale_func",
                      type=str, help="scale_func")

    parser.add_option('--ref_file', action="store", dest="ref_file",
                      type=str, help="ref_file")



    (options, args) = parser.parse_args()

    run_job(options.tile, int(options.year), options.file_mosaic, options.ncd_file_pattern,\
            options.hdf_dir, options.out_file, options.latlon_file,\
            options.create_func, options.copy_func,
            options.scale_func, options.ref_file)

    print 'make mosaic is done!!!'