# system libraries
import numpy as np
import os, sys
import glob
import gdal, osr
import stat
import importlib
import kernel_opt_regional as kernel_opt
import sar_opt_regional as sar

# import matplotlib.pyplot as plt
# gds = gdal.Open('/home/ucfamc3/max_python/data/baci_sar/20120305T212530_VV_modis_h18v04.tif')
# img = gds.GetRasterBand(1).ReadAsArray()
# plt.imshow(img)
# plt.show()

# Max libraries
# import baci_fs
# import max_utils

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

# ****************************************************************************************

if __name__ == "__main__":

    cols = 300#50#1#40#1
    rows = 300#50#1#40#1
    N = 1200#1
    M = 1200#1
    ul_px = 0#927#0#1033#0#399#0
    ul_py = 0#365#0#893#0#578#0

    step = 7

    # tile = 'h18v04'
    # tile = 'h17v04'
    # tile = 'h17v05'

    # tiles = ['h17v03', 'h18v03', 'h19v03', 'h20v03']
    # tiles = ['h19v03']
    # tiles = ['h18v03']
    # tiles = ['h19v04']
    # tiles = ['h19v11']
    # tiles = ['h21v08']
    # tiles = ['h17v07']

    # tiles = ['h19v04', 'h19v05', 'h20v03', 'h20v04']
    # tiles = ['h19v11']

    # tiles = ['h17v03',\
    #          'h19v11', 'h19v12', 'h20v11', 'h20v12', \
    #          'h21v08', 'h21v09', 'h22v08', 'h22v09', \
    #          'h17v07', 'h17v08', 'h18v07', 'h18v08']

    #tiles = ['h17v03', 'h17v04', 'h17v05', 'h18v03', 'h18v04', 'h18v05', 'h19v03', 'h19v04', 'h19v05', 'h20v03', 'h20v04']
    #tiles = ['h19v11', 'h19v12', 'h20v11', 'h20v12']
    # tiles = ['h21v08', 'h21v09', 'h22v08', 'h22v09']
    # tiles = ['h17v07', 'h17v08', 'h18v07', 'h18v08']
    tiles = ['h19v11', 'h19v12', 'h20v11', 'h20v12', 'h21v08', 'h21v09', 'h22v08', 'h22v09', 'h17v07', 'h17v08', 'h18v07', 'h18v08']
    # year = 2012

    # sensor = 'asar'
    sensor = 'sent'
    mode = 'desc'
    # mode = 'asc'
    com = ''

    # modes = ['descending', 'ascending']
    # sensors = ['asar', 'sentinel-1']
    modes = ['ascending', 'descending']
    sensors = ['asar', 'sentinel-1']

    n_jobs = 0

    for sensor in sensors:
        print sensor
        for mode in modes:
            print mode
            for tile in tiles:
                print tile

                if 'geog' in os.uname()[1]:
                    root = '/home/ucfamc3/max_python/data/baci_sar/'
                    out_dir = '/home/ucfamc3/max_python/data/baci_sar_2/'
                    # root_air = '/home/ucfamc3/DATA/air_temperature/1_km/'
                else:
                    # root = '/group_workspaces/jasmin2/baci/sigil/sentinel-1_modis/%s/'%tile
                    # root = '/group_workspaces/jasmin2/baci/sigil/sentinel-1_ascending_modis/%s/' % tile
                    root = '/group_workspaces/jasmin2/baci/sigil/sar_modis/%s_%s_modis/%s/' % (sensor, mode, tile)
                    # root = '/group_workspaces/jasmin2/baci/sigil/asar_descending_modis/%s/' % tile
                    # root = '/group_workspaces/jasmin2/baci/sigil/asar_ascending_modis/%s/' % tile
                    # out_dir = '/group_workspaces/jasmin2/baci/sigil/baci_sar_1/'
                    out_dir = '/group_workspaces/jasmin2/baci/sigil/baci_sar/'
                    # root_air = '/group_workspaces/jasmin2/baci/sigil/air_temperature/1_km/'

                for year in xrange(2016, 2018):
                    # pattern = 'h18v04_%d_*_1km.nc' % year
                    # pattern = 'MOD11A1.A%d*.hdf' % year
                    # pattern = '%d*_modis_%s.tif' % (year, tile)
                    pattern = '%d*_modis_%s.tif'
                    # pattern = '%d*_modis.tif' % year
                    opt_py_file = 'sar_opt_regional.py'

                    product_name = 'sar_%d' % year

                    # Write a mosaic file with information about produced blocks
                    # more convenient to merge all block together
                    file_mosaic = out_dir + 'mosaic_%s_%s_%s_%d_%dday.txt' % (sensor, mode, tile, year, step)

                    f_list = np.array([file for file in glob.glob(root + pattern % (year, tile))])
                    if len(f_list) == 0:
                        print 'files for year %d do not exist' % year
                        continue


                    f = open(file_mosaic, 'w')
                    f.write('N M ulx uly cols rows\n')
                    f.close()

                    for ulx in xrange(ul_px, ul_px + N, cols):
                        for uly in xrange(ul_py, ul_py + M, rows):

                            # add a string to a mosaic file
                            fm = open(file_mosaic, 'a+')
                            fm.write('%d %d %d %d %d %d\n' % (N, M, ulx, uly, cols, rows))
                            fm.close()

                            ncd_out = out_dir + '%s_%s_%s_%d_%d_%d_%dday.nc' % (sensor, mode, tile, year, ulx, uly, step)

                            brf = kernel_opt.opt_img()

                            lrx = ulx + cols
                            lry = uly + rows

                            # when we on the local UCL machine we can overwrite previous files
                            if np.logical_or(os.path.isfile(ncd_out) == False, 'geog' in os.uname()[1]):

                                file_sh = 'sh/%s_%s_%dx%d_%d,%d.sh' % (product_name, tile, N, M, ulx, uly)
                                #print file_sh
                                f = open(file_sh, 'w')

                                com_str = 'python sar_opt_regional.py --ncd_out %s --root %s --pattern %s --year %d --tile %s --ulx %d --uly %d --cols %d --rows %d --step %d' % \
                                          (ncd_out, root, pattern, year, tile, ulx, uly, cols, rows, step)
                                #print com_str
                                if 'geog' in os.uname()[1]:
                                    sar_opt = sar.opt_img()

                                    # ncd_out, downsample_dir, sent1_mod_file, year, tile,  ulx, uly, lrpx, lrpy, cc, rr, geo, proj

                                    sar_opt.do_job(ncd_out, root, pattern, year, tile, ulx, uly, 0, 0, cols, rows, step=step)
                                    # os.system(com_str)
                                else:
                                    f.write('source ~/virtualenv/bin/activate\n')
                                    bsub_out = out_dir + "output/%s_%s_%dx%d_%d,%d.o\n" % (product_name, tile, N, M, ulx, uly)
                                    #print bsub_out
                                    f.write("#BSUB -o " + bsub_out)
                                    bsub_err = out_dir + "output/%s_%s_%dx%d_%d,%d.e\n" % (product_name, tile, N, M, ulx, uly)
                                    #print  bsub_err
                                    f.write("#BSUB -e " + bsub_err)
                                    # f.write("#BSUB -J %dx%d_%d,%d_%s_%s\n" %
                                    #         (N, M, ulx, uly, product_name, tile))
                                    f.write("#BSUB -J %d,%d_%s\n" % (ulx, uly, product_name))
                                    f.write("#BSUB -R 'rusage[mem=8000]' -M 8000\n")
                                    f.write("#BSUB -W 24:00\n")
                                    f.write("#BSUB -n 1\n")
                                    f.write("#BSUB -q short-serial\n")
                                    com = "bsub < " + file_sh
                                f.write(com_str)
                                f.write('\n')
                                f.close()

                                os.system(com)
                                n_jobs += 1
                            else:
                                print 'file %s already exists' % ncd_out

    print 'Number of jobs: ', n_jobs
