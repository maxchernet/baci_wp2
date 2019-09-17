# BACI WP2
# This code cleans defect blocks after first calculation by
# kernel_opt_regional.py or therm_opt_regional.py
# The code deletes blocks which for some reasons were
# not properly calculated.
# This allows to recalculate only these bad blocks


import numpy as np
import os
#import gdal, osr
#import kernel_opt_regional as kernel_opt
import netCDF4 as nc

cols = 40
rows = 40
N = 1200
M = 1200
ul_px = 0
ul_py = 0

product = 'rho'

tile = 'h22v08'
for year in [2008, 2009, 2010, 2014, 2015]:#xrange(2010, 2014):

    if 'geog' in os.uname()[1]:
        # root = '/home/ucfamc3/max_python/data/mod09_1km/'
        out_dir = '/home/ucfamc3/max_python/data/baci_mod09_' + tile + '/'
    else:
        # root = '/group_workspaces/cems/baci/sigil/mod09_1km/'
        if product == 'rho':
            out_dir = '/group_workspaces/cems/baci/sigil/baci_mod09_' + tile + '/'
            pattern = 'mod09_%s_%d_%d_%d.nc'
        if product == 'sar':
            out_dir = '/group_workspaces/cems/baci/sigil/baci_sar/'
            pattern = 'sar_%s_%d_%d_%d.nc'

    # pattern = 'h18v04_%d_*_1km.nc' % year
    # opt_py_file = 'kernel_opt_regional.py'
    #
    # product_name = 'mod09_%d' % year

    for ulx in xrange(ul_px, ul_px + N, cols):
        for uly in xrange(ul_py, ul_py + M, rows):

            ncd_out = out_dir + pattern % (tile, year, ulx, uly)

            lrx = ulx + cols
            lry = uly + rows

            if os.path.isfile(ncd_out):
                try:
                    ds = nc.Dataset(ncd_out)
                    if product == 'rho':
                        data = ds.groups['reflectance'].variables['refl_b1'][:]
                    if product == 'sar':
                        data = ds.variables['bs'][:]
                    print year, type(data), np.max(data), np.min(data)
                except:
                    print 'dataset can not be opened, del:', ncd_out
                    os.system('rm ' + ncd_out)
                try:
                    if np.max(data.filled(-99)) == -99:
                        print 'del:', ncd_out
                        os.system('rm ' + ncd_out)
                except:
                    print 'wrong type'
                print ''
	    else:
		print 'file %s does not exist' % ncd_out
