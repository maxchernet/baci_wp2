# system libraries
import numpy as np
import os
import glob
import gdal, osr
import stat
import importlib

# User lib
import kernel_opt_regional as kernel_opt
import therm_opt_regional as therm



cols = 60#300#120#240#120#100#1#40#1
rows = 60#300#120#240#120#100#1#40#1
N = 1200#300#1200#1
M = 1200#300#1200#1
ul_px = 0#755#250#0#1027#1033#0#1033#0#210
ul_py = 0#208#100#0#877#893#85#893#0#210

# Number of days per step
step = 7

# tile = 'h17v03'
# tile = 'h17v04'
# tile = 'h17v05'
# tile = 'h18v03'
# tile = 'h18v04'
# tile = 'h18v05'
# tile = 'h19v03'
# tile = 'h19v04'
# tile = 'h19v05'

tile = 'h20v03'
# tile = 'h20v04'
# tile = 'h19v11'
# tile = 'h19v12'
# tile = 'h20v11'
# tile = 'h20v12'
# tile = 'h22v08'

# tiles = ['h20v03', 'h20v04']
# tiles = ['h20v04', 'h19v11', 'h19v12', 'h20v11', 'h20v12', 'h22v08']
# tiles = ['h17v07', 'h17v08', 'h18v07', 'h18v08', 'h21v08', 'h21v09', 'h22v09']
# tiles = ['h17v08', 'h18v07', 'h18v08', 'h21v08', 'h21v09', 'h22v09']
# tiles = ['h21v08', 'h21v09', 'h22v09']
# tiles = ['h18v08']


# tiles = ['h17v07', 'h17v08', 'h18v08']

# 18.04.2018
# tiles = ['h18v07']
# tiles = ['h22v09']
# tiles = ['h21v08', 'h21v09']
# tiles = ['h19v04']
# tiles = ['h19v12', 'h20v12']

# tiles = ['h17v03', 'h17v04', 'h17v05', 'h18v03', 'h18v04', 'h18v05', 'h19v03', 'h19v04', 'h19v05', 'h20v03', 'h20v04']
tiles = ['h19v11', 'h19v12', 'h20v11', 'h20v12']

# product_name = 'mcd11a1'
product_name = 'mod11a1'

for tile in tiles:

    print 'tile:', tile

    if 'geog' in os.uname()[1]:
        # root = '/space/ucfamc3/viterbo/MODIS/MOD11A1/'
        root = '/home/ucfamc3/max_python/data/baci_mod11/'
        out_dir = '/home/ucfamc3/max_python/data/baci_mod11/'
        root_air = '/home/ucfamc3/DATA/air_temperature/1_km/'
    else:
        root = '/group_workspaces/jasmin2/baci/sigil/%s_one_hdf/%s/' % (product_name, tile)
        out_dir = '/group_workspaces/jasmin2/baci/sigil/baci_mod11/'
        root_air = '/group_workspaces/jasmin2/baci/sigil/air_temperature/1_km/%s/' % tile

    opt_py_file = 'therm_opt_regional.py'
    com = ''



    for year in xrange(2015, 2016):
        print 'year: %d, tile: %s' % (year, tile)
        # Pattern for input filename
        pattern = product_name.upper() + '_%d.hdf'
        # product_name = '%s_%d' % (product_name, year)

        # Write a mosaic file with information about produced blocks
        # more convenient to merge all block together
        file_mosaic = out_dir + product_name + '_%d_mosaic_%s_%dday.txt' % (year, tile, step)

        f = open(file_mosaic, 'w')
        f.write('N M ulx uly cols rows\n')
        f.close()

        for ulx in xrange(ul_px, ul_px + N, cols):
            for uly in xrange(ul_py, ul_py + M, rows):


                # add a string to a mosaic file
                fm = open(file_mosaic, 'a+')
                fm.write('%d %d %d %d %d %d\n' % (N, M, ulx, uly, cols, rows))
                fm.close()

                # Output filename
                ncd_out = out_dir + '%s_%s_%d_%d_%d_%dday.nc' % (product_name, tile, year, ulx, uly, step)

                brf = kernel_opt.opt_img()

                lrx = ulx + cols
                lry = uly + rows

                # if np.logical_or(os.path.isfile(ncd_out) == False, 'geog' in os.uname()[1]):
                if np.logical_or(os.path.isfile(ncd_out) == False, os.path.isfile(ncd_out) == True):

                    therm_opt = therm.opt_img()

                    print root + pattern % year

                    # Check whether block is empty. I.e. water, etc. If it's true doesn't necessary to submit it.
                    # img0, img_qc0, doys = therm_opt.get_img(root + pattern % year, ulx, uly, cols, rows, sub_ds_n=0)
                    # if np.max(img0) == 0:
                    #     print 'Empty block'
                    #     continue

                    file_sh = 'sh/%s_%d_%s_%dx%d_%d,%d.sh' % (product_name, year, tile, N, M, ulx, uly)
                    f = open(file_sh, 'w')

                    com_str = 'python lst_regional.py --ncd_out %s --root %s --root_air %s --pattern %s --year %d --tile %s --ulx %d --uly %d --cols %d --rows %d --step %d' % \
                              (ncd_out, root, root_air, pattern, year, tile, ulx, uly, cols, rows, step)

                    if 'geog' in os.uname()[1]:

                        # do_job(self, ncd_file, mod11_dir, air_dir, mod11_file, site, year0, tile, ulpx, ulpy, lrpx, lrpy, cc, rr, geo, proj,do_cross=True, do_plot = True, step=1)
                        therm_opt.do_job(ncd_out, root, root_air, pattern, '', year,\
                                     tile, ulx, uly, 0, 0, cols, rows, '', '', step=step)
                        # os.system(com_str)
                    else:
                        f.write('source ~/virtualenv/bin/activate\n')
                        f.write("#BSUB -o " + out_dir + "output/%s_%d_%s_%dx%d_%d,%d.o\n" %
                                (product_name, year, tile, N, M, ulx, uly))
                        f.write("#BSUB -e " + out_dir + "output/%s_%d_%s_%dx%d_%d,%d.e\n" %
                                (product_name, year, tile, N, M, ulx, uly))
                        f.write("#BSUB -J %dx%d_%d,%d_%s_%d_%s\n" %
                                (N, M, ulx, uly, product_name, year, tile))
                        f.write("#BSUB -R 'rusage[mem=8000]' -M 8000000\n")
                        f.write("#BSUB -W 24:00\n")
                        f.write("#BSUB -n 1\n")
                        f.write("#BSUB -q short-serial\n")
                        com = "bsub < " + file_sh
                    f.write(com_str)
                    f.write('\n')
                    f.close()

                    os.system(com)
                else:
                    print 'file %s already exists' % ncd_out