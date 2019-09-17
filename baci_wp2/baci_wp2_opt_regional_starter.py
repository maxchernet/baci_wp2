import numpy as np
import os
import kernel_opt_regional as kernel_opt

print '*********************************************************************'
print '*********************************************************************'
print '*********************************************************************'
print '*********************************************************************'

cols = 120#1
rows = 120#1
N = 1200#100#1200#1
M = 1200#100#1200#1
ul_px = 0#1120#100#100#250#910 #250#0#280#250  #250 #890  #893  #914
ul_py = 0#60#200#250#100#1073 #100#0#115#100 #100 #1000 #1033 #100#1066
com = ''
step = 7

# African tiles:
# SA: h19v11, h19v12, h20v11, h20v12
# HA: h21v08, h21v09, h22v08, h22v09
# WA: h17v07, h17v08, h18v07, h18v08


# tile = 'h17v03'
# tile = 'h17v04'
# tile = 'h17v05'
# tile = 'h18v03'
# tile = 'h18v04'
# tile = 'h18v05'
# tile = 'h19v03'
# tile = 'h19v04'
# tile = 'h19v05'

# tile = 'h20v03'
# tile = 'h20v04'
# tile = 'h19v11'
# tile = 'h19v12'
# tile = 'h20v11'
# tile = 'h20v12'
# tile = 'h22v08'

# 'h19v11', 'h19v12', 'h20v11', 'h20v12'
# 'h21v08', 'h21v09', 'h22v08', 'h22v09'
# 'h17v07', 'h17v08', 'h18v07', 'h18v08'
# 'h17v07', 'h17v08', 'h18v07', 'h18v08', 'h21v08'

# Check!!: mod09_2014_h17v07_1200x1200_480,0.e

for tile in ['h20v11']:
    for year in xrange(2015, 2016):

        if 'geog' in os.uname()[1]:
            root = '/home/ucfamc3/max_python/data/mod09_1km/'
            out_dir = '/home/ucfamc3/max_python/data/baci_mod09/'
        else:
            root = '/group_workspaces/jasmin2/baci/sigil/mod09_1km/%s/' % tile
            # out_dir = '/group_workspaces/cems/baci/sigil/baci_mod09/'
            out_dir = '/group_workspaces/jasmin2/baci/sigil/baci_mod09_%s/' % tile

        pattern = '%s_%d_*_1km.nc' % (tile, year)
        opt_py_file = 'kernel_opt_regional.py'

        product_name = 'mod09_%d' % year

        if os.path.isdir(out_dir) == False:
            print 'Directory %s does not exist. Make it' % out_dir
            os.mkdir(out_dir)
            os.mkdir(out_dir + 'output/')

        # Write a mosaic file with information about produced blocks
        # more convenient to merge all block together
        file_mosaic = out_dir + product_name + '_mosaic_%s_%dday.txt' % (tile, step)

        # if os.path.isfile(file_mosaic) == False:

        f = open(file_mosaic, 'w')
        f.write('N M ulx uly cols rows\n')
        f.close()

        # else:
        #     f = open(file_mosaic, 'a+')
        #     f.write('\n')

        for ulx in xrange(ul_px, ul_px + N, cols):
            for uly in xrange(ul_py, ul_py + M, rows):

                print 'ulx=%d, uly=%d' % (ulx, uly)

                # add a string to a mosaic file
                fm = open(file_mosaic, 'a+')
                fm.write('%d %d %d %d %d %d\n' % (N, M, ulx, uly, cols, rows))
                fm.close()

                ncd_out = out_dir + 'mod09_%s_%d_%d_%d_%dday.nc' % (tile, year, ulx, uly, step)

                # get instance of the class kernel_opt
                brf = kernel_opt.opt_img()

                # lrx = ulx + cols
                # lry = uly + rows

                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if os.path.isfile(ncd_out) == False:
                # if ncd_out != '':
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    file_sh = 'sh/%s_%s_%dx%d_%d,%d.sh' % (product_name, tile, N, M, ulx, uly)
                    f = open(file_sh, 'w')

                    com_str =\
                    'python kernel_opt_regional.py --ncd_out %s --root %s --pattern %s --year %d --tile %s --ulx %d --uly %d --cols %d --rows %d --step %d' % \
                    (ncd_out, root, pattern, year, tile, ulx, uly, cols, rows, step)

                    if 'geog' in os.uname()[1]:
                        # os.system(com_str)
                        lrx = int(ulx) + int(cols)
                        lry = int(uly) + int(rows)
                        brf.do_job(ncd_out, root, pattern, year, tile, ulx, uly, lrx, lry, cols, rows, step=step, do_adj_years=True)
                    else:
                        f.write('source ~/virtualenv/bin/activate\n')
                        # out_dir_lotus = "/group_workspaces/cems/baci/sigil/baci_%s/output/" % product_name[0:5]
                        # if os.path.isdir(out_dir_lotus) == False:
                        #    os.mkdir(out_dir_lotus)
                        f.write("#BSUB -o " + out_dir + "output/%s_%s_%dx%d_%d,%d.o\n" %
                                (product_name, tile, N, M, ulx, uly))
                        f.write("#BSUB -e " + out_dir + "output/%s_%s_%dx%d_%d,%d.e\n" %
                                (product_name, tile, N, M, ulx, uly))
                        f.write("#BSUB -J %dx%d_%d,%d_%s_%s\n" %
                                (N, M, ulx, uly, product_name, tile))
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

                    # brf.do_job(ncd_out, root, pattern, year, tile, ulx, uly, lrx, lry, cols, rows, do_adj_years=False)

                    # com_str = 'python %s --ncd_out %s --root %s --pattern %s --ulx %d --uly %d --cols %d --rows %d --tile %s --year %s' % \
                    #           (opt_py_file, ncd_out, root, pattern, ulx, uly, cols, rows, tile, year)
                    # print com_str


# /group_workspaces/jasmin2/baci/sigil/baci_mod09_h21v09/mod09_h21v09_2001_0_360_7day.nc