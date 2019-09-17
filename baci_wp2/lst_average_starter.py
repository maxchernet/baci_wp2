import os
import numpy as np
import gdal
import glob

# *******************************************************************************************

step = 7

# tiles = ['h17v03', 'h17v04', 'h17v05',\
#          'h18v03', 'h18v04', 'h18v05',\
#          'h19v03', 'h19v04', 'h19v05',\
#          'h20v03', 'h20v04']

# tiles = ['h17v03', 'h17v04', 'h17v05', 'h18v03', 'h18v04', 'h18v05', 'h19v03', 'h19v04', 'h19v05', 'h20v03', 'h20v04', \
#          'h19v11', 'h19v12', 'h20v11', 'h20v12', \
#          'h21v08', 'h21v09', 'h22v08', 'h22v09', \
#          'h17v07', 'h17v08', 'h18v07', 'h18v08']

# tiles =  ['h17v07', 'h17v08', 'h18v07', 'h18v08']
# tiles = ['h21v08', 'h21v09', 'h22v08', 'h22v09']

tiles = ['h17v03', 'h17v04', 'h17v05', 'h18v03', 'h18v04', 'h18v05', 'h19v03', 'h19v04', 'h19v05', 'h20v03', 'h20v04', \
         'h19v11', 'h19v12', 'h20v11', 'h20v12']

for region in ['09_south_africa', '10_horn_of_africa', '11_west_africa', '13_europe']:

    path_out = '/group_workspaces/jasmin2/baci/sigil/baci_wp2_files/%s/lst/' % region

    for year in xrange(2000, 2017):

        for tile in tiles:
            path_in = '/group_workspaces/jasmin2/baci/sigil/mcd11a1/%s/' % tile
            file_in = 'MCD11A1_%d.hdf' % year
            file_out = path_out + 'lst_average_%s_%d_%dday.nc' % (tile, year, step)
            ref_dir = '/group_workspaces/jasmin2/baci/sigil/mod11a1/' + tile + '/'
            f_list = np.array([file for file in glob.glob(ref_dir + '*.hdf')])
            ref_file = gdal.Open(f_list[0]).GetSubDatasets()[0][0]
            meter_file = 'latlon/%s_2015_meters.npz' % tile

            if os.path.isfile(path_out + 'lst_%s_%d_%dday.nc' % (tile, year, step)) == False:
                continue

            if os.path.isfile(path_in + file_in):

                file_sh = 'sh/lst_%s_%d.sh' % (tile, year)
                f = open(file_sh, 'w')

                com_str = 'python lst_average.py --file_in %s --file_out %s --year %d --step %d --ref_file %s --meter_file %s' % \
                          (path_in + file_in, file_out, year, step, ref_file, meter_file)
                print com_str

                f.write("#BSUB -o output/lst_%s_%d.o\n" % (tile, year))
                f.write("#BSUB -e output/lst_%s_%d.e\n" % (tile, year))
                f.write("#BSUB -J lst_%s_%d\n" % (tile, year))
                f.write("#BSUB -R 'rusage[mem=8000]' -M 8000000\n")
                f.write("#BSUB -W 24:00\n")
                f.write("#BSUB -n 1\n")
                f.write("#BSUB -q short-serial\n")
                com = "bsub < " + file_sh
                f.write(com_str)
                f.write('\n')
                f.close()
                os.system(com)

                # do_average(path_in=path_in,\
                #            path_out = '/group_workspaces/cems/baci/sigil/baci_wp2_files/13_europe/lst/',\
                #            file_in=file_in,\
                #            year=year, \
                #            tile=tile,\
                #            step=7)
            else:
                print 'file %s does not exist' % file_in