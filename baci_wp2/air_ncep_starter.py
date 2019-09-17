# Starter of upscaling air temperature to MODIS
# Air temperature data are from NCEP reanalysis


import numpy as np
import netCDF4 as nc
import os

year = 2001

# tile = 'h17v07'
# tile = 'h17v08'
# tile = 'h18v07'
tile = 'h18v08'
# tile = 'h21v08'
# tile = 'h21v09'
# tile = 'h22v09'
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

ds = nc.Dataset('/group_workspaces/jasmin2/baci/sigil/air_temperature/2.5_deg/air.%d.nc' % year)

air_new = np.zeros((1200, 1200))

step = 10

for year in xrange(2000, 2017):
    for d in xrange(0, ds.variables['air'].shape[0], step):
        com_str = 'python air_ncep.py --tile %s --year %d --day_from %d --day_to %d' % \
                  (tile, year, d, d + step)

        file_sh = 'sh/air_%s_%d_%d.sh' % (tile, year, d)
        f = open(file_sh, 'w')

        f.write('source ~/virtualenv/bin/activate\n')
        f.write("#BSUB -o output/air_%s_%d_%d.o\n" % (tile, year, d))
        f.write("#BSUB -e output/air_%s_%d_%d.e\n" % (tile, year, d))
        f.write("#BSUB -J air_%s_%d_%d\n" % (tile, year, d))
        f.write("#BSUB -R 'rusage[mem=8000]' -M 8000000\n")
        f.write("#BSUB -W 24:00\n")
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -q short-serial\n")

        com = "bsub < " + file_sh
        f.write(com_str)
        f.write('\n')
        f.close()
        print com
        os.system(com)