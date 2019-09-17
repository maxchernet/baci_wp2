import numpy as np
#import matplotlib.pyplot as plt
#import netCDF4 as nc
#import gdal
import os
import glob


def starter(file_pattern_path, var, band):
        files = np.array(glob.glob(file_pattern_path))
        print 'files: ', files.shape
        for f in files:
                #        print f
                com_str = 'python baci_aggr.py --file_name %s --var %s --band %s'  %  (f, var, band)
                print com_str

                file_sh = 'sh/baci_aggr.sh'
                f = open(file_sh, 'w')

                f.write("#BSUB -o output/baci_aggr.o\n") 
                f.write("#BSUB -e output/baci_aggr.e\n")
                f.write("#BSUB -R 'rusage[mem=8000]' -M 8000\n")
                f.write("#BSUB -W 24:00\n")
                f.write("#BSUB -n 1\n")
                f.write("#BSUB -q short-serial\n")
                com = "bsub < " + file_sh
                f.write(com_str)
                f.write('\n')
                f.close()
                os.system(com)


dir_main = '/group_workspaces/jasmin2/baci/sigil/baci_wp2_files/'

for b in xrange(1, 8):
        file_pattern = '13_europe/optical/rho_h*_7day.nc'
        var = 'reflectance/refl_b%d' % b
        band = 'band%d' % b
        starter(dir_main + file_pattern, var, band)

file_pattern = '13_europe/lst/lst_h*_7day.nc'
var = 'lst'
band = 'band1'
starter(dir_main + file_pattern, var, band)


file_pattern = '13_europe/sar/*_7day_*.nc'
var = 'bs'
band = 'band1'
starter(dir_main + file_pattern, var, band)
