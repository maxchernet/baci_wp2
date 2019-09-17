import gp_emulator as gp
import brdf_reg
import numpy as np
import optparse
from collections import OrderedDict
import os

lad = 2
if 'geog' in os.uname()[1]:
    # Local emul. dir
    emul_dir = '/home/ucfamc3/DATA/emulators/semidiscrete_3/lad%d/' % lad
    # emul_dir = '/home/ucfamc3/DATA/emulators/semidiscrete/lad_01/'
else:
    # CEMS dir
    emul_dir = '/group_workspaces/cems/baci/sigil/emulators/semidiscrete_3/lad%d/' % lad


# python create_synth.py --save_path %s --vza %f --vaa %f --sza %f --saa %f --raa %f --year %d --doy %d --site %s

parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(),
                               usage=globals()['__doc__'])

parser.add_option('--save_path', action="store", dest="save_path",
                  type=str, help="save_path")

parser.add_option('--vza', action="store", dest="vza",
                  type=str, help="vza")

parser.add_option('--vaa', action="store", dest="vaa",
                  type=str, help="vaa")

parser.add_option('--sza', action="store", dest="sza",
                  type=str, help="sza")

parser.add_option('--saa', action="store", dest="saa",
                  type=str, help="saa")

parser.add_option('--raa', action="store", dest="raa",
                  type=str, help="raa")

parser.add_option('--year', action="store", dest="year",
                  type=str, help="year")

parser.add_option('--doy', action="store", dest="doy",
                  type=str, help="doy")

parser.add_option('--site', action="store", dest="site",
                  type=str, help="site")

(options, args) = parser.parse_args()

save_path = options.save_path
vza = float(options.vza)
vaa = float(options.vaa)
sza = float(options.sza)
saa = float(options.saa)
raa = float(options.raa)
year = int(options.year)
doy = int(options.doy)
site = options.site

wl_modis = np.array([469.0, 555.0, 645.0, 858.5, 1240.0, 2130.0])

x = OrderedDict()

x['xlai'] = np.exp(-2.5 / 2.)
x['xhc'] = 1.
x['rpl'] = 0.1
x['xkab'] = np.exp(-60. / 100.)
x['scen'] = 0.5
x['xkw'] = np.exp(-0.0176 * 50.)
x['xkm'] = np.exp(-0.002 * 100.)
x['xleafn'] = 1.9
x['xs1'] = 0#1.76
x['xs2'] = 0#0.65


ds = np.load('model/params_%d.npz' % year)
lai_img = ds['lai']
kab_img = ds['kab']
scen_img = ds['scen']

brf_mod = np.zeros((6, lai_img.shape[1], lai_img.shape[2]))
brf_mod_0 = np.zeros((6, lai_img.shape[1], lai_img.shape[2]))

# emul_dir = '/home/ucfamc3/DATA/emulators/semidiscrete_3/lad2/'

# load emulator for nadir view
file_emu_0, nad_sza_0, nad_vza_0, nad_raa_0 = brdf_reg.find_nad_emul(0, 0, 180, emul_dir)
em_0 = gp.MultivariateEmulator (dump=file_emu_0)

#for i, d in enumerate(doys):
# print 'doy:', doys[i], i
file_emu, nad_sza, nad_vza, nad_raa = brdf_reg.find_nad_emul(sza, vza, np.abs(saa - vaa), emul_dir)
# print file_emu
em = gp.MultivariateEmulator (dump=file_emu)

for px in xrange(lai_img.shape[1]):
    for py in xrange(lai_img.shape[2]):

        # x['xlai'] = lai_img[doy-1, px, py]
        # x['xkab'] = kab_img[doy-1, px, py]
        # x['scen'] = scen_img[doy-1, px, py]

        rtm_in = np.array([lai_img[doy-1, px, py], \
                           x['xhc'], x['rpl'], \
                           kab_img[doy - 1, px, py], \
                           scen_img[doy - 1, px, py], \
                           x['xkw'], x['xkm'], x['xleafn'], x['xs1'], x['xs2']])

        brf = em_0.predict(np.array(rtm_in))[0]
        brf_mod_0[:, px, py] = np.array(brf[np.round(wl_modis).astype(int) - 400])

        brf = em.predict(np.array(rtm_in))[0]
        brf_mod[:, px ,py] = np.array(brf[np.round(wl_modis).astype(int) - 400])

if os.path.isdir(save_path + 'ready') == False:
    os.mkdir(save_path + 'ready')
np.savez(save_path + 'ready/synth_%d_%d_%s' % (year, doy, site), brf_mod = brf_mod, brf_mod_0 = brf_mod_0, ang = [vza, vaa, sza, saa, raa])
print 'job', year, doy, site, 'has finished!'