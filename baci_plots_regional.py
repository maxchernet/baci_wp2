import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
# import baci_fs
# from mpl_toolkits.basemap import Basemap
import pdb
import julday
import datetime

print '*******************************************************************'

# data_dir = '../data/baci_wp2_files/'
data_dir = '/group_workspaces/cems/baci/sigil/baci_wp2_files/13_southern_europe/optical/'
# data_dir = '../data/'

# UK **********************

# Brattleby
lat = 53.316
lon = -0.578

# Wytham
lat2 = 51.775
lon2 = -1.339

# Italy ********************

# Viterbo (town)
lat = 42.42
lon = 12.108

# Civitavecchia
lat = 42.09
lon = 11.82

# Somalia ***********************

# Galkayo
lat = 6.76
lon = 47.45

# Denmark **********************

# Copenhagen
lat = 55.68
lon = 12.57

# Aarhus
lat = 56.16
lon = 10.20

# Hainich ********************************

# Muhlhausen
lat = 51.21
lon = 10.45

# Eisenach
lat = 50.98
lon = 10.32

y = np.arange(-1, 1, 1/366.)**2
print y[150]


def get_n_days(year):
    n_days = (datetime.date(year, 12, 31) - datetime.date(year, 1, 1)).days + 1
    date_range = np.array([(datetime.date(year, 1, 1) + datetime.timedelta(days=i)).strftime('%Y.%m.%d') \
                           for i in range(0, n_days)])
    return n_days, date_range

year = 2015
nd, dr = get_n_days(year)
# dates = np.array([(datetime.date(year, 1, 1) +\
#                                            datetime.timedelta(days=i)) for i in range(0, nd)])

# print data_dir + 'baci_mod09/rho_h18v04_2011_30day.nc'
# ds0 = nc.Dataset(data_dir + 'rho_h18v04_2015_7day.nc')
# print ds0.variables['julday'][:]
# dates = np.array([julday.jul2date2(d).date() for d in ds0['julday']])
# print dates.shape

# ds_fwd0 = nc.Dataset(data_dir + 'rho_h18v04_2015_7day_fwd.nc')

def get_rho(data_dir, lon_arr, lat_arr, year_range, band, file_pattern):
    """

    Parameters
    ----------
    lon_arr
    lat_arr
    year_range
    band
    file_pattern

    Returns
    -------

    """

    rho = []
    rho_sd = []
    unc = []
    dates_fwd = []
    dates_orig = []
    y_fwd = []
    y_orig = []
    dates = []

    errors = [5000, 500, 5000, 5000, 5000, 5000]

    for year in year_range:
        print file_pattern % year
        ds = nc.Dataset(data_dir + file_pattern % year)

        x = np.argmin(np.abs(ds['lon'][:] - lon_arr[0]))
        xy = np.unravel_index(x, ds['lon'][:].shape)
        print 'lon:', ds['lon'][xy[0], xy[1]]
        x = xy[1]
        y = np.argmin(np.abs(ds['lat'][:] - lat_arr[0]))
        xy = np.unravel_index(y, ds['lat'][:].shape)
        print 'lat:', ds['lat'][xy[0], xy[1]]
        y = xy[0]

        print 'x, y = ', x, y

        tmp = ds['reflectance/refl_b%d' % (band + 1)][:, y, x]

        # print 'ds[reflectance/refl_b%d][:, %d, %d]: ' % (band+1, y,x), ds['reflectance/refl_b%d' % (band + 1)][:, y, x].shape
        # print ds['reflectance/refl_b1']
        # print ds['julday']

        # try:
        tmp = ds['reflectance/refl_b%d' % (band + 1)][:, y, x]
        rho = np.append(rho, tmp)
        date_tmp = np.array([julday.jul2date2(d).date() for d in ds['julday']])
        dates = np.append(dates, date_tmp)
        # except:
        #     print 'year %d is corrupted' % year
        #     ds.close()
        #     ds = nc.Dataset(file_pattern % (year - 1))
        #
        #     rho = np.append(rho, ds['reflectance/refl_b%d' % (band + 1)][:, y, x])
        #     dates = np.append(dates, np.array([julday.jul2date2(d).date() for d in ds['julday']]))

        # print 'shapes:', rho.shape, dates.shape

        rho_sd0 = ds['reflectance/refl_b%d_sd' % (band + 1)][:, y, x]

        err = 5 * np.linspace(-1, 1, ds['reflectance/refl_b%d' % (band + 1)][:, y, x].shape[0]) ** 2

        unc_0 = 1 / (ds.groups['reflectance'].variables['refl_b%d' % (band + 1)][:, y, x] * errors[band]) * err

        ind = np.where(ds.groups['reflectance'].variables['refl_b%d' % (band + 1)][:, y, x] <= 0)
        unc_0[ind] = 1
        unc = np.append(unc, unc_0)
        rho_sd0[ind] = 1
        rho_sd = np.append(rho_sd, rho_sd0)

        # print data_dir + 'rho_h18v04_%d_7day_fwd.nc' % year

        # ds_fwd = nc.Dataset(data_dir + 'rho_h18v04_%d_7day_fwd.nc' % year)
        fname_fwd = file_pattern.split('.')[0] + '_fwd.nc'
        print file_pattern.split('/')[-1].split('.')
        print fname_fwd
        ds_fwd = nc.Dataset(data_dir + fname_fwd % year)

        # print ds_fwd['y_fwd']

        ind = np.where(ds_fwd['y_fwd'][:, band, y, x] > 0)
        y_fwd = np.append(y_fwd, ds_fwd['y_fwd'][:, band, y, x][ind])
        dates_fwd = np.append(dates_fwd, date_tmp[ind])

        # print 'y_fwd:', y_fwd.shape

        ind = np.where(ds_fwd['y_orig'][:, band, y, x] > 0)
        y_orig = np.append(y_orig, ds_fwd['y_orig'][:, band, y, x][ind])

        dates_orig = np.append(dates_orig, date_tmp[ind])

        # print dates_fwd

        # print year, rho.shape, dates.shape

    return rho, rho_sd, y_fwd, y_orig, dates, dates_fwd, dates_orig




tile = "h18v04"

site_name = 'South_Europe'

band_names = ['Red', 'NIR', 'Blue', 'Green', 'SWIR 1240 nm', 'SWIR 2130 nm']
colors = ['r', 'm', 'b', 'g', 'Brown', 'SaddleBrown']

towns = np.array(['Viterbo', 'Civitavecchia'])

lat_arr = np.array([42.42,  42.0922])
lon_arr = np.array([12.107, 11.796])

file_pattern = 'rho_' + tile + '_%d_7day.nc'

fig = plt.figure(figsize=(16, 7))
n_bands = 2
x = 100;
y = 200
up_lim = [0.2, 0.5]

plt.title('%s Red-NIR reflectance' % (site_name), fontsize=16)
for band in xrange(n_bands):
    # plt.subplot(n_bands, 1, band + 1)

    rho, rho_sd, y_fwd, y_orig, dates, dates_fwd, dates_orig = get_rho(data_dir,\
                                                                       lon_arr,\
                                                                       lat_arr, \
                                                                       np.arange(2001, 2016), \
                                                                       band=band, \
                                                                       file_pattern=file_pattern)

    # plt.title('%s, %s' % (site_name, band_names[band]), fontsize=16)
    plt.plot(dates, rho, c=colors[band], lw=2, label='Normalized (nadir) reflectance %s' % band_names[band])

    plt.fill_between(dates, rho + rho_sd, rho - rho_sd, color="0.8")



    plt.plot(dates_orig, y_orig, '.', ms=4, label='Original reflectance %s' % band_names[band])
    plt.plot(dates_fwd, y_fwd, '.', ms=4, label='Fwd reflectance %s' % band_names[band])
    plt.ylim(0, up_lim[band])
plt.subplots_adjust(bottom=0.3)
plt.legend(bbox_to_anchor=(1.0, -0.11), fontsize=14, loc=0, ncol=3)


plt.ylabel('Reflectance', fontsize=16)
plt.xticks(fontsize=14)
plt.grid()

# fig.tight_layout()
fig.savefig('fig/%s_%s_red-nir.png' % (tile, towns[0]), bbox_inches='tight')
# fig.savefig('fig/%s_%s_red-nir.png' % (tile, towns[0]))

plt.show()
print dates.shape, rho.shape