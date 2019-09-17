import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pdb
import julday
import matplotlib.animation as animation



def get_n_days(year):
    n_days = (datetime.date(year, 12, 31) - datetime.date(year, 1, 1)).days + 1
    date_range = np.array([(datetime.date(year, 1, 1) + datetime.timedelta(days=i)).strftime('%Y.%m.%d') \
                           for i in range(0, n_days)])
    return n_days, date_range


def get_sar_ts(x, y, sar_file):
    ds = nc.Dataset(sar_file)
    ts = ds['bs'][y, x, :]
    ts_sd = ds['bs_sd'][y, x, :]
    date = [julday.jul2date2(d).date() for d in ds['date_jul']]

    return ts, ts_sd, date



def get_sar_img(sar_file, doy):
    ds = nc.Dataset(sar_file)
    return ds['bs'][:, :, doy - 1]



def get_lst_ts(x, y, lst_file):
    ds = nc.Dataset(lst_file)
    ts = ds['lst'][y, x, :]
    ts_sd = ds['lst_sd'][y, x, :]
    date = [julday.jul2date2(d).date() for d in ds['date_jul']]

    return ts, ts_sd, date



def get_lst_img(lst_file, doy):
    ds = nc.Dataset(lst_file)
    return ds['lst'][:, :, doy - 1]



def get_rho_ts(x, y, year_range, band, rho_file):
    rho = []
    rho_sd = []
    unc = []
    jd = []

    errors = [5000, 500, 5000, 5000, 5000, 5000]

    for year in year_range:
        print rho_file % year
        ds = nc.Dataset(rho_file % year)
        tmp = ds.groups['reflectance'].variables['refl_b%d' % (band)][:, y, x]
        try:
            tmp = ds.groups['reflectance'].variables['refl_b%d' % (band)][:, y, x]
            rho = np.append(rho, tmp)
        except:
            print 'year %d is corrupted' % year
            ds.close()
            ds = nc.Dataset(rho_file % (year - 1))
            rho = np.append(rho, ds.groups['reflectance'].variables['refl_b%d' % (band)][:, y, x])
        rho_sd0 = ds.groups['reflectance'].variables['refl_b%d_sd' % (band)][:, y, x]
        jd = np.append(jd, ds.variables['julday'])

        err = 5 * np.arange(-1, 1, 2 / float(ds.groups['reflectance']. \
                                             variables['refl_b%d' % (band)][:, y, x].shape[0])) ** 2

        unc_0 = 1 / (ds.groups['reflectance'].variables['refl_b%d' % (band)][:, y, x] * errors[band]) * err

        ind = np.where(ds.groups['reflectance'].variables['refl_b%d' % (band)][:, y, x] <= 0)
        unc_0[ind] = 1
        unc = np.append(unc, unc_0)
        rho_sd0[ind] = 1
        rho_sd = np.append(rho_sd, rho_sd0)
    dates = np.array([julday.jul2date2(d).date() for d in jd])

    return rho, unc, dates


def get_rho_img(rho_file, doy):
    # print rho_file % year
    ds = nc.Dataset(rho_file)
    # print ds['julday']

    img_rgb = np.zeros((ds.groups['reflectance'].variables['refl_b1'].shape[1], \
                        ds.groups['reflectance'].variables['refl_b1'].shape[2], \
                        3))
    # img_rgb_fwd = np.zeros((ds.groups['reflectance'].variables['refl_b1'].shape[1],\
    #                        ds.groups['reflectance'].variables['refl_b1'].shape[2],\
    #                        3))

    # print ds.groups['reflectance'].variables['refl_b1']
    img_rgb[:, :, 0] = ds.groups['reflectance'].variables['refl_b2'][doy - 1, :, :]
    img_rgb[:, :, 1] = ds.groups['reflectance'].variables['refl_b1'][doy - 1, :, :]
    img_rgb[:, :, 2] = ds.groups['reflectance'].variables['refl_b4'][doy - 1, :, :]

    f_name_fwd = rho_file.split('.nc')[0] + '_fwd.nc'

    ind = np.where(img_rgb > 1)
    img_rgb[ind] = 1
    ind = np.where(img_rgb < 0)
    img_rgb[ind] = 0

    for xx in xrange(1, 1190):
        for yy in xrange(1, 1190):
            tmp = np.zeros((5, 5, 3)) * 1.

            tmp[:, :, :] = img_rgb[xx:xx + 5, yy:yy + 5, :]

            tmp[1:4, 1:4] = 1

            if np.mean(tmp) == 1.:
                img_rgb[xx:xx + 5, yy:yy + 5, :] = 1
    return img_rgb




def get_latlon(fname):
    ds = nc.Dataset(fname)
    return ds['lat'][:], ds['lon'][:]


def get_xy(lat, lon, lat_ref, lon_ref, grid_step):
    lon_str = ['%.1f' % s for s in lon[0, :][::grid_step]]
    lat_str = ['%.1f' % s for s in lat[:, 0][::grid_step]]

    x = np.argmin(np.abs(lon[:] - lon_ref))
    xy = np.unravel_index(x, lon[:].shape)
    print 'lon:', lon[xy[0], xy[1]]
    x = xy[1]
    y = np.argmin(np.abs(lat - lat_ref))
    xy = np.unravel_index(y, lat.shape)
    print 'lat:', lat[xy[0], xy[1]]
    y = xy[0]

    return x, y


def do_img(rho_file, lst_file, sar_file, site_name, year, doy_range, out_name, lat_arr, lon_arr, towns, color, out_dir):
    grid_step = 200
    y_cut = 400
    fs = 20

    fig = plt.figure(figsize=(20, 12))

    # *******************************
    # SAR
    # *******************************

    ax1 = plt.subplot(231)
    # sar_file = '/home/ucfamc3/max_python/data/asar_asc_h18v04_2011_vv_30days.nc'
    bs_img = get_sar_img(sar_file, 0)
    print 'bs_img:', bs_img.shape
    plt.imshow(bs_img[y_cut:, :], vmin=-15, vmax=-5, cmap=plt.cm.Greys)
    cb_sar = plt.colorbar(shrink=0.8, orientation='horizontal')
    cb_sar.set_label(label='Backscatter, Db', fontsize=fs)

    # *******************************
    # LST
    # *******************************

    ax2 = plt.subplot(232)
    lst_img = get_lst_img(lst_file, 0)
    plt.imshow(lst_img[y_cut:, :], vmin=255, vmax=290, cmap=plt.cm.jet)
    cb_lst = plt.colorbar(shrink=0.8, orientation='horizontal')
    cb_lst.set_label(label='LST, K', fontsize=fs)

    # *******************************
    # RHO
    # *******************************

    ax3 = plt.subplot(233)

    img_rgb = get_rho_img(rho_file % year, doy_range[0])
    img_rgb = img_rgb[y_cut:, :, :]
    print 'img_rgb:', img_rgb.shape

    plt.imshow(img_rgb, cmap=plt.cm.Greys)
    cb = plt.colorbar(shrink=0.8, orientation='horizontal')
    cb.set_label(label='Reflectance. False color (NIR,R,G)', fontsize=fs)
    plt.imshow(img_rgb * 5)

    xt = np.arange(img_rgb.shape[1])[::grid_step]
    yt = np.arange(img_rgb.shape[0])[::grid_step]

    lat, lon = get_latlon(rho_file % 2011)

    print 'lat lon:', lat.shape, lon.shape

    lon_str = ['%.1f' % s for s in lon[0, :][::grid_step]]
    lat_str = ['%.1f' % s for s in lat[:, 0][::grid_step]]
    for sub in [1, 2, 3]:
        plt.subplot(2, 3, sub)
        plt.xticks(xt, lon_str, rotation='vertical', fontsize=16)
        plt.yticks(yt, lat_str, fontsize=16)

        for i in [0]:  # xrange(lat_arr.shape[0]):
            x, y = get_xy(lat[y_cut:, :], lon[y_cut:, :], lat_arr[i], lon_arr[i], grid_step)
            plt.plot(x, y, 'o', color=color[i], ms=20, alpha=0.3, markeredgecolor='c')
            plt.plot(x, y, 'o', color=color[i], ms=7, markeredgecolor='c')
            plt.text(x, y, towns[i], fontsize=18)

    x, y = get_xy(lat, lon, lat_arr[i], lon_arr[i], grid_step)

    # *******************************
    # SAR
    # *******************************

    ax4 = plt.subplot(234)

    sar_ts, sar_ts_sd, sar_dates = get_sar_ts(x, y, sar_file)
    print 'xy=', x, y
    print sar_ts

    plt.title('SAR', fontsize=16)
    plt.fill_between(sar_dates, sar_ts + sar_ts_sd, sar_ts - sar_ts_sd, color="0.8")
    plt.plot(sar_dates, sar_ts, c='Navy', lw=2, label='Backscatter')
    plt.ylabel('Db', fontsize=16)
    plt.xticks(fontsize=16, rotation=45)
    plt.grid()
    plt.legend(fontsize=14, loc=0)

    # *******************************
    # LST
    # *******************************

    ax5 = plt.subplot(235)

    lst_ts, lst_ts_sd, lst_dates = get_lst_ts(x, y, lst_file)

    plt.title('LST', fontsize=16)
    plt.plot(lst_dates, lst_ts, c='OrangeRed', lw=2, label='LST')

    plt.fill_between(lst_dates, lst_ts + lst_ts_sd, lst_ts - lst_ts_sd, color="0.8")

    plt.ylabel('T, K', fontsize=16)
    plt.xticks(fontsize=16, rotation=45)

    plt.grid()
    plt.legend(fontsize=14, loc=0)

    # *******************************
    # RHO
    # *******************************

    ax6 = plt.subplot(236)

    plt.title('Optical', fontsize=16)

    rho_ts, rho_ts_sd, rho_dates = get_rho_ts(x, y, np.arange(2011, 2012), band=2, rho_file=rho_file)
    plt.plot(rho_dates, rho_ts, c='m', lw=2, label='Normalized (nadir) reflectance NIR')
    plt.fill_between(rho_dates, rho_ts + rho_ts_sd, rho_ts - rho_ts_sd, color="0.8")

    rho_ts, rho_ts_sd, rho_dates = get_rho_ts(x, y, np.arange(2011, 2012), band=1, rho_file=rho_file)
    plt.plot(rho_dates, rho_ts, c='Crimson', lw=2, label='Normalized (nadir) reflectance Red')
    plt.fill_between(rho_dates, rho_ts + rho_ts_sd, rho_ts - rho_ts_sd, color="0.8")

    plt.ylabel('Reflectance', fontsize=16)
    plt.xticks(fontsize=16, rotation=45)

    plt.grid()
    plt.legend(fontsize=14, loc=0)

    # fig.savefig(out_name, bbox_inches='tight')

    fig.tight_layout()

    img = []
    for doy in doy_range:
        # SAR

        bs_img = get_sar_img(sar_file, doy)
        # print 'bs_img:', bs_img.shape
        im_sar = ax1.imshow(bs_img[y_cut:, :], vmin=-15, vmax=-5, cmap=plt.cm.Greys)

        im_sar1 = ax4.axvline(x=sar_dates[doy - 1], color='DodgerBlue', lw=7, alpha=0.2)
        im_sar11 = ax4.axvline(x=sar_dates[doy - 1], color='DodgerBlue', lw=2)

        # LST

        lst_img = get_lst_img(lst_file, doy)
        im_lst = ax2.imshow(lst_img[y_cut:, :], vmin=255, vmax=314, cmap=plt.cm.jet)

        im_lst1 = ax5.axvline(x=lst_dates[doy - 1], color='DodgerBlue', lw=7, alpha=0.2)
        im_lst11 = ax5.axvline(x=lst_dates[doy - 1], color='DodgerBlue', lw=2)

        # RHO

        img_rgb = get_rho_img(rho_file % year, doy)
        img_rgb = img_rgb[y_cut:, :, :]
        im_rho = ax3.imshow(img_rgb * 5)

        im_rho1 = ax6.axvline(x=rho_dates[doy - 1], color='DodgerBlue', lw=7, alpha=0.2)
        im_rho11 = ax6.axvline(x=rho_dates[doy - 1], color='DodgerBlue', lw=2)

        img.append([im_sar, im_sar1, im_sar11, im_lst, im_lst1, im_lst11, im_rho, im_rho1, im_rho11])

    im_ani = animation.ArtistAnimation(fig, img, interval=10000, repeat_delay=3000, blit=True)
    # Set up formatting for the movie files
    # Writer = animation.writers['ffmpeg']
    Writer = animation.writers['imagemagick']
    writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    im_ani.save(out_dir + 'baci_ani_tst.mp4', writer=writer)



site_name = 'South_Europe'

# band_names = ['Red', 'NIR', 'Blue', 'Green', 'SWIR 1240 nm', 'SWIR 2130 nm']
# data_dir = '/home/ucfamc3/max_python/data/'
data_dir = '/group_workspaces/cems/baci/sigil/baci_wp2_files/'
#data_dir = '../data/'
# out_dir = '/home/ucfamc3/max_python/data/video/'
out_dir = 'video/'

# rho_file = data_dir + 'baci_mod09/rho_h18v04_%d_30day.nc'
# band = 1

# fig = plt.figure(figsize=(16, 6))
# n_bands = 2
# x=100; y=200
# up_lim = [0.2, 0.6]


# Hainich
towns = np.array(['Eisenach', 'Muhlhausen', ''])
lat_arr = np.array([50.98, 51.21, 51.2])
lon_arr = np.array([10.32, 10.45, 10.5])

# Denmark
towns = np.array(['Aarhus', ''])
lat_arr = np.array([56.16, 56.5])
lon_arr = np.array([10.20, 10.0])

# Italy

# Viterbo (town)
#lat = 42.42
#lon = 12.108

# Civitavecchia
#lat = 42.09
#lon = 11.82

# 13_southern_europe/lst/
# 13_southern_europe/sar/
towns = np.array(['Viterbo', 'Civitavecchia'])
# lat_arr = np.array([42.421,  42.0928])
# lon_arr = np.array([12.108, 11.85])

lat_arr = np.array([42.42,  42.0922])
lon_arr = np.array([12.107, 11.796])

# convert first letter of site name to upper case
site_name = site_name[0].upper() + site_name[1:]

do_img(rho_file = data_dir + 'rho_h18v04_%d_30day.nc',\
       lst_file = data_dir + '13_southern_europe/lst/lst_h18v04_2011_30day.nc',\
       sar_file = data_dir + '13_southern_europe/sar/asar_asc_h18v04_2011_vv_30days.nc',\
       site_name = site_name,\
       year = 2011,\
       doy_range = np.arange(1, 12),\
       out_name = 'fig/%s_img.png' % site_name,\
       lat_arr = lat_arr,\
       lon_arr = lon_arr,\
       towns = towns,\
       color = ['c', 'Orange'],\
       out_dir = out_dir)