# -*- coding: utf-8 -*-
"""
Created on Tue Nov 08 14:28:12 2016

@author: MaX

upscale MOD11 to 500m and do temporal regularization.
"""
import numpy as np
import matplotlib.pyplot as plt
import osr
import gdal
import os
import sys
# sys.path.append('/home/ucfamc3/install/modis_data/modis_data')
# import modis_reader
import struct
import julday
import datetime as dt
import scipy.optimize
import scipy.special as ssp
import scipy.sparse as sp
import jose_cost as jc
from baci_fs import *
import glob
import optparse
import time
import netCDF4 as nc
import h5py

import pdb

class opt_img(object):
    def __init__(self):
        """
        Constructor
        """

        self.state_grid = np.arange(1, 365, 1)
        
        self.resolution = 1000



    def create_regulariser(self, nx, lag=1):
        """Creates a regulariser with a default lag of 1. Clearly, other lags can be
        specified by changing ``lag``.

        Parameters
        -----------
        nx: int
            The size of the state grid (e.g. the number of time steps)
        lag: int
            The regularisation lag. By default, it is set to 1 (consecutive
            timesteps)

        Returns
        ---------
        A regularisation constraint matrix like in Quaife & Lewis (2010). This is a
        sparse matrix of type "LiL".
        """
        # I = np.diag(np.ones(nx))
        # D = (I - np.roll(I, -lag)).T
        # D2 = D.T.dot(D)
        # Z = np.zeros_like(I)
        # DD0 = np.array([D2, Z, Z]).reshape(nx * 3, nx).T
        # DD1 = np.array([Z, D2, Z]).reshape(nx * 3, nx).T  # HACK! 10*
        # DD2 = np.array([Z, Z, D2]).reshape(nx * 3, nx).T  # HACK! 10*
        # DD = np.array([DD0, DD1, DD2]).reshape((nx * 3, nx * 3))
        #DD = sp.lil_matrix(DD)

        I = np.diag(np.ones(nx))
        DD = np.matrix(I - np.roll(I, 1))

        return DD



    def solve_lin(self, A, B, t, sigma):
        D = self.create_regulariser(A.shape[0])
        regularisation_term = 1 * (D.T.dot(D))
        C = np.diag(sigma.flatten())
        # print 'A.dot(C).dot(A.T)', A.dot(C).dot(A.T) + regularisation_term
        # print 'A.dot(C).dot(B)', A.dot(C).dot(B)
        # aa = A.dot(C).dot(A.T) + regularisation_term
        # bb = B.dot(C).dot(B.T)

        aa = (A*C).dot(A.T) + regularisation_term
        bb = (B*C).dot(B.T)

        x = np.linalg.solve(aa, bb)
        # plt.plot(x)
        print ''



    def get_n_days(self, year):
        n_days = (dt.date(year, 12, 31) - dt.date(year, 1, 1)).days + 1
        date_range = np.array([(dt.date(year, 1, 1) + dt.timedelta(days=i)).strftime('%Y.%m.%d')\
                               for i in range(0, n_days)])
        return n_days, date_range
        

    # *************************************************************************


    def mean_img(self, img, t_in, t_out):
        # print 'img:', img.shape
        # print 'doys', doys.shape
        lst_mean = np.zeros((t_out.shape[0], img.shape[1], img.shape[2]))
        # lst_mean = []
        # print 'lst_mean.shape', lst_mean.shape
        # lst_average = []
        lst_std = []
        # doys_average = []
        for i in xrange(t_out.shape[0]-1):
            ind_doy = np.logical_and(t_in >= t_out[i], t_in < t_out[i + 1])
            if True in ind_doy:
                mm = np.mean(img[ind_doy, :, :], axis=0)
                lst_mean[i, :, :] = mm
                # print mm.shape
                # lst_mean = lst_mean.append(mm)
                # lst_mean = np.dstack([lst_mean, mm])
                # lst_mean = np.concatenate((lst_mean, mm))

                # lst_mean = np.append(lst_mean, mm, axis=0)
                # print lst_mean.shape
                # lst_std = np.append(lst_std, np.std(img[ind_doy]))
                # doys_average = np.append(doys_average, i)
        lst_mean[-1, :, :] = mm
        return lst_mean #np.delete(lst_mean, 0, axis=2)


    # *************************************************************************

    def get_img(self, f_name, px, py, cols, rows, sub_ds_n=0):
        """
        get an array from a geo image.

        Parameters
        -----------
        f_name : str
            file name
        px : int
            x pixel
        py : int
            y pixel
        cols : int
            number of columns
        rows : int
            number of rows
        sub_ds_n : int
            number of a subdataset
        """

        ds = h5py.File(f_name, 'r')
        # img = ds['lst'][:, px:px+cols, py:py+rows]
        # img_qc = ds['lst_qc'][:, px:px + cols, py:py + rows]
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #plt.imshow(ds['lst'][2,:,:])
        #plt.show()
        img = ds['lst'][:, py:py + rows, px:px + cols]
        img_qc = ds['lst_qc'][:, py:py + rows, px:px + cols]

        # ind = img[:,0,0]!=0
        # plt.plot(img[:,0,0][ind])
        # plt.show()

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        doys = ds['doy'][:]

        print 'img.shape:', img.shape
        # print 'doys:', doys

        # gg = gdal.Open(f_name, gdal.GA_ReadOnly)
        # subs = gg.GetSubDatasets()
        #
        # # daily LST is a first subdataset
        # gg = gdal.Open(subs[sub_ds_n][0])
        #
        # s = gg.ReadRaster(xoff=px, yoff=py, xsize=cols, ysize=rows)
        #
        # # all possible format characters for conversion between C and python values
        # c_type = ["h", "H", "B", "b", "i", "I", "l", "L", \
        #              "q", "Q", "f", "d", "s", "p", "P", "c", "?"]
        # c = "h"
        # try:
        #    # first try to open data with the short data type (h)
        #    x = np.array(struct.unpack(cols * rows * c, s))
        # except:
        #    # if above doesnt work try all other datatypes
        #    for c in c_type:
        #        try:
        #            x = np.array(struct.unpack(cols * rows * c, s))
        #            break
        #        except:
        #            continue
        #
        # if c == '?':
        #    print 'type is not found'
        #    return -1
        #
        # img = np.reshape(x, (cols, rows))
        return img, img_qc, doys

        # *************************************************************************



    def get_img_air(self, f_name, px, py, cols, rows, sub_ds_n=0):
            """
            get an array from a geo image.

            Parameters
            -----------
            f_name : str
                file name
            px : int
                x pixel
            py : int
                y pixel
            cols : int
                number of columns
            rows : int
                number of rows
            sub_ds_n : int
                number of a subdataset
            """
            # print 'f_name:', f_name
            # print np.load(f_name).keys()
            img = np.load(f_name)['air_new_interp'][px:px+cols, py:py+rows]

            d = int(f_name.split('.')[0].split('_')[-1][1:])
            y = int(f_name.split('.')[0].split('_')[-2])
            jd = julday.date2jul(julday.doy2date(y, d))

            return img, jd


    # *************************************************************************





    def do_opt(self, obs, obs_sd, date_jul, gamma=70, t = np.arange(1,365,1)):
        """
        Do optimization for given vector of observation obs

        :param obs: vector of observations
        :param obs_sd: observational uncertainties
        :param date_jul: vector of julian days
        :param gamma: dynamical model error (smoothness)
        :param t: time grid
        :return:
        """

        # Define output dates.
        # We want to work with whole years
        year = julday.jul2date2(date_jul[0]).year
        jd1 = julday.date2jul(dt.date(year, 1, 1))
        year = julday.jul2date2(date_jul[-1]).year
        jd2 = julday.date2jul(dt.date(year, 12, 31))
        # t = np.arange(jd1, jd2 + 1, 4)
        t_obs = np.arange(jd1, jd2 + 1, 1)
        x_prior = np.ones_like(t) * np.mean(obs)

        y_obs = np.zeros(t_obs.shape)
        sigma_obs = np.zeros(t_obs.shape)

        for i in xrange(t_obs.shape[0]-1):
            for j in xrange(date_jul.shape[0]):
                if t_obs[i] == date_jul[j]:
                    y_obs[i] = obs[j]
                    sigma_obs[i] = obs_sd[j]
        
        obs_mask = np.zeros(t.shape[0]).astype(bool)

        obs_mask = np.where(sigma_obs > 0, True, False)

        sigma_obs_t = np.zeros(t.shape)

        sigma_prior = 100.

        retval = np.zeros(x_prior.shape[0])
        post_sd = np.zeros(x_prior.shape[0])

        # Sometimes fmin_bfgs() can give "input/output error". The reasons are unknown for me. So if this happens we can
        # simply return zeros what flags this pixel as a bad one.
        try:
            retval = scipy.optimize.fmin_bfgs(jc.cost_function_t, x_prior, fprime=jc.der_cost_function_t, \
                                          args=(y_obs, x_prior, sigma_obs, sigma_prior, gamma, obs_mask, t, t_obs),
                                          disp=1, retall=1)[0]
        except:
            print 'An error in fmin_bfgs():', sys.exc_info()[0]
            return retval, post_sd, sigma_prior, sigma_obs_t

        # Calculate observational uncertainty in grid of output.
        # This is requared for calculation of posteriour uncertainty
        for i in xrange(t.shape[0] - 1):
            ind = np.logical_and(t_obs >= t[i], t_obs < t[i+1])
            k = 0
            for j in xrange(t_obs[ind].shape[0]):
                if y_obs[ind][j] != 0:
                    sigma_obs_t[i] = sigma_obs_t[i] + sigma_obs[ind][j]
                    k += 1
            if k !=0:
                sigma_obs_t[i] = sigma_obs_t[i] / k
                    # der_j_obs[i] = (x[i] - y[ind][j]) ** 2 / sigma_obs[ind][j] ** 2
        # We assume that uncertainty of no observation is very large
        sigma_obs_t[sigma_obs_t == 0] = 999

        # pdb.set_trace()
        # Do uncertainty
        post_sd = jc.posterior_covariance(retval, sigma_prior, sigma_obs_t, gamma)

        return retval, post_sd, sigma_prior, sigma_obs_t

# ***************************************************************************************

    def get_lst_unc (self, img_qc):
        """
        Get uncertainty for LST MOD11

        img_qc : array
            QC array MOD11
        :return: img_err - array with errors in Kelvins
        """

        img_err = np.zeros(img_qc.shape)

        # firs find pixels which were not produced
        ind0 = np.logical_and(img_qc != 2, img_qc != 3)
        # pixels where error <= 1K
        k1 = np.bitwise_and(img_qc, 0b11000000)
        # pixels where error <= 2K
        k2 = np.bitwise_and(img_qc, 0b01000000)
        k3 = np.bitwise_and(img_qc, 0b10000000)
        k4 = np.bitwise_and(img_qc, 0b11000000)
        ind = k1 <= 63
        img_err[ind] = 1
        ind = k2 == 64
        img_err[ind] = 2
        ind = k3 == 128
        img_err[ind] = 3
        ind = k3 == 192
        img_err[ind] = 5
        # assign very high values for pixels which were not produced
        img_err[~ind0] = 999

        return img_err

# ***********************************************************************************************

    def save_netcdf(self, lst, lst_sd, lst_orig, lst_orig_sd, date_jul_orig, date_jul, out_name):#, geo_str, proj_str):
        """
        Save Sentinel-1 backscatter and uncertainty into a netCDF file

        :param bs:
        :param bs_sd:
        :return:
        """

        nc_out = nc.Dataset(out_name, 'w')
        nc_out.description = 'Land Surface temperature with uncertainties'

        nc_out.createDimension('x', lst.shape[2])
        nc_out.createDimension('y', lst.shape[1])
        nc_out.createDimension('time', lst.shape[0])

        nc_out.createDimension('x_obs', lst_orig.shape[2])
        nc_out.createDimension('y_obs', lst_orig.shape[1])
        nc_out.createDimension('time_obs', lst_orig.shape[0])

        nc_out.createVariable('lst', 'f4', ('time', 'x', 'y'), zlib=True)
        nc_out.createVariable('lst_sd', 'f4', ('time', 'x', 'y'), zlib=True)

        nc_out.createVariable('lst_orig', 'f4', ('time_obs', 'x_obs', 'y_obs'), zlib=True)
        nc_out.createVariable('lst_orig_sd', 'f4', ('time_obs', 'x_obs', 'y_obs'), zlib=True)

        nc_out.createVariable('date_jul', 'i', 'time')
        nc_out.createVariable('date_jul_obs', 'i', 'time_obs')

        nc_out.variables['lst'][:] = lst
        nc_out.variables['lst_sd'][:] = lst_sd

        nc_out.variables['lst_orig'][:] = lst_orig
        nc_out.variables['lst_orig_sd'][:] = lst_orig_sd

        print 'nc_out.variables[date_jul][:]: ', nc_out.variables['date_jul'][:].shape
        print 'date_jul.shape: ', date_jul.shape
        nc_out.variables['date_jul'][:] = date_jul
        nc_out.variables['date_jul_obs'][:] = date_jul_orig


# ******************************************************************************************


    def jul_range(self, jd01, jd02, step=1):
        """
        Define output dates.
        We want to work with whole years

        :param jd01:
        :param jd02:
        :return: t
        """
        year = julday.jul2date2(jd01).year
        jd1 = julday.date2jul(dt.date(year, 1, 1))
        year = julday.jul2date2(jd02).year
        jd2 = julday.date2jul(dt.date(year, 12, 31))
        t = np.arange(jd1, jd2 + 1, step)

        return t


# **************************************************************************************


    def do_job(self, ncd_file, mod11_dir, air_dir, mod11_file, site, year0, tile, ulpx, ulpy, lrpx, lrpy, cc, rr,\
               do_cross=True, do_plot = True, step=1):
        """
        Run regularization for a block with size (cc x rr) of a MOD11 image

        :param mod11_dir: str; Directory of MOD11 files
        :param mod11_file: str; Pattern for files search
        :param ulpx: int; upper left x pixel of a block
        :param ulpy: int; upper left y pixel of a block
        :param cc: int; number of  columns
        :param rr: int; number o rows
        :return:
        """

        print 'ulpx, ulpy:', ulpx, ulpy

        # Get lists of files of a current year and +/-


        # f_list = []
        # f_list1 = np.array([file for file in glob.glob(mod11_dir + 'MOD11A1.A%d*.hdf' % (year0 - 1))])
        # f_list2 = np.array([file for file in glob.glob(mod11_dir + mod11_file)])
        # f_list3 = np.array([file for file in glob.glob(mod11_dir + 'MOD11A1.A%d*.hdf' % (year0 + 1))])
        #
        # f_list = np.concatenate((f_list1, f_list2, f_list3))
        # f_list = f_list2

        # get file lists of air temperature
        # air_dir = '/group_workspaces/cems/baci/sigil/air_temperature/1_km/'
        f_list_air = []
        f_list_air_1 = np.array([file for file in glob.glob(air_dir + 'air_new_%s_%d_d*.npz' % (tile, (year0 - 1)))])
        f_list_air_2 = np.array([file for file in glob.glob(air_dir + 'air_new_%s_%d_d*.npz' % (tile, year0))])
        f_list_air_3 = np.array([file for file in glob.glob(air_dir + 'air_new_%s_%d_d*.npz' % (tile, (year0 + 1)))])


        def sort_files(flst):
            """
            sort filesnames in flst in ascending order
            :param flst:
            :return:
            """
            days = []
            for ss in flst:
                days = np.append(days, ss.split('/')[-1].split('_')[-1].split('.')[0].split('d')[-1]).astype(int)
            ind = np.argsort(days)
            flst = flst[ind]
            return flst

        f_list_air_1 = sort_files(f_list_air_1)
        f_list_air_2 = sort_files(f_list_air_2)
        f_list_air_3 = sort_files(f_list_air_3)

        # if we don't have first or last year, assert an empty list to middle year
        f_list_air = f_list_air_2
        if f_list_air_1.shape[0] == 0:
            f_list_air_1 = f_list_air_2
        if f_list_air_3.shape[0] == 0:
            f_list_air_3 = f_list_air_2

        # make a range of julian days from first to last day of year
        jd0 = julday.date2jul(julday.doy2date(year0, 1))
        n_days, date_range = self.get_n_days(year0)
        jd1 = julday.date2jul(julday.doy2date(year0, n_days))
        # input temporal grid (jul days)
        t = self.jul_range(jd0, jd1, step=1)
        # output  temporal grid with specific step
        t_out = self.jul_range(jd0, jd1, step=step)


        fname = mod11_dir + mod11_file % year0
        fname1 = mod11_dir + mod11_file % (year0 - 1)
        fname2 = mod11_dir + mod11_file % (year0 + 1)

        img, img_qc, doys = self.get_img(fname, ulpx, ulpy, cc, rr, sub_ds_n=0)

        # LST array
        # img = np.zeros((t.shape[0], cc, rr))

        # img_qc = np.zeros(img.shape).astype(int)

        img_out = np.zeros((t_out.shape[0], img.shape[1], img.shape[2]))
        img_out_sd = np.zeros(img_out.shape)

        # Go to Kelvins
        img = img * 0.02
        img_qc = img_qc.astype(int)
        # img[doys-1, :, :] = img0 * 0.02
        # img_qc[doys - 1, :, :] = img_qc0

        # jul_days = np.zeros(t.shape[0])
        jul_days = []

        for doy in doys:
            # jul_days[doy - 1] = julday.date2jul(julday.doy2date(year0, doy))
            jul_days = np.append(jul_days, julday.date2jul(julday.doy2date(year0, doy)))

        jul_days[0] = jd0
        jul_days[-1] = jd1

        # array for air temperature
        # img_air = np.zeros((t.shape[0], cc, rr))
        img_air = np.zeros((f_list_air.shape[0], cc, rr))
        jd_air = []


        for i, fname in enumerate(f_list_air):
            data_tmp, dd = self.get_img_air(fname, ulpx, ulpy, cc, rr, sub_ds_n=0)
            # img_air[dd - 1, :, :] = data_tmp
            img_air[i, :, :] = data_tmp
            jd_air = np.append(jd_air, dd)


        # copy temperature from last day of previous year where QC is good

        if os.path.isfile(fname1):
            print 'copy temperature from last day of previous year where QC is good'
            ii = -1
            qc0 = 2
            # tmp_qc = self.get_img(fname1, ulpx, ulpy, cc, rr, sub_ds_n=1).astype(int)
            tmp_img, tmp_qc, tmp_doys = self.get_img(fname1, ulpx, ulpy, cc, rr, sub_ds_n=0)

            while np.logical_or(qc0 == 2, qc0 == 3):
                # tmp = self.get_img(f_list1[ii], ulpx, ulpy, cc, rr, sub_ds_n=1).astype(int)
                qc0 = tmp_qc[ii, cc / 2, rr / 2]
                # img[:, :, 0] = self.get_img(f_list1[ii], ulpx, ulpy, cc, rr, sub_ds_n=0) * 0.02
                img[0, :, :] = tmp_img[ii, :, :] * 0.02
                img_air[0, :, :], jd = self.get_img_air(f_list_air_1[ii], ulpx, ulpy, cc, rr, sub_ds_n=0)
                img_qc[0, :, :] = tmp_qc[ii, :, :]
                ii = ii - 1
                if ii * (-1) >= tmp_img.shape[0]:
                    break


        img_orig = img[:].copy()


        # f_list1 = np.array([])
        # f_list3 = np.array([])

        # LST uncertainty
        img_err = self.get_lst_unc(img_qc)

        # Make an array of uncertainties
        # ind = img_air == 0
        for ii in xrange(jd_air.shape[0]):
            if np.logical_and(0 in img_air[ii,:,:], jd_air[ii] in jul_days):
                img_err[jul_days == jd_air[ii],:,:] = 999

        ind = img_orig == 0
        img_err[ind] = 999

        # ind = img != 0
        ind = img_err != 999

        img_diff = np.zeros(img.shape)
        ind_air = []
        for ii, jj in enumerate(jul_days):
            if jj in jd_air:
                img_diff[ii, :, :] = img_air[jd_air == jj, :, :] - img[ii, :, :]
                ind_air = np.append(ind_air, jj)
            else:
                print jj

        # Go to difference between LST and air temperature
        # img = img_air - img

        # remove water etc from input data
        ind_err = img_err == 999
        img_diff[ind_err] = 0


        # find uncertaities of LST and air difference as percentage of LST uncertainties.
        perc = (img_err * 100) / (np.max(img_orig[ind]) - np.min(img_orig[ind]))
        # img_err_air = ((np.max(img[ind]) - np.min(img[ind])) * perc) / 100.
        img_err_air = img_err.copy()


        #if np.min(img_err) == 999:
        #if img.shape[0] > 0:
            #img_err = np.zeros(img.shape)
            # print 'no data in this block\n'
            # self.save_netcdf(np.moveaxis(img_out, 0, 2), \
            #                  np.moveaxis(img_out_sd, 0, 2), \
            #                  np.moveaxis(img_orig, 0, 2), \
            #                  np.moveaxis(img_err, 0, 2), \
            #                  t, jul_days, ncd_file, geo, proj)
            # return -1

        cx = img_diff.shape[1] / 2
        cy = img_diff.shape[2] / 2

        if do_cross:
            # Cross validation

            img_tmp = img_diff.copy()
            img_tmp[img_err == 999] = np.NaN
            img_mean = np.nanmean(img_tmp, axis=(1,2))
            ind = ~np.isnan(img_mean)

            #ind = img_err[:, cx, cy] != 999
            #y_obs = img[ind, cx, cy]
            y_obs = img_mean[ind]
            # sigma_obs = img_err[ind, cx, cy]
            # sigma_obs = img_err[ind, cx, cy]
            sigma_obs_air = img_err_air[ind, cx, cy]
            date_obs = jul_days[ind]
            # num_off = y_obs.shape[0] * (1 - (30 / float(y_obs.shape[0])))
            num_off = int(y_obs.shape[0] * 0.7)
            ind_cv = np.unique(np.round(np.random.rand(num_off) * (y_obs.shape[0] - 1)).astype(int))
            diff = []
            # gammas = [0.1, 0.5, 1, 10, 20, 50, 70, 100, 200, 500, 1000, 1500, 2000, 3000]
            gammas = [1, 10, 20, 50, 70, 100, 200, 500, 1000, 1500, 2000, 3000]
            for gamma in gammas:

                mask = np.ones(y_obs.shape).astype(bool)
                mask[ind_cv] = False

                print 'gamma=', gamma

                retval_0, post_sd_0, sigma_prior, sigma_obs_t = self.do_opt(y_obs[ind_cv], sigma_obs_air[ind_cv],\
                                                                            date_obs[ind_cv], gamma=gamma, t=t_out)

                ind_in = np.in1d(t, date_obs[mask])
                try:
                    diff = np.append(diff, np.sum((retval_0[ind_in] - y_obs[mask])**2 / 2.))
                except:
                    print 'retval_0[ind_in] and y_obs[mask] are not equal'
                    diff = np.append(diff, 100000)
                print 'difference:', diff[-1], gamma

                # if do_plot:
                #     plt.title('%d, %d, %d, %.4f'%(gamma, y_obs[ind_cv].shape[0], date_obs[mask].shape[0], diff))
                #     plt.fill_between(t, retval_0 + post_sd_0, retval_0 - post_sd_0, color="0.8")
                #     plt.plot(t, retval_0)
                #     plt.plot(t[ind_in], retval_0[ind_in], 'o', c='y')
                #     plt.plot(date_obs[mask], y_obs[mask], '.', c = 'r')
                #     plt.plot(date_obs[ind_cv], y_obs[ind_cv], '.', c='g')
                #     plt.show()

            gamma = gammas[np.argmin(diff)]
            print 'OPtimal gamma = ', gamma
        else:
            gamma = 1


        # retval = np.zeros(t.shape) * -1
        # post_sd = np.zeros(t.shape) * -1
        for i in xrange(img_diff.shape[1]):
            for j in xrange(img_diff.shape[2]):

                ind = np.logical_and(img_err[:, i, j] != 999, jul_days[:] != 0)
                if img_diff[ind, i, j] != []:

                    print 'gamma:', gamma
                    print 'ind:', len(ind)

                    # self.solve_lin(img[ind, i, j], img[ind, i, j], t, img_err[ind, i, j])
                    if np.sum(jul_days[ind]) !=0:
                        retval_0, post_sd_0, sigma_prior, sigma_obs_t = self.do_opt(img_diff[ind, i, j],\
                                                                                    img_err_air[ind, i, j],\
                                                                                    jul_days[ind],\
                                                                                    gamma=gamma,\
                                                                                    t=t_out)
                        retval = retval_0
                        post_sd = post_sd_0

                        img_out[:, i, j] = retval
                        img_out_sd[:, i, j] = post_sd
                    else:
                        img_out[:, i, j] = 0
                        img_out_sd[:, i, j] = 0



        def fill_zeros(arr):
            for ii in xrange(1, arr.shape[0]):
                if np.sum(arr[ii, :, :]) == 0:
                    arr[ii, :, :] = arr[ii - 1, :, :]
            if np.sum(arr[0, :, :]) == 0:
                arr[0, :, :] = arr[1, :, :]
            return arr


        img_air = fill_zeros(img_air)




        #!!!!!!!!!!!
        # Let's get back to LST
        img_air_mean = self.mean_img(img_air, t, t_out)

        img_air_mean = self.mean_img(img_air, jd_air, t_out)
        print 'img_air_mean: ', img_air_mean.shape

        img_air_mean = fill_zeros(img_air_mean)

        img_air_mean[img_out==0] = 0
        img_out = img_air_mean - img_out

        #!!!!!!!!!!!




        # Save retrieved data to netCDF. Change dimensions by np.moveaxis() for geoserver.
        # self.save_netcdf(img_out, img_out_sd, img, img_err, t, jul_days, ncd_file, geo, proj)
        # self.save_netcdf(np.moveaxis(img_out, 0, 2),\
        #                  np.moveaxis(img_out_sd, 0, 2),\
        #                  np.moveaxis(img_orig, 0, 2),\
        #                  np.moveaxis(img_err, 0, 2),\
        #                  jul_days, t_out, ncd_file)


        self.save_netcdf(img_out, img_out_sd, img_orig, img_err, jul_days, t_out, ncd_file)


        # get indices where we have observations.
        # ind = np.logical_and(img_err[:, 0, 0] != 999, jul_days[:] != 0)
        ind = img_err[:, cx, cy] != 999

        do_plot = False
        if np.logical_and(do_plot == True, np.max(ind) != False):
            # plt.fill_between(t, retval+post_sd, retval-post_sd, color="0.8")
            # plt.plot(t, retval, color='Crimson')

            # ind = img_orig == 0
            plt.title('ulpx=%d, ulpy=%d:' % (ulpx, ulpy))
            # plt.fill_between(t_out,\
            #                  img_out[:, 0, 0] + img_out_sd[:, 0, 0], \
            #                  img_out[:, 0, 0] - img_out_sd[:, 0, 0], color="0.8")

            post_sd = img_out_sd[:, cx, cy]
            ci_5 = np.sqrt(2) * ssp.erfinv(0.05)
            ci_25 = np.sqrt(2) * ssp.erfinv(0.25)
            ci_75 = np.sqrt(2) * ssp.erfinv(0.75)

            plt.fill_between(t_out, img_out[:, cx, cy] - ci_75 * post_sd, img_out[:, cx, cy] + ci_75 * post_sd,
                             color='0.9')
            plt.fill_between(t_out, img_out[:, cx, cy] - ci_25 * post_sd, img_out[:, cx, cy] + ci_25 * post_sd,
                             color='0.8')
            plt.fill_between(t_out, img_out[:, cx, cy] - ci_5 * post_sd, img_out[:, cx, cy] + ci_5 * post_sd,
                             color='0.6')



            plt.plot(t_out, img_out[:, cx, cy], color='Crimson', label='Solution and associated uncertainties')


            plt.errorbar(jul_days[ind], img_orig[ind, cx, cy], yerr = img_err[ind, cx, cy],\
                         ecolor='g', marker='.', ls='', label='Original LST and associated uncertainties')

            # plt.errorbar(jul_days[ind], img[ind, 0, 0], yerr=img_err[ind, 0, 0], ecolor='g', marker='o', c='g', ls='')

            plt.plot(jd_air, img_air[:, cx, cy], color='c', marker='.', label='Air Temperature')

            plt.legend(loc=0)

            plt.xlabel('Julian day')
            plt.ylabel(('LST, K'))

            # plt.ylim(-4.35 , -4.3)
            # plt.ylim(280, 310)

            plt.show()
            plt.savefig('fig/therm_opt_test.png')

        print img_diff.shape


# **************************************************************************************************

def run_this_block(ncd_file, mod11_dir, mod11_air_dir, mod11_file, ulpx, ulpy, cc, rr, step):
    t_reg = opt_img()
    t_reg.do_job(ncd_file, mod11_dir, mod11_air_dir, mod11_file, ulpx, ulpy, cc, rr, step)


# **************************************************************************************************
        
if __name__ =="__main__":

    parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(),
                                   usage=globals()['__doc__'])

    parser.add_option('--site', action="store", dest="site",
                      type=str, help="test site")

    parser.add_option('--ncd_out', action="store", dest="ncd_out",
                      type=str, help="output netCDF file")

    parser.add_option('--root', action="store", dest="root",
                      type=str, help="dir where source files are located")

    parser.add_option('--root_air', action="store", dest="root_air",
                      type=str, help="dir where air temperature source files are located")

    parser.add_option('--pattern', action="store", dest="pattern",
                      type=str, help="mod11_file")

    parser.add_option('--ulx', action="store", dest="ulx",
                      type=str, help="upper left x")

    parser.add_option('--uly', action="store", dest="uly",
                      type=str, help="upper left y")

    parser.add_option('--cols', action="store", dest="cols",
                      type=str, help="number of columns")

    parser.add_option('--rows', action="store", dest="rows",
                      type=str, help="number of rows")
                      
    # parser.add_option('--geo', action="store", dest="geo",
    #                   type=str, help="geo transform")
    #
    # parser.add_option('--proj', action="store", dest="proj",
    #                   type=str, help="projection")

    parser.add_option('--tile', action="store", dest="tile",
                      type=str, help="MODIS tile")

    parser.add_option('--year', action="store", dest="year",
                      type=str, help="year")

    parser.add_option('--step', action="store", dest="step",
                      type=str, help="temporal step")



    (options, args) = parser.parse_args()

    print 'site:', options.site
    print 'root:', options.root
    print 'pattern:', options.pattern
    # year = int(options.pattern[9:13])
    # print 'year str:', options.pattern % year[8:12]
    year = int(options.year)
    print 'year:', year
    print 'ulx:', options.ulx
    print 'uly', options.uly
    print 'cols:', options.cols
    print 'rows:', options.rows
    # print 'geo:',  options.geo
    # print 'proj:', options.proj
    print 'step:', options.step

    timeBefore = time.clock()

    # Create an intance of the opt_img class
    t_reg = opt_img()

    # t_reg.do_job(options.ncd_out, options.root, options.root_air, options.pattern, options.site, year, options.tile,\
    #              int(options.ulx), int(options.uly), 0, 0, int(options.cols),
    #                  int(options.rows), options.geo, options.proj, step=int(options.step))

    t_reg.do_job(options.ncd_out, options.root, options.root_air, options.pattern, options.site, year, options.tile, \
                 int(options.ulx), int(options.uly), 0, 0, int(options.cols),
                 int(options.rows), step=int(options.step))

    timeAfter = time.clock()
    print 'total elapsed time (h): %.5f' % ((timeAfter - timeBefore) / (60. * 60.))

    print 'therm_opt_regional done!!!'