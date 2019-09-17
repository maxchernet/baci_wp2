# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 13:45:05 2016

@author: Maxim Chernetskiy

Make time series of images of Sentinel-1 for the BACI project
We want to find optimal smoothing of Sentinel-1 backscatter
taking uncertainties into account.
Sentinel images are downscaled to MODIS 500m resolution.
"""


import gdal
import osr
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import glob
import os
import datetime as dt
import julday
import scipy.optimize
import scipy.special as ssp
import jose_cost as jc

from baci_fs import *
import netCDF4 as nc
import re
import datetime

import pdb

#import threading
import time

import struct

import optparse
import max_utils


class opt_img(object):

    # def __init__(self):
    #     """
    #     Constructor of the Sentinel-1 time series class
    #     :param site:
    #     :param year:
    #     :return:
    #     """
    #
    #     d1 = julday.date2jul(dt.date(year, 1, 1))
    #     d2 = julday.date2jul(dt.date(year, 12, 31))
    #     self.t = np.arange(d1, d2 + 1)
    #     self.site = site

    def get_back(self, loc, sent1_name, ulpx0, ulpy0, cols, rows):
        """
        Get backscatter from files from specified directory (loc).
        Files must have two bands: backscatter and uncertainty
         Return ndarray arrays of backscatter and uncertainty

        :param loc: directory
        :param sent1_name: name pattern for names to search
        :return: sent1_in: backscatter; sent1_in_sd: uncertanty;
                 date_jul: julian day, date

        """
        #pol = 'vv'
        # fs = baci_fs(site)
        # We need the same size for all images in time series.
        # For this we use one reference image
        # All images must have the same projection and geo. units!
        # if fs.site['ref_file'] == '':
        # if we don't have a reference for this site use as
        # a reference first file in list
        # f_list = np.array([file for file in glob.glob(loc + sent1_name)])
        # else:
        #     f_list = np.array([file for file in glob.glob(loc + fs.site['ref_file'])])

        f_list = np.array([file for file in glob.glob(loc + sent1_name)])

        # Get geo. coord of a reference image
        gds = gdal.Open(f_list[0])
        geo = gds.GetGeoTransform()
        mx0 = ulpx0 * geo[1] + geo[0]
        my0 = ulpy0 * geo[5] + geo[3]

        sent1_in = np.zeros((f_list.shape[0], cols, rows))
        # print  'sent1_in.shape:', sent1_in.shape
        sent1_in_sd = np.zeros((f_list.shape[0], cols, rows))
        sent1_in[:] = -99
        date = []
        date_jul = []
        for i,fname in enumerate(f_list):
            print fname
            gds = gdal.Open(fname)

            # get a pixel of current image by geo. coords of ref. image
            geo = gds.GetGeoTransform()
            ulpx = np.round((mx0 - geo[0]) / geo[1]).astype(int)
            ulpy = np.round((my0 - geo[3]) / geo[5]).astype(int)

            s = gds.ReadRaster(xoff=ulpx, yoff=ulpy, xsize=cols, ysize=rows,
                             band_list=[1, 2])
            # print 'g.RasterCount', g.RasterCount
            tmp = np.array(struct.unpack(gds.RasterCount * cols * rows * "f", s))
            ds = np.reshape(tmp, (gds.RasterCount, cols, rows))

            #if np.max(ds) == -99:
            #    continue

            sent1_in[i, :, :] = ds[0, :, :]
            sent1_in_sd[i, :, :] = ds[1, :, :]


            # sent1_in[:,:,i] = gds.GetRasterBand(1).ReadAsArray()
            # sent1_in_sd[:, :, i] = gds.GetRasterBand(2).ReadAsArray()

            f = fname.split('/')[-1]
            date = np.append(date, dt.date(int(f[:4]), int(f[4:6]), int(f[6:8])))
            date_jul = np.append(date_jul, julday.date2jul(date[-1]))
        # print date
        # print date_jul
        # print 'sent1_in', sent1_in.shape
        # print 'sent1_in_sd', sent1_in_sd.shape

        if date != []:
            ind = np.argsort(date_jul)
            date_jul = date_jul[ind]
            date = date[ind]
            sent1_in = sent1_in[ind,:,:]
            sent1_in_sd = sent1_in_sd[ind,:,:]

        #return 10 * np.log10(sent1_in), 10 * np.log10(sent1_in_sd), date_jul, date
        return sent1_in, sent1_in_sd, date_jul, date


    def do_opt(self, bs, bs_sd, dj_opt, d_opt):
        """
        Find optimal smooth solution for time series of backsctter data
        taking uncertainties into account.

        :param back_scat:
        :param back_scat_sd:
        :param date_jul:
        :param date:
        :return:
        """

        # bs = np.concatenate((sent1_in[i, j, :], sent1_in[i, j, :], sent1_in[i, j, :]))
        # bs_sd = np.concatenate((sent1_in_sd[i, j, :], sent1_in_sd[i, j, :], sent1_in_sd[i, j, :]))
        # dj_opt = np.concatenate((date_jul_1, date_jul, date_jul_2))
        # d_opt = np.concatenate((date_1, date, date_2))


        # date0 = []

        #x0 = np.arange(dj_opt[0] - 1, dj_opt[-1] + 1)

        dj_opt1 = []
        dj_opt2 = []
        d_opt1 = []
        d_opt2 = []
        for jd in dj_opt:
            d = julday.jul2date2(jd).date()
            # print d
            d0 = datetime.date(d.year - 1, d.month, d.day)
            d_opt1 = np.append(d_opt1, d0)
            # print date0[-1]
            dj_opt1 = np.append(dj_opt1, julday.date2jul(d0))

            d0 = datetime.date(d.year + 1, d.month, d.day)
            d_opt2 = np.append(d_opt2, d0)
            dj_opt2 = np.append(dj_opt2, julday.date2jul(d0))

            # date_1 = np.append(date_1, datetime.date())


        back_scat = np.concatenate((bs, bs, bs))
        back_scat_sd = np.concatenate((bs_sd, bs_sd, bs_sd))
        date_jul = np.concatenate((dj_opt1, dj_opt, dj_opt2))
        date = np.concatenate((d_opt1, d_opt, d_opt2))

        # back_scat = bs
        # back_scat_sd = bs_sd
        # date_jul = dj_opt
        # date = d_opt


        year = julday.jul2date2(date_jul[0]).year
        d1 = julday.date2jul(dt.date(year, 1, 1))
        year = julday.jul2date2(date_jul[-1]).year
        d2 = julday.date2jul(dt.date(year, 12, 31))
        x = np.arange(d1, d2 + 1)

        #x = np.arange(date_jul[0] - 1, date_jul[-1] + 1)

        ind = ~np.isnan(back_scat)
        back_scat = back_scat[ind]
        back_scat_sd = back_scat_sd[ind]
        date_jul = date_jul[ind]
        date = date[ind]

        x_prior = np.ones_like(x) * np.mean(back_scat)
        # print x
        y_obs = np.zeros(x.shape)
        sigma_obs = np.zeros(x.shape)
        for i in xrange(x.shape[0]):
            for j in xrange(date_jul.shape[0]):
                if x[i] == date_jul[j]:
                    y_obs[i] = back_scat[j]
                    sigma_obs[i] = back_scat_sd[j]

        obs_mask = np.ones(x.shape[0]).astype(int)
        obs_mask = np.where(y_obs > 0, True, False)

        # print 'obs_mask\n', obs_mask
        # print 'y_obs\n', y_obs

        sigma_prior = 100.
        gamma = 1
        # sigma_obs = 0.015

        # os.chdir('/home/ucfamc3/max_python/baci_brdf/')

        # Do optimization of sentinel data in order to restore gaps
        print 'starting fmin_bfgs with y_obs ', y_obs[obs_mask].shape
        retval = scipy.optimize.fmin_bfgs(jc.cost_function, x_prior, fprime=jc.der_cost_function, \
                                          args=(y_obs[obs_mask], x_prior, sigma_obs[obs_mask], sigma_prior, gamma,
                                                obs_mask), disp=1, retall=15, maxiter=1500)

        bounds = bounds = np.zeros((x_prior.shape[0], 2))
        bounds[:, 0] = 0
        bounds[:, 1] = 0.5



        # retval = scipy.optimize.fmin_l_bfgs_b(jc.cost_function, x_prior, fprime=jc.der_cost_function, \
        #                                   args=(y_obs[obs_mask], x_prior, sigma_obs[obs_mask], sigma_prior, gamma,
        #                                         obs_mask),\
        #                                       bounds = bounds, factr=10000000, pgtol=1e-05, maxiter=1500, disp=1)



        # retval = scipy.optimize.minimize(jc.cost_function, x_prior, method="L-BFGS-B", \
        #                         jac=True, bounds=the_bounds, options=self.optimisation_options)


        # retval = scipy.optimize.minimize(jc.cost_function, x_prior, method="L-BFGS-B", \
        #                                  args=(y_obs[obs_mask], x_prior, sigma_obs[obs_mask], sigma_prior, gamma, obs_mask),\
        #                                  jac=False)
        #
        # self.optimisation_options = {"factr": 1000, \
        #                              "m": 400, "pgtol": 1e-12, "maxcor": 200, \
        #                              "maxiter": 1500, "disp": True}



        # Do uncertainty
        post_sd = jc.posterior_covariance(retval[0], sigma_prior, sigma_obs, gamma, obs_mask)


        ind = []
        for jj in xrange(x.shape[0]):
            if x[jj] in self.t:#x0:
                ind = np.append(ind, jj).astype(int)

        return retval[0][ind], post_sd[ind]


            #def get_ts(self, px, py):




    def do_opt2(self, back_scat, back_scat_sd, date_jul, date, t = np.arange(1, 365, 1)):
        """
        Find optimal smooth solution for time series of backsctter data
        taking uncertainties into account.

        :param back_scat:
        :param back_scat_sd:
        :param date_jul:
        :param date:
        :return:
        """

        # back_scat = bs
        # back_scat = bs_sd
        # date_jul = dj_opt
        # date = d_opt

        year = julday.jul2date2(date_jul[0]).year
        d1 = julday.date2jul(dt.date(year, 1, 1))
        year = julday.jul2date2(date_jul[-1]).year
        d2 = julday.date2jul(dt.date(year, 12, 31))
        x = np.arange(d1, d2 + 1)

        ind = np.logical_and(~np.isnan(back_scat), back_scat > -90)
        back_scat = back_scat[ind]
        back_scat_sd = back_scat_sd[ind]
        date_jul = date_jul[ind]
        date = date[ind]

        x_prior = np.ones_like(x) * np.mean(back_scat)
        y_obs = np.zeros(x.shape)
        sigma_obs = np.zeros(x.shape)
        sigma_obs[:] = 99
        for i in xrange(x.shape[0]):
            for j in xrange(date_jul.shape[0]):
                if x[i] == date_jul[j]:
                    y_obs[i] = back_scat[j]
                    sigma_obs[i] = back_scat_sd[j]

        obs_mask = np.ones(x.shape[0]).astype(int)
        obs_mask = np.where(y_obs != 0, True, False)

        sigma_prior = 100.
        gamma = 1

        # Do optimization of sentinel data in order to restore gaps
        # print 'starting fmin_bfgs with y_obs ', y_obs[obs_mask].shape
        # plt.plot(y_obs[obs_mask])
        # plt.plot(x_prior[obs_mask])
        # plt.show()

        # retval = scipy.optimize.fmin_bfgs(jc.cost_function, x_prior, fprime=jc.der_cost_function, \
        #                                   args=(y_obs[obs_mask], x_prior, sigma_obs[obs_mask], sigma_prior, gamma,
        #                                         obs_mask), disp=1, retall=15, maxiter=1500)



        bounds = np.zeros((x_prior.shape[0], 2))
        bounds[:, 0] = -200
        bounds[:, 1] = 0.5
        if y_obs[obs_mask] != []:
            retval = scipy.optimize.fmin_l_bfgs_b(jc.cost_function, x_prior, fprime=jc.der_cost_function, \
                                              args=(y_obs[obs_mask], x_prior, sigma_obs[obs_mask], sigma_prior, gamma,
                                                    obs_mask),\
                                                  bounds = bounds, factr=10000000, pgtol=1e-05, maxiter=1500,\
                                                  iprint=0)[0]
                                                  # disp=0,

            # Do uncertainty
            try:
                post_sd = jc.posterior_covariance(retval, sigma_prior, sigma_obs, gamma, obs_mask)
            except:
                print 'jc.posterior_covariance exception'
                post_sd = np.ones(x.shape[0])

        else:
            retval = np.ones(x.shape[0])
            post_sd = np.ones(x.shape[0])


        ind = []
        for jj in xrange(x.shape[0]):
            if x[jj] in self.t:#x0:
                ind = np.append(ind, jj).astype(int)

        return retval[ind], post_sd[ind]

# ***********************************************************************************************

    def do_opt3(self, obs, obs_sd, date_jul, gamma=70, t=np.arange(1, 365, 1)):
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
        x_prior = np.ones_like(t) * np.mean(obs[obs != -99])

        y_obs = np.zeros(t_obs.shape)
        sigma_obs = np.zeros(t_obs.shape)

        for i in xrange(t_obs.shape[0] - 1):
            for j in xrange(date_jul.shape[0]):
                if t_obs[i] == date_jul[j]:
                    y_obs[i] = obs[j]
                    sigma_obs[i] = obs_sd[j]

        obs_mask = np.zeros(t.shape[0]).astype(bool)

        obs_mask = np.where(sigma_obs > 0, True, False)

        sigma_prior = 100.
        y_obs[y_obs == -99] = 0
        # pdb.set_trace()
        retval = scipy.optimize.fmin_bfgs(jc.cost_function_t, x_prior, fprime=jc.der_cost_function_t, \
                                          args=(y_obs, x_prior, sigma_obs, sigma_prior, gamma, obs_mask, t, t_obs),
                                          disp=1, retall=1)[0]

        # Calculate observational uncertainty in grid of output.
        # This is requared for calculation of posteriour uncertainty
        sigma_obs_t = np.zeros(t.shape)
        for i in xrange(t.shape[0] - 1):
            ind = np.logical_and(t_obs >= t[i], t_obs < t[i + 1])
            k = 0
            for j in xrange(t_obs[ind].shape[0]):
                if y_obs[ind][j] != 0:
                    sigma_obs_t[i] = sigma_obs_t[i] + sigma_obs[ind][j]
                    k += 1
            if k != 0:
                sigma_obs_t[i] = sigma_obs_t[i] / k
                # der_j_obs[i] = (x[i] - y[ind][j]) ** 2 / sigma_obs[ind][j] ** 2
        # We assume that uncertainty of no observation is very large
        sigma_obs_t[sigma_obs_t == 0] = 999

        # pdb.set_trace()
        # Do uncertainty
        post_sd = jc.posterior_covariance(retval, sigma_prior, sigma_obs_t, gamma)

        return retval, post_sd, sigma_prior, sigma_obs_t



# ***********************************************************************************************



    def save_netcdf(self, bs_orig, bs_orig_sd, date_jul_orig,  bs, bs_sd, date_jul, out_name):
        """
        Save Sentinel-1 backscatter and uncertainty into a netCDF file

        :param bs:
        :param bs_sd:
        :return:
        """

        print 'netCDF name:'
        print out_name
        print ''

        nc_out = nc.Dataset(out_name, 'w')
        nc_out.description = 'Smoothed Sentinel-1 backscatter with uncertainties'

        nc_out.createDimension('x', bs.shape[1])
        nc_out.createDimension('y', bs.shape[2])
        nc_out.createDimension('date', bs.shape[0])

        nc_out.createDimension('x_orig', bs_orig.shape[1])
        nc_out.createDimension('y_orig', bs_orig.shape[2])
        nc_out.createDimension('date_orig', bs_orig.shape[0])

        #nc_out.createDimension('date_orig', bs_orig.shape[2])

        nc_out.createVariable('bs', 'f4', ('date', 'x', 'y'), zlib=True)
        nc_out.createVariable('bs_sd', 'f4', ('date', 'x', 'y'), zlib=True)

        nc_out.createVariable('bs_orig', 'f4', ('date_orig', 'x_orig', 'y_orig'), zlib=True)
        nc_out.createVariable('bs_orig_sd', 'f4', ('date_orig', 'x_orig', 'y_orig'), zlib=True)

        nc_out .createVariable('date_jul', 'i', 'date')
        nc_out.createVariable('date_jul_orig', 'i', 'date')

        nc_out.variables['bs'][:] = bs

        # print 'nc_out.variables[bs][:]: ', nc_out.variables['bs'][:].shape
        # print nc_out.variables['bs'][:] >= 10000

        nc_out.variables['bs'][:][nc_out.variables['bs'][:] >= 10000] = ma.masked
        nc_out.variables['bs_sd'][:] = bs_sd
        nc_out.variables['bs_sd'][:][nc_out.variables['bs_sd'][:] >= 10000] = ma.masked
        nc_out.variables['date_jul'][:] = date_jul

        nc_out.variables['bs_orig'][:] = bs_orig
        nc_out.variables['bs_orig'][:][nc_out.variables['bs_orig'][:] >= 10000] = ma.masked
        nc_out.variables['bs_orig_sd'][:] = bs_orig_sd
        nc_out.variables['bs_orig_sd'][:][nc_out.variables['bs_orig_sd'][:] >= 10000] = ma.masked

        #nc_out.variables['date_jul_orig'][:] = date_jul_orig

        # save geo information and projection
        # nchar_proj = len(proj)
        # nchar_geo = len(geo)
        # nc_out.createDimension('nchar_proj', nchar_proj)
        # nc_out.createDimension('nchar_geo', nchar_geo)
        # nc_out.createVariable('proj', 'S1', ('nchar_proj'))
        # nc_out.createVariable('geo_transform', 'S1', 'nchar_geo')
        # proj_char = nc.stringtochar(np.array([proj], 'S%d' % nchar_proj))
        # geo_char = nc.stringtochar(np.array([geo], 'S%d' % nchar_geo))
        # nc_out.variables['proj'][:] = proj_char
        # nc_out.variables['geo_transform'][:] = geo_char



    '''
    ncd_out, root, pattern, fs.site['name'], int(pattern[9:13]), fs.site['tile'],\
                                       int(ulx), int(uly),\
                                       int(ulx) + int(cc), \
                                       int(uly) + int(rr), int(cc), int(rr), \
                                       new_geo, new_proj


    do_job(self, ncd_file, mod11_dir, mod11_file, site, year0, tile, ulpx, ulpy, lrpx, lrpy, cc, rr, geo, proj,\
               do_cross=False, do_plot = False):

    '''


# **********************************************************************************************************************



    def do_job(self, ncd_out, downsample_dir, file_pattern, year, tile,  ulx, uly, lrpx, lrpy, cc, rr, step=1):

        print 'Starting doing main job with ' + ncd_out

        site =''

        #f_list = np.array([file for file in glob.glob(downsample_dir + sent1_mod_file)])

        file_in_1 = file_pattern % (year, tile)
        # get backscatter for the current year for block (ulx, uly) with size cc/rr (columns rows)
        sent1_in, sent1_in_sd, date_jul, date = self.get_back(downsample_dir, file_in_1, ulx, uly, cc, rr)

        # get first and last julian days for 'year'
        year = julday.jul2date2(date_jul[0]).year
        d1 = julday.date2jul(dt.date(year, 1, 1))
        d2 = julday.date2jul(dt.date(year, 12, 31))
        # get temporal range
        self.t = np.arange(d1, d2 + 1, 1)
        t_out = np.arange(d1, d2 + 1, step)

        head_flag = False
        tail_flag = False

        # get filename for previous year
        file_in_0 = file_pattern % (year - 1, tile)
        f_list = np.array([file for file in glob.glob(downsample_dir + file_in_0)])
        # if data exist get backscatter for previous year
        if len(f_list) != 0:
            sent1_in_0, sent1_in_0_sd, date_jul_0, date_0 = self.get_back(downsample_dir, file_in_0, ulx, uly, cc, rr)
            # get data if we have any good pixels
            if np.max(sent1_in_0[:]) != -99:
                # we need last day for this year therefore looping backward
                for i in xrange(sent1_in_0.shape[0]-1, 0, -1):
                    if np.max(sent1_in_0[i, :, :]) != -99:
                        tmp = np.zeros((sent1_in.shape[0]+1, sent1_in.shape[1], sent1_in.shape[2]))
                        tmp[1:,:, :] = sent1_in[:]
                        tmp[0, :, :] = sent1_in_0[i, :, :]
                        sent1_in = tmp.copy()

                        tmp[1:,:, :] = sent1_in_sd[:]
                        tmp[0, :, :] = sent1_in_0_sd[i,:,:]
                        sent1_in_sd = tmp.copy()

                        head_flag = True

                        break
            # add first date
            date_jul = np.append(d1, date_jul)
            date = np.append(dt.date(year, 1, 1), date)


        # get filename for next year
        file_in_2 = file_pattern % (year + 1, tile)
        f_list = np.array([file for file in glob.glob(downsample_dir + file_in_2)])
        if len(f_list) != 0:
            sent1_in_2, sent1_in_2_sd, date_jul_2, date_2 = self.get_back(downsample_dir, file_in_2, ulx, uly, cc, rr)

            if np.max(sent1_in_2[:]) != -99:
                for i in xrange(0, sent1_in_2.shape[2]):
                    if np.max(sent1_in_2[i, :, :]) != -99:
                        tmp = np.zeros((sent1_in.shape[0]+1, sent1_in.shape[1], sent1_in.shape[2]))
                        tmp[:-1,:, :] = sent1_in[:]
                        tmp[-1, :, :] = sent1_in_2[i, :, :]
                        sent1_in = tmp.copy()

                        tmp[:-1, :, :] = sent1_in_sd[:]
                        tmp[-1, :, :] = sent1_in_2_sd[i, :, :]
                        sent1_in_sd = tmp.copy()

                        tail_flag = True

                        break
            # add last date
            date_jul = np.append(date_jul, d2)
            date = np.append(date, dt.date(year, 12, 31))


        #t = np.arange(date_jul[0] - 1, date_jul[-1] + 1)

        sent1_out = np.zeros((t_out.shape[0], sent1_in.shape[1], sent1_in.shape[2]))
        sent1_out_sd = np.zeros((t_out.shape[0], sent1_in.shape[1], sent1_in.shape[2]))

        if np.max(sent1_in) != -99:
            for i in xrange(sent1_in.shape[1]):
                for j in xrange(sent1_in.shape[2]):
                    #print 'pixel %d %d' % (i,j)
                    # retval, post_sd = self.do_opt2(sent1_in[i, j, :], sent1_in_sd[i, j, :], date_jul, date, t_out)
                    retval, post_sd, sigma_prior, sigma_obs_t = self.do_opt3(sent1_in[:, i, j],\
                                                                             sent1_in_sd[:, i, j],\
                                                                             date_jul, gamma=1,  t=t_out)

                    # sent1_out[i, j, :] = 10 * np.log10(retval)
                    # sent1_out_sd[i, j, :] = 10 * np.log10(post_sd)

                    sent1_out[:, i, j] = retval
                    sent1_out_sd[:, i, j] = post_sd


        if head_flag == True:
            sent1_in = np.delete(sent1_in, 0, axis=0)
            sent1_in_sd = np.delete(sent1_in_sd, 0, axis=0)
            date_jul = np.delete(date_jul, 0)
            date = np.delete(date, 0)

        if tail_flag == True:
            sent1_in = np.delete(sent1_in, -1, axis=0)
            sent1_in_sd = np.delete(sent1_in_sd, -1, axis=0)
            date_jul = np.delete(date_jul, -1)
            date = np.delete(date, -1)

        print 'sent1_in:', sent1_in.shape

        # self.save_netcdf(sent1_out, sent1_out_sd, self.t, 'sentinel/sent1_%d_%d_%d_%d.nc' % (ulx, uly, cc, rr))
        #ind = np.in1d(self.t, date_jul)

        ind = []
        for ii, d0 in enumerate(self.t):
            for d1 in date_jul:
                if d0 == d1:
                    ind = np.append(ind, ii).astype(int)

        bs_orig = np.ones(sent1_in.shape) * 10000
        bs_orig_sd = np.ones(sent1_in.shape) * 10000
        # bs_orig[:,:,ind] = sent1_in
        # bs_orig_sd[:, :, ind] = sent1_in_sd
        bs_orig[:] = sent1_in
        bs_orig_sd[:] = sent1_in_sd

        bs_orig[np.isnan(bs_orig)] = 10000
        bs_orig_sd[np.isnan(bs_orig_sd)] = 10000

        # Save result to netCDF file
        self.save_netcdf(bs_orig, bs_orig_sd, date_jul, sent1_out, sent1_out_sd, t_out, ncd_out)

        print 'run_job series is done!!!'


        do_plot = False
        if do_plot:
            i = 0
            j = 0
            # t = np.arange(date_jul[0] - 1, date_jul[-1] + 1)

            ci_5 = np.sqrt(2) * ssp.erfinv(0.05)
            ci_25 = np.sqrt(2) * ssp.erfinv(0.25)
            ci_75 = np.sqrt(2) * ssp.erfinv(0.75)

            plt.fill_between(t_out, sent1_out[:, i, j] - ci_75 * post_sd, sent1_out[:, i, j] + ci_75 * post_sd,
                             color='0.9')
            plt.fill_between(t_out, sent1_out[:, i, j] - ci_25 * post_sd, sent1_out[:, i, j] + ci_25 * post_sd,
                             color='0.8')
            plt.fill_between(t_out, sent1_out[:, i, j] - ci_5 * post_sd, sent1_out[:, i, j] + ci_5 * post_sd,
                             color='0.6')

            # plt.fill_between(t, retval[0])
            print i,j
            plt.plot(t_out, sent1_out[:, i, j], label='Interpolated backscatter')

            # plt.errorbar(date_jul, 10 * np.log10(sent1_in[i, j, :]), yerr=10 * np.log10(sent1_in_sd[i, j, :]),
            #              marker='o', ls='')

            ind = sent1_in[:, i, j] > -90
            try:
                plt.errorbar(date_jul[ind], sent1_in[ind, i, j], yerr=sent1_in_sd[ind, i, j],
                             marker='o', ls='', label='Observations')
            except:
                print 'can not plot errorbar. sorry.'

            #plt.ylim(-5, -15)
            plt.xlabel('Julian day', fontsize='large')
            plt.ylabel('Backscatter', fontsize='large')
            plt.legend(loc=0, fontsize='large')
            plt.show()




if __name__ == "__main__":


    # downsample_dir, sent1_mod_file, ulx, uly, cc, rr
    #'python sar_opt_regional.py --ncd_out %s --root %s --pattern %s --year %d --tile %s --ulx %d --uly %d --cols %d --rows %d'

    parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(),
                                   usage=globals()['__doc__'])

    parser.add_option('--ncd_out', action="store", dest="ncd_out",
                      type=str, help="output netCDF file")

    parser.add_option('--root', action="store", dest="root",
                      type=str, help="dir where source files are located")

    parser.add_option('--pattern', action="store", dest="pattern",
                      type=str, help="pattern for sent1_mod files")

    parser.add_option('--year', action="store", dest="year",
                      type=str, help="year")

    parser.add_option('--tile', action="store", dest="tile",
                      type=str, help="MODIS tile")

    parser.add_option('--ulx', action="store", dest="ulx",
                      type=str, help="upper left x")

    parser.add_option('--uly', action="store", dest="uly",
                      type=str, help="upper left y")

    parser.add_option('--cols', action="store", dest="cols",
                      type=str, help="number of columns")

    parser.add_option('--rows', action="store", dest="rows",
                      type=str, help="number of rows")

    parser.add_option('--step', action="store", dest="step",
                      type=str, help="temporal step")

    (options, args) = parser.parse_args()

    print 'ddir:', options.root
    print 'smod:', options.pattern
    print 'ulx:', options.ulx
    print 'uly', options.uly
    print 'cols:', options.cols
    print 'rows:', options.rows
    print 'step:', options.step

    timeBefore = time.clock()

    # sent1_ts = opt_img(site=options.site)
    sent1_ts = opt_img()
    # sent1_ts.do_job(options.ncd_out, options.root, options.pattern,\
    #                  int(options.ulx), int(options.uly), int(options.cols), int(options.rows),\
    #                  options.geo, options.proj)

    sent1_ts.do_job(options.ncd_out, options.root, options.pattern, \
               int(options.year), options.tile,\
               int(options.ulx), int(options.uly),\
               int(options.ulx) + int(options.cols),\
               int(options.uly) + int(options.rows),\
               int(options.cols), int(options.rows), int(options.step))

    timeAfter = time.clock()

    print 'total elapsed time (h): %.5f' % ((timeAfter - timeBefore) / (60.*60.))