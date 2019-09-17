#
# **************************************************************************
# Create time series of tiff images of MODIS normalized reflectance
# using Jose's modis_data tool
# **************************************************************************
#
import os
import glob
import modis_data
import gdal, ogr, osr
import numpy as np
import struct
import sys
import h5py
import netCDF4 as nc
import datetime
#import os
import gp_emulator as gp
import scipy.stats as ss
import scipy.sparse as sp
import scipy.sparse.linalg as sl
import cPickle as pkl
from collections import OrderedDict
from datetime import datetime as dt
import julday
import optparse
import threading

import pdb
import time

import matplotlib.pyplot as plt

host_name = os.uname()[1]
print 'Host is ', host_name
if 'geog' in host_name:
    sys.path.append('/home/ucfamc3/install/brdf_jose/brdf_filter_max/brdf_filter/')
else:
    #if '.rl.' in host_name:
    sys.path.append('/home/users/sigil/install/brdf_jose/brdf_filter_max/brdf_filter/')

    
from kernels import Kernels


# ********************************************************************************

class opt_img(object):
    def     __init__(self):
        """
        Contructor of the class
        """

        # self.threshold = 3.
        self.gamma_min = 3
        self.gamma_max = 12
        self.n_samples = 40
        # self.do_plots = False
        # self.do_albedo = True
        # self.verbose = True

        self.nbands = 7
        self.bu = np.array([0.004, 0.015, 0.003, 0.004, 0.013, 0.010, 0.006])

        # Determine 250 or 500 meters product
        # self.resolution = 500

        # self.pixelWidth = 500
        # self.pixelHeight = 500

    # ********************************************************************************

    def outlier_filter(self, doys, refl, y_fwd, bu, fname, unc=None, threshold=20.,
                       do_plots=True):
        """

        :param doys:
        :param refl:
        :param y_fwd:
        :param bu:
        :param fname:
        :param unc:
        :param threshold:
        :param do_plots:
        :return:
        """
        n_bands = refl.shape[1]
        if do_plots:
            fig, axs = plt.subplots(figsize=(14, 12), nrows=n_bands, ncols=1,
                                    sharex=True)
            axs = axs.flatten()
        pasar = np.ones(doys.shape[0])#y_fwd * 0.
        pasar = pasar.astype(np.bool)
        rpasar = (y_fwd * 0).astype(np.float)
        for band in [0,1,2,3,4,6]: #xrange(n_bands):

            rpasar[doys-1, band] = np.abs(y_fwd[doys-1, band] - refl[:, band]) / self.bu[band]
            # pasar[:, band] = (rpasar[:, band] <= threshold)
            for i in xrange(pasar.shape[0]):
                pasar[i] = np.invert(rpasar[doys-1, band][i] > threshold)

            if do_plots:
                #sel = pasar[:, band]
                sel = pasar
                axs[band].vlines(doys, refl[:, band] - threshold * self.bu[band],
                                 refl[:, band] + threshold * self.bu[band], color="0.8",
                                 alpha=0.8)
                axs[band].plot(doys[sel], refl[sel, band], 'x')
                axs[band].plot(doys[~sel], refl[~sel, band], 'o', mfc="none")

                axs[band].plot(doys[sel], y_fwd[sel, band], '-')
                #pretty_axes(axs[band])

        if do_plots:
            # fig.savefig('pdf/' + fname.replace(".txt", "_outliers.pdf"),
            #             dpi=150, bbox_inches="tight")
            fig.savefig(fname, dpi=150, bbox_inches="tight")

        return rpasar, pasar

    # *******************************************************************************************************


    def create_regulariser (self, nx, lag=1):
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
        I = np.diag (np.ones(nx))
        D = (I - np.roll ( I, -lag)).T
        D2 = D.T.dot(D)
        Z = np.zeros_like (I)
        DD0 = np.array([D2, Z, Z]).reshape ( nx*3, nx).T
        DD1 = np.array([Z, D2, Z]).reshape ( nx*3, nx).T # HACK! 10*
        DD2 = np.array([Z, Z, D2]).reshape ( nx*3, nx).T # HACK! 10*
        DD = np.array( [ DD0, DD1, DD2]).reshape ((nx*3, nx*3))
        DD = sp.lil_matrix( DD )
        return DD

    def predict_refl (self, xsol, doys, vza, sza, raa, rho, do_plots=True, ndoys=365 ):
        """Having an estimate of the kernel weights over an entire period
        of time, this function will simulate observations acquired for the
        same spectral band, stored in rho and with geometries given by
        vza, sza and raa. The doys are the dates of the observations and
        refer to the position of the time step in the xsol array.
        """
        passer = (doys - 1).astype(np.int)
        # Max: ndoys
        f0 = xsol[:ndoys][passer]
        f1 = xsol[ndoys:(ndoys*2)][passer]
        f2 = xsol[(ndoys*2):][passer]



        KK =  Kernels( vza, sza, raa, \
                LiType='Sparse', doIntegrals=False, \
                normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                RossType='Thick' )

        rho_fwd = f0 + KK.Ross[:]*f1 + KK.Li[:]*f2
        if do_plots:
            plt.figure(figsize=(15,10))
            plt.subplot(2,2,1)
            plt.plot(doys, rho_fwd, 'o', label='rho_fwd')
            plt.plot(doys, rho, 'o', label='rho')
            plt.grid()
            plt.legend()

            #plt.figure(figsize=(15,10))
            plt.subplot(2, 2, 2)
            plt.plot(doys, rho- rho_fwd, 'o', label='rho - rho_fwd')
            plt.grid()
            plt.legend()
            #plt.figure(figsize=(15,10))
            plt.subplot(2, 2, 3)
            plt.plot(rho, rho_fwd, 'o', label='rho vs fwd')
            bb = [ rho.min()*0.95, rho.max()*1.1]
            plt.plot( bb, bb, 'k-')
            plt.grid()
            plt.legend()
            plt.tight_layout()
            plt.subplot(2, 2, 4)
            plt.plot(f0, marker='o', label='f0')
            plt.plot(rho, marker='o', label='rho')
            plt.plot(rho_fwd, marker='o', label='rho_fwd')
            plt.legend()
            plt.tight_layout()
            #plt.show()
        # rmse = (rho-rho_fwd).std()
        # rmse = (rho - f0).std()
        # rmse = np.sqrt(((rho - rho_fwd)**2).mean()/rho.shape[0])
        rmse = np.sqrt(((rho - f0) ** 2).mean() / rho.shape[0])
        return rmse, rho_fwd

    # **********************************************************************************************************************

    # **********************************************************************************************************************

    def solve_regularised_problem (self, lambdas, band_unc, doys, vza, sza, raa, rho,
                                   doy_range=None, do_unc=False):

        """Solve the regularised linear kernel BRDF problem. The function solves
        for the kernel weights (f0, f1 and f2) for one or more spectral bands,
        taking into account the inherent uncertainty in the observations, the
        acquisition timing and geometries. The inversion is made more robust with a
        regularisation term (1 lag at the moment, easy to change), and thus the
        value of the regularisation scaling term needs to be provided (you may
        estimate this through e.g. cross-validation). In case you want to
        extrapolate from the observational period, you can also define your temporal
        grid using the ``doy_range`` optional parameter. Finally, the full posterior
        covariance matrix can be calculated (and returned) if the ``do_unc`` option
        is set.

        Parameters
        -----------
        lambdas: float
            The value of the 1 timestep lag regularisation scalar. 0 effecitvely
            shuts of the regularisation, and a very large value will lead to fitting
            a constant value over the entire time series.
        band_unc: array
            The per band uncertainty. This can either be a 1D array, where the per
            band uncertainty is assumed identical for all measurements in that band,
            or you can have a 2D array, where the uncertainty is given by
            observation.
        doys: array
            The time of the observations. Make sure this is an integer.
        vza: array
            The array of view zenith angles (in degrees).
        sza: array
            The array of solar zenith angles (in degrees).
        raa: array
            The array of relative azimuth angles (in degrees).
        rho: array
            The reflectance data, a 2D array, where the first dimension is time, and
            the second is spectral band (e.g. 365*7).
        doy_range: array, optional
            A time array to report the kernel weights on. If not set, the function
            will do this from the first time step to the last in the doys.
        do_unc: boolean, optional
            Whether you want the function to calculate the posterior uncertainty.

        Returns
        ----------
        We return at least three arrays: the solution array ``xsol``, the RMSE value
        for fitting the observations, and the forward modelled observations (these
        are useful diagnostics). If ``do_unc`` is set, the posterior uncertainty
        matrix will also be returned. The MAP value of the state (size
        ``3*len(doy_range)``), where the first chunk is the Isotropic weights, the
        second the volumetric and the third the geometric weights. The second
        parameter is the root mean squared error (RMSE) with the observations per
        band, and finally, the predicted observation reflectance is also returned.
        Additionally, if ``do_unc`` was set to ``True``, a ``n_bands`` element list
        is returned with the posterior covariance  uncertainty matrices.
        """
        # The number of kernels is 3
        n_kernels = 3
        # Figure out how many spectral bands we'll be using...
        if np.ndim ( band_unc ) < 2:
            n_bands = len ( band_unc )
        else:
            n_bands = band_unc.shape[1]
        # The state grid is calculated by default, or the user provided one is used.

        if doy_range is None:
            doy_range = np.arange ( doys.min(), doys.max() + 1 )
            first_day = doys.min()
        else:
            first_day = doy_range.min()

        # Some information that'll be used for array size calculations.
        # Some of observations can be missed
        # num_doys = len(doy_range)
        num_doys = np.max(doy_range)

        all_obs = [list(doys).count(x) \
                        for x in doy_range ]
        all_obs = np.array( all_obs )

        # Calculate the relevant kernels.
        K_obs =  Kernels( vza, sza, raa, \
                LiType='Sparse', doIntegrals=False, \
                normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                RossType='Thick' )
        # Define some matrices (and vectors)
        K = sp.lil_matrix ( ( all_obs.sum(), n_kernels*num_doys ) )
        y = sp.lil_matrix ( (all_obs.sum(), n_bands) )

        cov_mat = np.zeros ( ( all_obs.sum(), n_bands ) )
        # This coming loop iterates over the observations and places the kernels,
        # observations and everything else in their positions inside the sparse
        # matrices defined above. A nice candidate for JITing....
        i_obs = -1

        for ( i, doy ) in enumerate ( doy_range ):
            iloc = np.where(doys==doy)[0]
            if doy in doys: # We have an obs
                for ix in iloc:
                    i_obs = i_obs + 1

                    K [ i_obs, doy - first_day ]  = 1.0 # Isotropic
                    K [ i_obs, doy - first_day + num_doys ] = \
                                        K_obs.Ross[ix]
                    # print i
                    K [ i_obs, doy - first_day + num_doys*2 ] = \
                                        K_obs.Li[ix]

                    y[i_obs, :] = rho[ix]
                    if np.ndim ( band_unc ) < 2:
                        cov_mat[i_obs, :] = (1./band_unc )**2
                    elif np.ndim ( band_unc ) == 2:
                        # Sigma is defined per observation
                        cov_mat[ i_obs, : ] = \
                                (1./band_unc [i_obs, :] )**2

        # We now create the regulariser matrix...
        D = self.create_regulariser ( num_doys )

        # Define the outputs...
        x_sol = np.zeros ( ( 3*num_doys, n_bands ) )
        y_fwd = np.zeros ( ( num_doys, n_bands ) )
        y_orig = np.zeros((num_doys, n_bands))

        if np.ndim ( y_fwd ) == 1:
            y_fwd = y_fwd [ :, np.newaxis]
        rmse = np.zeros ( n_bands )

        # Do a per band inversion in this next loop

        # print 'lambdas:', lambdas.shape
        # print '(D.T.dot(D)):', (D.T.dot(D)).shape
        regularisation_term = lambdas*(D.T.dot(D))
        # print 'lambdas:', lambdas.shape, type(lambdas)
        # print np.diag(lambdas).shape
        # regularisation_term = np.diag(lambdas) * (D.T.dot(D))
        # print 'regularisation_term:', regularisation_term.shape

        unc = []

        # !!! Max
        cov_mat_full = []

        for band in xrange ( n_bands ):
            covariance_matrix = sp.lil_matrix ( ( all_obs.sum(),
                        all_obs.sum()))
            covariance_matrix.setdiag ( cov_mat[:,band] )
            # print 'K.T*covariance_matrix*K:', (K.T*covariance_matrix*K).shape
            A = (K.T*covariance_matrix*K + regularisation_term).tocsc()
            #A = (K.T*covariance_matrix*K).tocsc()
            b = ((K.T*covariance_matrix)*y[ :, band]).tocsc()

            retval = sp.linalg.spsolve ( A, b )
            
            # We want only positive values
            # Here I do it by brute force but should be an other way of doing this.
            ind = np.where(retval < 0)[0]
            retval[ind] = 0.001

            timeBefore = time.clock()
            if do_unc:
                # sometimes rather in random way we can get a error "Factor is exactly singular"
                try:
                    unc.append ( sp.linalg.inv ( A ) )
                except:
                    print "Factor is exactly singular????"
                    return np.array([[-1],[-1]]), [-1], [-1], [-1], [-1], [-1]

                # !!! Max
                # I commented it in hope that it can speed the process up
                # cov_mat_full.append(A)

            timeAfter = time.clock()
            #print 'price of unc. calc.:', timeAfter - timeBefore

            passer = (doys - 1).astype(np.int)

            x_sol [ :, band ] = retval
            y_fwd[ passer, band] = np.array ( K*retval ).squeeze()

            # ndoys = doy_range.shape[0]
            # f0 = x_sol[:ndoys, band][passer]
            # f1 = x_sol[ndoys:(ndoys * 2), band][passer]
            # f2 = x_sol[(ndoys * 2):, band][passer]

            f0 = x_sol[:num_doys, band][passer]
            f1 = x_sol[num_doys:(num_doys * 2), band][passer]
            f2 = x_sol[(num_doys * 2):, band][passer]

            y_fwd[ passer, band] = f0 + K_obs.Ross[:] * f1 + K_obs.Li[:] * f2

            rmse[band] = ( y[ :, band ].todense().squeeze() -
                        y_fwd[ passer, band ]).std()

            y_orig[passer, band] = y[ :, band ].todense().squeeze()

        if do_unc:
        # !!! Max: return cov_mat_full
            return x_sol, rmse, y_fwd, y_orig, unc, cov_mat_full
        else:
            return x_sol, rmse, y_fwd, y_orig



    # **********************************************************************************************************************



    def solve_regularised_problem_edge(self, lambdas, band_unc, doys, vza, sza, raa, rho,
                                  doy_range=None, do_unc=False):

        # The number of kernels is 3
        n_kernels = 3
        # Figure out how many spectral bands we'll be using...
        if np.ndim(band_unc) < 2:
            n_bands = len(band_unc)
        else:
            n_bands = band_unc.shape[1]
        # The state grid is calculated by default, or the user provided one is used.

        if doy_range is None:
            doy_range = np.arange(doys.min(), doys.max() + 1)
            first_day = doys.min()
        else:
            first_day = doy_range.min()

        # Some information that'll be used for array size calculations.
        # Some of observations can be missed
        # num_doys = len(doy_range)
        num_doys = np.max(doy_range)

        all_obs = [list(doys).count(x) \
                   for x in doy_range]
        all_obs = np.array(all_obs)

        # Calculate the relevant kernels.
        K_obs = Kernels(vza, sza, raa, \
                        LiType='Sparse', doIntegrals=False, \
                        normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                        RossType='Thick')
        # Define some matrices (and vectors)
        K = sp.lil_matrix((all_obs.sum(), n_kernels * num_doys))
        y = sp.lil_matrix((all_obs.sum(), n_bands))

        cov_mat = np.zeros((all_obs.sum(), n_bands))
        # This coming loop iterates over the observations and places the kernels,
        # observations and everything else in their positions inside the sparse
        # matrices defined above. A nice candidate for JITing....
        i_obs = -1

        for (i, doy) in enumerate(doy_range):
            iloc = np.where(doys == doy)[0]
            if doy in doys:  # We have an obs
                for ix in iloc:
                    i_obs = i_obs + 1

                    K[i_obs, doy - first_day] = 1.0  # Isotropic
                    K[i_obs, doy - first_day + num_doys] = \
                        K_obs.Ross[ix]
                    # print i
                    K[i_obs, doy - first_day + num_doys * 2] = \
                        K_obs.Li[ix]

                    y[i_obs, :] = rho[ix]
                    if np.ndim(band_unc) < 2:
                        cov_mat[i_obs, :] = (1. / band_unc) ** 2
                    elif np.ndim(band_unc) == 2:
                        # Sigma is defined per observation
                        cov_mat[i_obs, :] = \
                            (1. / band_unc[i_obs, :]) ** 2

        # We now create the regulariser matrix...
        D = self.create_regulariser(num_doys)

        # Define the outputs...
        x_sol = np.zeros((3 * num_doys, n_bands))
        y_fwd = np.zeros((num_doys, n_bands))
        y_orig = np.zeros((num_doys, n_bands))

        if np.ndim(y_fwd) == 1:
            y_fwd = y_fwd[:, np.newaxis]
        rmse = np.zeros(n_bands)

        # Do a per band inversion in this next loop

        # print 'lambdas:', lambdas.shape
        # print '(D.T.dot(D)):', (D.T.dot(D)).shape
        # regularisation_term = lambdas*(D.T.dot(D))
        # print 'lambdas:', lambdas.shape, type(lambdas)
        # print np.diag(lambdas).shape


        unc = []

        # !!! Max
        cov_mat_full = []

        for band in xrange(n_bands):

            regularisation_term = sp.lil_matrix(np.diag(lambdas[:,band]) * (D.T.dot(D)))
            # print 'regularisation_term:', regularisation_term.shape

            covariance_matrix = sp.lil_matrix((all_obs.sum(),
                                               all_obs.sum()))
            covariance_matrix.setdiag(cov_mat[:, band])
            # print 'K.T*covariance_matrix*K:', (K.T * covariance_matrix * K).shape
            A = (K.T * covariance_matrix * K + regularisation_term).tocsc()
            # A = (K.T*covariance_matrix*K).tocsc()
            b = ((K.T * covariance_matrix) * y[:, band]).tocsc()

            retval = sp.linalg.spsolve(A, b)

            # We want only positive values
            # Here I do it by brute force but should be an other way of doing this.
            ind = np.where(retval < 0)[0]
            retval[ind] = 0.001

            timeBefore = time.clock()
            if do_unc:
                # sometimes rather in random way we can get a error "Factor is exactly singular"
                try:
                    unc.append(sp.linalg.inv(A))
                except:
                    print "Factor is exactly singular????"
                    return np.array([[-1], [-1]]), [-1], [-1], [-1], [-1], [-1]

                    # !!! Max
                    # I commented it in hope that it can speed the process up
                    # cov_mat_full.append(A)

            timeAfter = time.clock()
            # print 'price of unc. calc.:', timeAfter - timeBefore

            passer = (doys - 1).astype(np.int)

            x_sol[:, band] = retval
            y_fwd[passer, band] = np.array(K * retval).squeeze()

            # ndoys = doy_range.shape[0]
            # f0 = x_sol[:ndoys, band][passer]
            # f1 = x_sol[ndoys:(ndoys * 2), band][passer]
            # f2 = x_sol[(ndoys * 2):, band][passer]

            f0 = x_sol[:num_doys, band][passer]
            f1 = x_sol[num_doys:(num_doys * 2), band][passer]
            f2 = x_sol[(num_doys * 2):, band][passer]

            y_fwd[passer, band] = f0 + K_obs.Ross[:] * f1 + K_obs.Li[:] * f2

            rmse[band] = (y[:, band].todense().squeeze() -
                          y_fwd[passer, band]).std()

            y_orig[passer, band] = y[:, band].todense().squeeze()

        if do_unc:
            # !!! Max: return cov_mat_full
            return x_sol, rmse, y_fwd, y_orig, unc, cov_mat_full
        else:
            return x_sol, rmse, y_fwd, y_orig


    # **********************************************************************************************************************


    def guess_gamma_opt (self, rho, sza, vza, raa, bu, doys,\
                         gamma_min=1, gamma_max=12, n_samples=41,\
                         do_plots=True, verbose=True, doy_range=np.arange(1, 366) ):

        #terra_doys, terra_sza, terra_vza, terra_raa, terra_refl, terra_bu =  \
        #        read_refl_data ( fname,  sensor="terra" )

        #aqua_doys, aqua_sza, aqua_vza, aqua_raa, aqua_refl, aqua_bu =  \
        #        read_refl_data ( fname, sensor="aqua" )
        # Max: ( 1, 366 ) ->  ( 1, ndoys+1 )
        # doy_range = np.arange ( 1, ndoys+1 )

        # Make random samples for cross validation
        # ii = np.unique((np.random.rand(doys.shape[0]/0.9) * (doys.shape[0]-1)).astype(int))

        dist = ss.uniform(loc=0, scale=doys.shape[0])
        ii = np.unique(np.squeeze(gp.lhd(dist=dist, size=int(doys.shape[0]/1.5))).astype(int))

        ind = np.zeros(doys.shape[0]).astype(bool)
        ind[ii] = True
        terra_refl = rho[ind]
        terra_doys = doys[ind]
        terra_sza = sza[ind]
        terra_vza = vza[ind]
        terra_raa = raa[ind]
        terra_bu = self.bu
        ind = np.invert(ind)
        aqua_refl = rho[ind]
        aqua_doys = doys[ind]
        aqua_sza = sza[ind]
        aqua_vza = vza[ind]
        aqua_raa = raa[ind]
        aqua_bu = self.bu

        # plt.plot(terra_doys, terra_refl[:, 1], marker='o', ls='--', c='r')
        # plt.plot(aqua_doys, aqua_refl[:, 1], marker='o', ls='--', c='b')
        # plt.show()


        #nr, nc = n_subplots ( n_samples )
        n_bands = 7
        xx = np.zeros((8,n_samples))
        # if do_plots:
        #     fig1, axs1 = plt.subplots ( figsize=(14, 12), nrows=nr, ncols=nc,
        #             sharex=True )
        #     fig2, axs2 = plt.subplots ( figsize=(14, 12), nrows=nr, ncols=nc,
        #             sharex=True )
        #     axs1 = axs1.flatten()
        #     axs2 = axs2.flatten()
        if verbose:
            print "%-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s" % ( "gamma",
                    "B01", "B02", "B03","B04", "B05", "B06", "B07" )

        for i, gammas in enumerate( np.logspace( gamma_min, gamma_max, n_samples) ):

            #try:
            x_sol, rmse, y_fwd, y_orig = self.solve_regularised_problem(gammas, terra_bu,
                                                    terra_doys, terra_vza, terra_sza,
                                                    terra_raa, terra_refl,
                                                    doy_range=doy_range)
            #except:
            #    print 'regularised problem in guess_gamma_opt has not solved'
            #    continue

            # if do_plots:
            #     # Max: x_sol[:365, 1] -> x_sol[:ndoys, 1]
            #     axs1[i].plot(doy_range, x_sol[:doy_range.shape[0], 1], '-')
            #     axs1[i].set_title ( "%8.2E" % gammas )
            #     pretty_axes ( axs1[i])
            sys.stdout.flush()
            this_line = "%8.2E" % gammas
            xx[0,i] = gammas
            for band in xrange(7):
                # Max: ndoys=ndoys
                rmse, rho_fwd = self.predict_refl( x_sol[:,band], aqua_doys,
                                   aqua_vza, aqua_sza, aqua_raa,
                                   aqua_refl[:,band], do_plots=False, ndoys=np.max(doy_range))#doy_range.shape[0] )
                this_line = this_line + " %5.2E" % rmse
                xx[band+1, i] = rmse
                # if do_plots and band == 1:
                #     axs2[i].plot ( aqua_doys, aqua_refl[:,band], 'o')
                #     axs2[i].plot ( aqua_doys, rho_fwd, '-')
                #     axs2[i].set_title ( "%8.2E (RMSE: %5.1E)" % (gammas, rmse),
                #                       fontsize=8)
                #     pretty_axes ( axs2[i])
            if verbose:
                print this_line


        gamma_opt = xx[0, xx[np.array([1,2,3,4,5,7]),:].argmin(axis=1)].min()
        # print "Chose %g" % ( gamma_opt )

        opt_rmse = np.sum(xx[np.array([1,2,3,4,5,7]),:].min(axis=1))

        return gamma_opt, opt_rmse

    # **********************************************************************************************************************



    #
    # ************************************************************************************
    #
    def save_m_years(self, fnames, site='somalia'):
        #
        # Save a file which consists of multiple years
        #
        dd = np.zeros((1, 15))
        for i in range(fnames.shape[0]):
                d = np.loadtxt(fnames[i])
                d[:, 0] = d[:, 0] + dd[-1, 0]
                dd = np.concatenate((dd, d))
        fn = 'txt/MODIS_%s_%s-%s.txt' % (site, os.path.splitext(fnames[0])[0][-4:],\
                                        os.path.splitext(fnames[-1])[0][-4:])
        dd = np.delete(dd, 0, axis=0)
        np.savetxt(fn, dd, fmt='%.4f',\
                    header=open(fnames[0], 'r').read().split('\n')[0])
        ndoys = dd[-1, 0].astype(int)
        doy_range = np.arange(dd[0, 0].astype(int), ndoys + 1)
        return fn, doy_range

    #
    # ************************************************************************
    #
    def get_n_days(self, year):
        n_days = (datetime.date(year, 12, 31) - datetime.date(year, 1, 1)).days + 1
        date_range = np.array([(datetime.date(year, 1, 1) + datetime.timedelta(days=i)).strftime('%Y.%m.%d')\
                               for i in range(0, n_days)])
        return n_days, date_range

    #
    # ************************************************************************
    #
    def make_vrt(self, tile, year, loc, dir_vrt):

        # Create virtual datasets based of MODIS images

        m = modis_data.MODISFinder( tile=tile, year=year, head_loc=loc)
        m.create_virtual_data(dir_vrt)
    #
    # ************************************************************************
    #
    def get_time_series(self, tile, year, loc, px, py, cols=1, rows=1, step=1, save_dir='txt/'):
        """
        Read time series of reflectance data and save to a pkl file

        Parameters
        ------------------
        tile: MODIS tile
        loc: directory of data
        """

        # keys = ['sza', 'saa', 'vza', 'vaa', 'qa', 'b01', 'b02', 'b03', 'b04', 'b05', 'b06', 'b07']

        timeBefore = time.clock()

        f_list = np.array(glob.glob(loc + '%s_%d_*_1km.nc' % (tile, year)))

        # if there are no files for this year
        if f_list.shape[0] == 0:
            print 'Data for year %d and tile %s not found' % (year, tile)
            return -1

        # get DOYs from file name
        doys0 = np.sort([int(s.split('_')[4]) for s in f_list])
        # print 'doys0', doys0

        # arrange DOYs according to step. I.e. we assume that...
        # step = 7
        doys = np.zeros(doys0.shape).astype(int)
        doys_real = np.zeros(doys0.shape).astype(int)

        for jj, ii in enumerate(xrange(0, doys0.shape[0]+1, step)):
            # doys = np.append(doys, doys0[ii:ii+7])
            doys[ii:ii+step] = jj + 1#doys0[ii]
            doys_real[ii:ii+step] = doys0[ii]

        ind = np.argsort([int(s.split('_')[4]) for s in f_list])
        f_list = f_list[ind]

        # print 'loc:', loc
        # print 'tile:', tile
        # print 'year:', year
        # print 'f_list:', f_list[0]

        output = {}
        for f  in f_list:
            print 'f:', f
            try:
                ds = nc.Dataset(f)
            except:
                print 'something wrong with %s' % f
                continue

            # if a dataset is empty
            # i.e. if it has less than 2 bands
            # i.e. vza, vaa, sza, saa, proj, geo, qa, b1, b2
            if len(ds.variables.keys()) < 9:
                print 'Dataset exists but empty:'
                print 'loc:', loc
                print 'tile:', tile
                print 'year:', year
                print 'f_list:', f_list[0]
                print ''
                # return -1
            else:
                break
        if f == f_list[-1]:
            print 'all datasets are empty'
            return -1

        for key in ds.variables.keys():
            if len(ds.variables[key].shape) == 2:
                # output[key] = np.zeros((f_list.shape[0], ds.variables[key].shape[0], ds.variables[key].shape[1]))
                output[key] = np.zeros((f_list.shape[0], cols, rows))
            if len(ds.variables[key].shape) == 1:
                output[key] = np.zeros(ds.variables[key].shape[0]).astype(str)

        for i, fname in enumerate(f_list):
            #print fname
            ds = nc.Dataset(fname)
            for key in ds.variables.keys():
                if len(ds.variables[key].shape) == 2:

                    try:
                        # print ds.variables[key][px:px+cols, py:py+rows]
                        output[key][i, :, :] = ds.variables[key][px:px+cols, py:py+rows]
                    except:
                        print 'something wrong in output[%s][%d, :, :]' % (key, i)
                        print 'output:', output[key][i, :, :].shape
                        print 'ds.variables:', ds.variables[key][px:px+cols, py:py+rows].shape

                if len(ds.variables[key].shape) == 1:
                    output[key][:] = ds.variables[key][:]

        # print 'output.keys:', output.keys()

        QA_OK = np.array([8, 72, 136, 200, 1288, 2056, 2120, 2184, 2248])
        qa_passer = np.logical_or.reduce([output['qa'] == x for x in QA_OK])

        for b in [1,2,3,4,5,7]:

            qa_passer[output["b0%d" % b] <= 0] = 0.
            qa_passer[output["b0%d" % b] >= 10000] = 0.

            # if, for this pixel, we have just a few observations we don't need them
            if np.sum(qa_passer) < 2:
                qa_passer[:] = 0

            output['qa_passer'] = qa_passer
            [bin(b) for b in QA_OK]

            #doys = np.array([int(g.GetRasterBand(b+1).GetMetadata()['DoY']) for b in xrange(g.RasterCount)])
            #years = np.array([int(g.GetRasterBand(b + 1).GetMetadata()['Year']) for b in xrange(g.RasterCount)])
        output['doys'] = doys
        output['doys_real'] = np.unique(doys_real)
        output['years'] = np.ones(doys.shape) * year

        output["sza"] = output["sza"] / 100.
        output["saa"] = output["saa"] / 100.
        output["vza"] = output["vza"] / 100.
        output["vaa"] = output["vaa"] / 100.
        output['b01'] = output['b01'] / 10000.
        output['b02'] = output['b02'] / 10000.
        output['b03'] = output['b03'] / 10000.
        output['b04'] = output['b04'] / 10000.
        output['b05'] = output['b05'] / 10000.
        output['b06'] = output['b06'] / 10000.
        output['b07'] = output['b07'] / 10000.

        # print 'qa_passer:', output['qa_passer']
        #print output['b01']

        #pkl.dump(output, open(f_out, 'wb'))

        timeAfter = time.clock()
        elapsed_time = timeAfter - timeBefore
        print 'read time series time (s): ', elapsed_time

        print 'Read MODIS for year %d data is done' % year
        return output

    #
    #****************************************************************************
    #

    def get_rho(self, d, x, y):
        """
        Gets
        :param fname:
        :param x:
        :param y:
        :return:
        """

        rho = np.zeros((1, 7))
        doys = []; years = []; sza = []; vza = []; raa = []

        #for i, fname in enumerate(fnames):
            #d = pkl.load(open(fname, 'r'))

        sza = d['sza'][:, y, x][d['qa_passer'][:, y, x]]
        vza = d['vza'][:, y, x][d['qa_passer'][:, y, x]]
        raa = d['saa'][:, y, x][d['qa_passer'][:, y, x]] - d['vaa'][:, y, x][d['qa_passer'][:, y, x]]
        rho = np.zeros((sza.shape[0], 7))
        for b in xrange(7):
            rho[:, b] = d['b0%d' % (b+1)][:, y, x][d['qa_passer'][:, y, x]]
        doys = d['doys'][d['qa_passer'][:, y, x]]
        years = d['years'][d['qa_passer'][:, y, x]]

        #print 'doys:', doys.shape
        #print 'years: ', years.shape
        #print 'y, x:', y, x
        #print 'qa_passer: ', d['qa_passer'][:, y, x].shape
        ind = np.array(np.where(d['qa_passer'][:, y, x] == 1))
        # if len(ind[0]) > 0:
        #     print 'qa_passer true ', ind[0].shape
        # else:
        #     print 'all pixels are not good', ind

        # Make DOYs continuous
        # doys[years == years[0] + 1] = doys[years == years[0] + 1] + 365
        # doys[years == years[0] + 2] = doys[years == years[0] + 2] + 365 * 2
        # doys[years == years[0] + 1] = doys[years == years[0] + 1] + np.max(d['doys'])
        # doys[years == years[0] + 2] = doys[years == years[0] + 2] + np.max(d['doys']) * 2
        doys[years == d['years'][0] + 1] = doys[years == d['years'][0] + 1] + np.max(d['doys'])
        doys[years == d['years'][0] + 2] = doys[years == d['years'][0] + 2] + np.max(d['doys']) * 2


        # vza = np.append(vza, vza0)
        # sza = np.append(sza, sza0)
        # raa = np.append(raa, raa0)
        # # doys = np.append(doys, doys0 + 365 * i)
        # doys = np.append(doys, doys0 + 365)
        # years = np.append(years, years0)
        # rho = np.append(rho, rho0, axis=0)
        # rho = np.delete(rho, 0, 0)

        return rho, doys, years, sza, vza, raa


    # ***************************************************************************


    def plot_results(self, x_sol, y_fwd, bhr_spectral_nbar, bhr_spectral_nbar_unc, bhr_spectral, rho, doy_range, doys):
        n_doys = x_sol.shape[0] / 3

        doy_range = np.arange(np.min(doy_range), np.max(doy_range) + 1)

        doy_range2 = doy_range[365 : 365 + bhr_spectral_nbar[:, 1].shape[0]]

        plt.figure(figsize=(15,10))
        plt.subplot(2,2,1)
        plt.title('NIR')
        # plt.fill_between(doy_range2,\
        #                  bhr_spectral_nbar[:, 1] - bhr_spectral_nbar_unc[:, 1], \
        #                  bhr_spectral_nbar[:, 1] + bhr_spectral_nbar_unc[:, 1], color="0.9")
        # plt.plot(doy_range2, bhr_spectral_nbar[:, 1], c='c', lw=2)
        # plt.plot(doy_range2, bhr_spectral[:, 1], c='g', lw=2)

        ind = np.logical_and(doy_range >= self.min_doy, doy_range < self.max_doy)
        plt.plot(doy_range[ind], x_sol[:n_doys, 1][ind], c='c', lw=5)


        plt.plot(doy_range, x_sol[:n_doys, 1], c='r', lw=2)

        plt.plot(doys, rho[:, 1], marker='.', ls='')
        plt.ylim(-0.1, np.max(y_fwd[:,1])+0.15)
        plt.grid()

        plt.subplot(2,2,2)
        plt.title('Vis')
        # plt.fill_between(doy_range2, \
        #                  bhr_spectral_nbar[:, 0] - bhr_spectral_nbar_unc[:, 0], \
        #                  bhr_spectral_nbar[:, 0] + bhr_spectral_nbar_unc[:, 0], color='r', alpha=0.2)
        # plt.plot(doy_range2, bhr_spectral_nbar[:, 0], c='r', lw=2)
        # plt.fill_between(doy_range2, \
        #                  bhr_spectral_nbar[:, 2] - bhr_spectral_nbar_unc[:, 2], \
        #                  bhr_spectral_nbar[:, 2] + bhr_spectral_nbar_unc[:, 2], color='b', alpha=0.2)
        # plt.plot(doy_range2, bhr_spectral_nbar[:, 2], c='b', lw=2)
        # plt.fill_between(doy_range2, \
        #                  bhr_spectral_nbar[:, 3] - bhr_spectral_nbar_unc[:, 3], \
        #                  bhr_spectral_nbar[:, 3] + bhr_spectral_nbar_unc[:, 3], color='g', alpha=0.2)
        # plt.plot(doy_range2, bhr_spectral_nbar[:, 3], c='g', lw=2)

        plt.plot(doy_range[ind], x_sol[:n_doys, 0][ind], c='c', lw=5)
        plt.plot(doy_range, x_sol[:n_doys, 0], c='Orange', lw=2)

        plt.plot(doys, rho[:, 0], marker='.', ls='', c='r')
        # plt.plot(doys, rho[:, 2], marker='x', ls='', c='b')
        # plt.plot(doys, rho[:, 3], marker='x', ls='', c='g')


        plt.ylim(-0.01, np.max(y_fwd[:, np.array([0,2,3])])+0.06)
        plt.grid()

        plt.subplot(2,2,3)
        plt.title('NDVI')
        # plt.plot(doy_range, (x_sol[:(n_doys), 1] - x_sol[:(n_doys), 0])/(x_sol[:(n_doys), 1] + x_sol[:(n_doys), 0]), c='k', lw=2)
        # plt.plot(doy_range2,\
        #          (bhr_spectral_nbar[:, 1] - bhr_spectral_nbar[:, 0]) / (bhr_spectral_nbar[:, 1] + bhr_spectral_nbar[:, 0]),
        #          c='k', lw=2)
        #plt.plot(doy_range, bhr_spectral[:, 3], c='g')
        # plt.plot(doy_range, bhr_spectral[:, 1], c='m')
        #plt.plot(doys, rho[:, 1], marker='o', ls='--', c='r')
        plt.plot(doys, (rho[:,1] - rho[:,0])/(rho[:,1] + rho[:,0]), '--')
        #plt.plot(doys, rho[:, 3], marker='o', ls='--', c='g')
        plt.ylim(0, 1)
        plt.grid()

        plt.subplot(2,2,4)
        plt.title('Kernels weights NIR')
        plt.plot ( doy_range, x_sol[:(n_doys), 1], '-', lw=2, c='c', label = 'Isometric')
        plt.plot ( doy_range, x_sol[(n_doys):(2 * n_doys), 1], '--', lw=2, label='Volumetric')
        plt.plot ( doy_range, x_sol[(2 * n_doys):, 1], ':', lw=2, label='Geometric')
        plt.legend()

        #plt.ylim(-0.3,1.3)
        plt.grid()
        plt.show()
        # plt.savefig('fig/' + fname.split('/')[-1].split('.')[0] + '.png')


    # *********************************************************************************************************************


    def create_netcdf(self, netcdf_file, cols, rows, doy_range, do_zlib=True):
        """
        Create an netCDF file where we store the results

        Parameters
        ----------
        netcdf_file: str
            file to write

        Returns
        --------
        rootgrp: netCDF4
            a pointer to a netcdf file
        """

        # least_significant_digit for createVariable(...)
        lsd=None

        # create a new netCDF file
        rootgrp = nc.Dataset(netcdf_file, "w")
        rootgrp.description = "NBAR reflectance and BB albedo"

        # Create groups for NBAR reflectance and BB albedo
        refl_grp = rootgrp.createGroup("reflectance")
        refl_grp.description = "Normalized at nadir (NBAR) reflectance"
        albedo_grp = rootgrp.createGroup("albedo")
        albedo_grp.description = "Broad band (BB) albedo"

        # Create dimensions for reflectance data. I.e. time series of reflectance
        # are determined by x, y and a day since 2000
        x = rootgrp.createDimension("x", cols)
        y = rootgrp.createDimension("y", rows)
        day = rootgrp.createDimension("day", doy_range.shape[0])
        band = rootgrp.createDimension("band", 7)
        str_dim = rootgrp.createDimension("str_dim", 10)

        # We can set zlib=True for compression as an argument
        # of the createVariable function
        date_str = rootgrp.createVariable("date_str", "S1", ("day", "str_dim"), zlib=do_zlib, least_significant_digit=lsd)
        date_str.units = "string representation of date: yyyy.mm.dd"
        # date = rootgrp.createVariable("julday", "i")
        date = rootgrp.createVariable("julday", "i", "day", zlib=do_zlib, least_significant_digit=lsd)
        date.units = "Julian day"

        # Create variables of NBAR reflectance: 7 MODIS bands
        for i in range(1, 8):
            refl = refl_grp.createVariable("refl_b%d" % i, "f4", ("day", "x", "y"), zlib=do_zlib,\
                                           least_significant_digit=lsd)
            refl.units = "surface bidirectional reflectance, band %d" % i

            # reflectance uncertainty
            refl_sd = refl_grp.createVariable("refl_b%d_sd" % i, "f4", ("day", "x", "y"), zlib=do_zlib,\
                                              least_significant_digit=lsd)
            refl_sd.units = "uncertainty of surface bidirectional reflectance, band %d" % i
            
        y_fwd = refl_grp.createVariable("y_fwd", "f4", ("day", "band", "x", "y"), zlib=do_zlib,\
                                        least_significant_digit=lsd)

        y_orig = refl_grp.createVariable("y_orig", "f4", ("day", "band", "x", "y"), zlib=do_zlib, \
                                        least_significant_digit=lsd)
        
        albedo_vis = albedo_grp.createVariable("albedo_vis", "f4", ("day", "x", "y"), zlib=do_zlib,\
                                               least_significant_digit=lsd)
        albedo_vis.units = "broad band albedo"
        albedo_nir = albedo_grp.createVariable("albedo_nir", "f4", ( "day", "x", "y"), zlib=do_zlib,\
                                               least_significant_digit=lsd)
        albedo_nir.units = "broad band albedo"
        albedo_swir = albedo_grp.createVariable("albedo_swir", "f4", ("day", "x", "y"), zlib=do_zlib,\
                                                least_significant_digit=lsd)
        albedo_swir.units = "broad band albedo"

        # albedo uncertainty
        albedo_vis_sd = albedo_grp.createVariable("albedo_vis_sd", "f4", ("day", "x", "y"), zlib=do_zlib,\
                                                  least_significant_digit=lsd)
        albedo_vis_sd.units = "albedo standard deviation"
        albedo_nir_sd = albedo_grp.createVariable("albedo_nir_sd", "f4", ("day", "x", "y"), zlib=do_zlib,\
                                                  least_significant_digit=lsd)
        albedo_nir_sd.units = "albedo standard deviation"
        albedo_swir_sd = albedo_grp.createVariable("albedo_swir_sd", "f4", ("day", "x", "y"), zlib=do_zlib,\
                                                   least_significant_digit=lsd)
        albedo_swir_sd.units = "albedo standard deviation"

        # latitude and longitude arrays
        lat = rootgrp.createVariable('lat', "f4", ("x", "y"), zlib=do_zlib)
        lat.units = "latitude"
        lon = rootgrp.createVariable('lon', "f4", ("x", "y"), zlib=do_zlib)
        lon.units = "longitude"

        # save geo information and projection
        # nchar_proj = len(proj_str)
        # nchar_geo = len(geo_str)
        rootgrp.createDimension('nchar_proj', 400)
        rootgrp.createDimension('nchar_geo', 100)
        rootgrp.createVariable('proj', 'S1', ('nchar_proj'))
        rootgrp.createVariable('geo_transform', 'S1', 'nchar_geo')

        # proj_char = nc.stringtochar(np.array([proj_str], 'S%d' % nchar_proj))
        # geo_char = nc.stringtochar(np.array([geo_str], 'S%d' % nchar_geo))
        # rootgrp.variables['proj'][:] = proj_char
        # rootgrp.variables['geo_transform'][:] = geo_char

        return rootgrp


    # **********************************************


    def make_albedo(self, x_sol, unc, nbands=7):
        """
        Calculate broad band albedo and reflectance for given BRDF kernel solution

        :param x_sol:
        :return: bhr_spectral, bhr_spectral_unc, bhr_bb, bhr_bb_unc
        """

        n_doys = x_sol.shape[0] / 3
        n_bands = x_sol.shape[1]
        bhr_spectral = np.zeros((n_doys, n_bands))
        bhr_spectral_unc = np.zeros((n_doys, n_bands))
        bhr_spectral_nbar = np.zeros((n_doys, n_bands))
        bhr_spectral_nbar_unc = np.zeros((n_doys, n_bands))
        bhr_bb = np.zeros((n_doys, 3))
        bhr_bb_unc = np.zeros((n_doys, 3))
        to_vis = np.array([0.3265, 0., 0.4364, 0.2366, 0, 0, 0])
        a_to_vis = -0.0019
        to_nir = np.array([0., 0.5447, 0, 0, 0.1363, 0.0469, 0.2536])
        a_to_nir = -0.0068
        to_sw = np.array([0.3973, 0.2382, 0.3489, -0.2655, 0.1604, -0.0138, 0.0682])
        a_to_sw = 0.0036
        for band in xrange(n_bands):
            u1 = np.sqrt(unc[band].diagonal()[:n_doys])
            u2 = np.sqrt(unc[band].diagonal()[n_doys:(n_doys * 2)])
            u3 = np.sqrt(unc[band].diagonal()[2 * n_doys:])

            bhr_spectral[:, band] = (x_sol[:(n_doys), band] +
                                     0.189184 * x_sol[(n_doys):(2 * n_doys), band] +
                                     1.377622 * x_sol[(2 * n_doys):, band])
            bhr_spectral_unc[:, band] = (u1 + 0.189184 * u2 + 1.377622 * u3)

        bhr_bb[:, 0] = np.sum(bhr_spectral * to_vis, axis=1) + a_to_vis
        bhr_bb[:, 1] = np.sum(bhr_spectral * to_nir, axis=1) + a_to_nir
        bhr_bb[:, 2] = np.sum(bhr_spectral * to_sw, axis=1) + a_to_sw
        bhr_bb_unc[:, 0] = np.sum(bhr_spectral_unc * to_vis, axis=1) + a_to_vis
        bhr_bb_unc[:, 1] = np.sum(bhr_spectral_unc * to_nir, axis=1) + a_to_nir
        bhr_bb_unc[:, 2] = np.sum(bhr_spectral_unc * to_sw, axis=1) + a_to_sw

        kr = Kernels(0, 20, 0, \
                     LiType='Sparse', doIntegrals=False, \
                     normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                     RossType='Thick')

        # n_doys = x_sol.shape[0] / 3
        for band in xrange(self.nbands):
            bhr_spectral_nbar[:, band] = (x_sol[:(n_doys), band] +
                                          kr.Ross[0] * x_sol[(n_doys):(2 * n_doys), band] +
                                          kr.Li[0] * x_sol[(2 * n_doys):, band])

            u1 = np.sqrt(unc[band].diagonal()[:n_doys])
            u2 = np.sqrt(unc[band].diagonal()[n_doys:(n_doys * 2)])
            u3 = np.sqrt(unc[band].diagonal()[2 * n_doys:])
            bhr_spectral_nbar_unc[:, band] = (u1 + kr.Ross[0] * u2 + kr.Li[0] * u3)

        return bhr_spectral[self.min_doy:self.max_doy], bhr_spectral_unc[self.min_doy:self.max_doy],\
               bhr_bb[self.min_doy:self.max_doy], bhr_bb_unc[self.min_doy:self.max_doy],\
               bhr_spectral_nbar[self.min_doy:self.max_doy], bhr_spectral_nbar_unc[self.min_doy:self.max_doy]


    # *******************************************************************************************


    def do_window(self, netcdf, output, year, cols, rows, cx, cy, n_days, doys, doy_range):
        """
        Solve problem for each pixel inside moving window.

        Parameters
        --------------------
        ff: dictionary
            hdf file dictionary where we want to save result
        fnames: array
            an array of file names
        cx: int
            x pixel coordinate of window centre
        cy: int
            y pixel coordinate of window centre
        """

        win_x = np.arange(cols)
        win_y = np.arange(rows)

        # ------------------------------------------------
        # Get optimal gamma
        # ------------------------------------------------
        # get random coordinates for estimation of gamma using only 'good' pixels

        # gx_coord = (np.random.rand(cols / 3 + 1) * cols).astype(int)
        # gy_coord = (np.random.rand(rows / 3 + 1) * rows).astype(int)

        ind = np.where(output['qa_passer'] == True)
        if ind[0].shape[0] == 0:
            print 'no good pixels for this window'
            return -1

        gx_coord = []
        gy_coord = []
        # must be better and siml
        for ii in xrange((cols / 3 + 1)):
            rand = int(np.random.rand() * ind[0].shape[0])
            gx_coord = np.append(gx_coord, np.array(ind)[:, rand][1:][1]).astype(int)
            gy_coord = np.append(gy_coord, np.array(ind)[:, rand][1:][0]).astype(int)

        # Estimate gamma for several random coordinates
        gamma_rand = []
        for ii in xrange(gx_coord.shape[0]):
            # print 'estimating gamma N ', ii, gx_coord[ii], gy_coord[ii]
            if np.array(output['doys'][output['qa_passer'][:, gy_coord[ii], gx_coord[ii]]]).shape[0] > 10:
                rho, doys, years, sza, vza, raa = self.get_rho(output, gx_coord[ii], gy_coord[ii])

                #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # doy_range = np.unique(doys)
                #!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                passer = (doys - 1).astype(np.int)
                # gamma = 10000000000
                gamma_arr = []
                opt_rmse_arr = []
                # We can evaluate gamma several time and choose the best one

                # print 'ready to guess gamma with rho', rho.shape
                timeBefore = time.clock()
                for i in xrange(3):
                    gamma0, opt_rmse0 = self.guess_gamma_opt(rho, sza, vza, raa, self.bu, doys,\
                                                        gamma_min=self.gamma_min, gamma_max=self.gamma_max,\
                                                        n_samples=self.n_samples, doy_range=doy_range)
                    gamma_arr = np.append(gamma_arr, gamma0)
                    opt_rmse_arr = np.append(opt_rmse_arr, opt_rmse0)
                timeAfter = time.clock()
                elapsed_time = timeAfter - timeBefore
                # print 'guessing gamma elapsed time (s): ', elapsed_time
                
                gamma_rand = np.append(gamma_rand, gamma_arr[opt_rmse_arr.argmin()])
                gamma_rand, g_ind = np.unique(gamma_rand, return_index=True)
            # else:
            #     print 'QA for gamma has not passed'
        # print 'opt gamma:', gamma_rand
        if gamma_rand == []:
            print 'no gamma'
            return -1

        # ----------------------------------------------------------------
        # Solve the edge-preserving task for each gamma from previous step
        # ----------------------------------------------------------------
        gamma_opt_arr = np.zeros((np.max(doy_range) * 3, 7, gamma_rand.shape[0]))
        for ii, gamma in enumerate(gamma_rand):
            gamma_arr = np.ones((np.max(doy_range) * 3, 7)) * gamma
            w = np.ones((np.max(doy_range) * 3, 7))
            for edge_i in xrange(0, 100):
                gamma_arr = gamma_arr * w
                x_sol, rmse, y_fwd, y_orig, unc, cov_mat = self.solve_regularised_problem_edge(gamma_arr, self.bu, doys, \
                                                                                          vza, sza, raa, rho, \
                                                                                          do_unc=True, \
                                                                                          doy_range=doy_range)
                # if ii == 0:
                #     plt.plot(x_sol[:,1])
                #     plt.show()
                I = np.diag(np.ones(np.max(doy_range) * 3))
                D = (I - np.roll(I, -1)).T
                for band in xrange(7):
                    try:
                        w[:, band] = 1 / (np.sqrt(1 + (x_sol[:, band].dot(D) ** 2)))
                    except:
                        print 'error in weights estimation'
                        print 'x_sol: ', x_sol.shape
                        print 'w: ', w.shape
                        print 'D: '. D.shape
                        print 'doy_range: ', doy_range, np.max(doy_range)
                        continue
            # plt.plot(x_sol[:, 1])
            # plt.show()
            # plt.plot(gamma_arr[:, 1])
            # plt.show()
            gamma_opt_arr[:, :, ii] = gamma_arr



        # ----------------------------------------------------
        # Do inversion for every pixel
        # ----------------------------------------------------
        for xi, xx in enumerate(win_x):
            for yi, yy in enumerate(win_y):

                # print 'pixel:', xx, yy
                timeBefore0 = time.clock()

                # Get reflectance and then estimate optimal gamma for a given pixel.
                # We do it by estimating difference of temporal profiles of
                # gamma-pixel and a current pixel

                # Get number of DOYs
                n_doys = np.array(output['doys']).shape[0]
                # and number of DOYs which have passed QA
                n_doys_good = np.array(output['doys'][output['qa_passer'][:, yi, xi]]).shape[0]
                # print 'n_doys=', n_doys
                # print 'n_doys_good=', n_doys_good

                # if number of 'good' days in pixel more than N go further
                if n_doys_good > 10:
                    # print 'get rho for xi=%d, yi=%d'%(xi, yi)
                    rho, doys, years, sza, vza, raa = self.get_rho(output, xi, yi)
                else:
                    # print 'pixel (%d,%d) is not good'%(xx, yy)
                    continue

                diff = []
                # print 'estimate gamma for current pixel'
                for ii, gamma in enumerate(gamma_rand):
                    if np.array(output['doys'][output['qa_passer'][:, gy_coord[g_ind][ii], gx_coord[g_ind][ii]]]).shape[0] > 10:
                        # print 'get rho for gamma=%.5f' % gamma
                        grho, gdoys, gyears, gsza, gvza, graa = self.get_rho(output, gx_coord[g_ind][ii], gy_coord[g_ind][ii])

                        # Make arrays of reflectance comparable
                        rho0 =  np.zeros((doy_range.shape[0], 7))
                        grho0 = np.zeros((doy_range.shape[0], 7))
                        for jj, dd in enumerate(doy_range):
                            if dd in doys:
                                # We can have more than one observation for one day.
                                # In this case observations are avereged
                                rho0[jj-1, :] = np.mean(rho[doys==dd, :], axis=0)
                            if dd in gdoys:
                                grho0[jj-1, :] = np.mean(grho[gdoys==dd, :], axis=0)

                        diff = np.append(diff, np.sum((rho0 - grho0)**2))
                    # else:
                        # print 'something wrong', gy_coord[g_ind][ii], gx_coord[g_ind][ii]
                        # print np.array(output['doys'][output['qa_passer'][:, gy_coord[g_ind][ii], gx_coord[g_ind][ii]]]).shape[0]
                if diff == []:
                    # if diff is empty break the loop and go to next pixel
                    # print 'diff is empty'
                    continue
                gamma = gamma_rand[np.argmin(diff)]

                # print 'ready to solve problem with rho ', rho.shape
                timeBefore = time.clock()
                # try:
                x_sol, rmse, y_fwd, y_orig = self.solve_regularised_problem(gamma, self.bu, doys, vza, sza, raa, rho,
                                                                                do_unc=False, doy_range=doy_range)
                # except:
                    # print 'regularised problem for pixel %d %d has not solved' % (xx, yy)
                    # continue
                timeAfter = time.clock()


                elapsed_time = timeAfter - timeBefore
                # print '1st solving elapsed time (s): ', elapsed_time
                                                             
                # print 'filter outliers'
                rpasar, pasar = self.outlier_filter(doys, rho, y_fwd, self.bu,\
                                                    'fig/outlier_%d_%d.png'%(xi, yi), threshold=10., do_plots=False)
                doys = doys[pasar]
                rho = rho[pasar, :]
                sza = sza[pasar]
                vza = vza[pasar]
                raa = raa[pasar]

                # print 'solve problem with outliers filtered with rho, ', rho.shape
                timeBefore = time.clock()




                #***********************************
                # Start the edge-preserving thing with gamma-arrays
                #***********************************
                x_sol, rmse, y_fwd, y_orig, unc, cov_mat = self.solve_regularised_problem_edge(gamma_opt_arr[:, :, np.argmin(diff)], self.bu, doys,\
                                                                                          vza, sza, raa, rho,\
                                                                                          do_unc=True,\
                                                                                          doy_range=doy_range)

					
					
                if x_sol[0,0] == -1:
                    # print '**********************************************'
                    # print 'go to the next pixel'
                    # print '**********************************************'
                    continue
                # except:
                #     print 'regularised problem for pixel %d %d has not solved' % (xx, yy)
                #     continue

                timeAfter = time.clock()
                elapsed_time = timeAfter - timeBefore

                bhr_spectral, bhr_spectral_unc, \
                bhr_bb, bhr_bb_unc, \
                bhr_spectral_nbar, bhr_spectral_nbar_unc = self.make_albedo(x_sol, unc)

                # self.plot_results(x_sol, y_fwd, bhr_spectral_nbar, bhr_spectral_nbar_unc, bhr_spectral, rho, doy_range,
                #                   doys)

                # print '2nd solving elapsed time (s): ', elapsed_time
                # print 'get an answer and ready to save result in netcdf'


                # Copy results of forward modeling to a hdf file
                # We go around of the central pixel of window
                # netcdf.groups['reflectance'].variables['y_fwd'][cy + yy, cx + xx, :, :] = y_fwd[self.min_doy:self.max_doy, :]

                ind = np.logical_and(doy_range >= self.min_doy, doy_range < self.max_doy)

                # netcdf.groups['reflectance'].variables['y_fwd'][yy, xx, :, :] = y_fwd[self.min_doy:self.max_doy, :]

                # print 'netcdf.groups[reflectance].variables[y_fwd][:, :, yy, xx]',\
                #       netcdf.groups['reflectance'].variables['y_fwd'][:, :, yy, xx].shape
                # print 'y_fwd[self.min_doy:self.max_doy, :]', y_fwd[self.min_doy:self.max_doy, :].shape

                netcdf.groups['reflectance'].variables['y_fwd'][:, :, yy, xx] = y_fwd[self.min_doy:self.max_doy, :]

                netcdf.groups['reflectance'].variables['y_orig'][:, :, yy, xx] = y_orig[self.min_doy:self.max_doy, :]

                # for bb in range(1, 8):
                #     np.savez('covariance_%s/unc_%d-%d_%d_%d_band_%d_m' % (site, year_range[0], year_range[-1], xx, yy, bb), \
                #              data=cov_mat[bb - 1].data, indices=cov_mat[bb - 1].indices, \
                #              indptr=cov_mat[bb - 1].indptr, shape=cov_mat[bb - 1].shape)



                # Write to array
                for band in xrange(7):
                    the_band = "band_%d" % (band + 1)

                    netcdf.groups["reflectance"].variables["refl_b%d" % (band + 1)][:, yy, xx] = bhr_spectral_nbar[:, band]
                    netcdf.groups["reflectance"].variables["refl_b%d_sd" % (band + 1)][:, yy, xx] = bhr_spectral_nbar_unc[:, band]

                    the_band = "band_%d_nbar" % (band + 1)

                bb_name = ['vis', 'nir', 'swir']
                for band in xrange(3):
                    the_band = "bb_%s" % bb_name[band]

                    netcdf.groups["albedo"].variables["albedo_%s" % bb_name[band]][:, yy, xx] = bhr_spectral_nbar[:, band]
                    netcdf.groups["albedo"].variables["albedo_%s_sd" % bb_name[band]][:, yy, xx] = bhr_spectral_nbar_unc[:, band]

                # print 'output[doys_real]:', output['doys_real']
                # print 'output[doys_real]:', output['doys_real'][self.min_doy:self.max_doy]
                # print 'self.min_doy, self.max_doy', self.min_doy, self.max_doy

                date = np.array([julday.doy2date(year, d) for d in output['doys_real'][self.min_doy:self.max_doy]])
                # date = np.array([julday.doy2date(year, d) for d in output['doys_real'][:]])
                date_str = np.array([d.strftime('%Y.%m.%d') for d in date])
                date_jul = np.array([julday.date2jul(d) for d in date])

                # Range of dates for one year \\\\\
                nd, dr = self.get_n_days(year)

                # date_str = np.array([(datetime.date(year, 1, 1) +\
                #                            datetime.timedelta(days=i)).strftime('%Y.%m.%d') \
                #                            for i in xrange(0, nd)])

                # print 'date_str:'
                # print date_str.shape
                # print date_str
                # print netcdf.variables['date_str'].shape
                # print date_jul

                date_jul = []
                for dd in xrange(date_str.shape[0]):
                    netcdf.variables['date_str'][dd,:] = [ss for ss in date_str[dd]]
                    tmp = np.array(date_str[dd].split(".")).astype(int)
                    # date_jul = np.append(date_jul, np.round(julday.date2jul(dt(tmp[0], tmp[1], tmp[2]))).astype(int))
                    date_jul = np.append(date_jul, julday.date2jul(dt(tmp[0], tmp[1], tmp[2])))

                netcdf.variables['julday'][:] = date_jul
                
                #netcdf.variables['lat'][cy + yy, cx + xx] = output['lat'][yi, xi]
                #netcdf.variables['lon'][cy + yy, cx + xx] = output['lon'][yi, xi]

                timeAfter0 = time.clock()
                elapsed_time = timeAfter0 - timeBefore0
                # print 'total time (s) for pixel: ', elapsed_time
        # ------------------------------------------------------------------------------------------

        # self.plot_results(x_sol, y_fwd, bhr_spectral_nbar, bhr_spectral_nbar_unc, bhr_spectral, rho, doy_range, doys)



    # *******************************************************************************************



    def do_job(self, nc_file, loc, pattern, year, tile, ul_px, ul_py, lr_px, lr_py, cols, rows, \
               step=1, do_adj_years=True, do_vrt = False):
        """
        Do the main job. Call the routine of hdf5 creation, estimate gamma,
        solve regularized problem and save a hdf5 file.

        Parameters
        -----------
        hdf_dir: string
        loc: array
        site: string
        year_range: array
        tile: string
        ul_x: float
        ul_y: float
        lr_x: float
        lr_y: float
        cols: int
            window size x
        rows: int
            window size y
        do_vrt: bool
        """

        timeBefore = time.clock()

        # Size of an image must be multiple of window size
        N = (lr_px-ul_px) / int(cols) * cols
        M = (lr_py-ul_py) / int(rows) * rows

        # print 'ul_px, ul_py, lr_px, lr_py  = ', ul_px, ul_py, lr_px, lr_py
        # print 'Image size: %dx%d' % (N, M)
        # print 'Window size: %dx%d' % (cols, rows)
        if (lr_px-ul_px) < cols:
            print 'x size of the image must be >= than x size of window'
            return -1
        if (lr_py-ul_py) < rows:
            print 'y size of the image must be >= than y size of window'
            return -1

        # Get dates
        n_days = 0
        nd, dr = self.get_n_days(year)
        n_days = n_days + nd

        # Central pixel incide a window in pixel coordinates of this window.
        # I.e. if window is size is 10x10 then coordinates will be 0..9, 0..9 and
        # cenral pixel is 5.
        cx = int(cols/2.)
        for px in xrange(ul_px, ul_px+N, cols):
            cy = int(rows/2.)
            for py in xrange(ul_py, ul_py+M, rows):
                n_days, dr = self.get_n_days(year - 1)

                # print 'Read time series\n'
                output = OrderedDict()
                output_2 = self.get_time_series(tile, year, loc, px, py, cols, rows, step=step, save_dir='data/')
                if do_adj_years:
                    output_1 = self.get_time_series(tile, year-1, loc, px, py, cols, rows, step=step, save_dir='data/')
                    output_3 = self.get_time_series(tile, year+1, loc, px, py, cols, rows, step=step, save_dir='data/')
                    if output_2 == -1:
                        continue
                    if output_1 == -1:
                        output_1 = output_2.copy()
                        output_1['years'] = output_1['years'] - 1
                    if output_3 == -1:
                        output_3 = output_2.copy()
                        output_3['years'] = output_3['years'] + 1
                    keys = output_2.keys()
                    for key in keys:
                        output[key] = np.append(output_1[key], output_2[key], axis=0)
                        output[key] = np.append(output[key], output_3[key], axis=0)
                    nd1, dr = self.get_n_days(year - 1)
                    nd2, dr = self.get_n_days(year)
                    nd3, dr = self.get_n_days(year + 1)
                    n_days = nd1 + nd2 + nd3
                    self.min_doy = nd1
                    self.max_doy = nd1 + nd2

                else:
                    n_days, dr = self.get_n_days(year)
                    self.min_doy = 0
                    self.max_doy = n_days
                    output = output_2

                doys = output['doys'].copy()
                years = output['years'].copy()

                doys[years == years[0] + 1] = doys[years == years[0] + 1] + np.max(output['doys'])
                doys[years == years[0] + 2] = doys[years == years[0] + 2] + np.max(output['doys']) * 2
                self.max_doy = np.max(doys[years == year])
                self.min_doy = np.min(doys[years == year])
                doy_range = np.unique(doys)

                # Create a netcdf file to store results
                netcdf = self.create_netcdf(nc_file, N, M, doy_range[self.min_doy : self.max_doy])

                # print 'do window for px=%d, py=%d\n' % (px, py)
                # Move central pixels and thus do processing inside a new window
                self.do_window(netcdf, output, year, cols, rows, cx, cy, n_days, doys, doy_range)

                cy = cy + rows

            cx = cx + cols

            timeAfter = time.clock()
            elapsed_time = timeAfter - timeBefore
            print 'total elapsed time for job (s): ', elapsed_time

        return 0


# *******************************************************************************************


# *******************************************************************************************
# MAIN
# *******************************************************************************************

if __name__ == "__main__":
    """
    An example of call:
    python brdf_filter_img_12.py -s 'viterbo' -t 'h18v04' -y 42.38041111 -x 12.02656111 -g 2014 
    """

    timeBefore = time.clock()

    parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(),
                                   usage=globals()['__doc__'])
    #python %s --ncd_out %s --root %s --pattern %s --ulx %d --uly %d --cols %d --rows %d --tile %s --year %s


    parser.add_option('--ncd_out', action="store", dest="ncd_out",
                      type=str, help="output netCDF file")

    parser.add_option('--root', action="store", dest="root",
                      type=str, help="dir where source files are located")

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

    parser.add_option('--tile', action="store", dest="tile",
                      type=str, help="MODIS tile")

    parser.add_option('--year', action="store", dest="year",
                      type=str, help="year")

    parser.add_option('--step', action="store", dest="step",
                      type=str, help="step in number of days")


    (options, args) = parser.parse_args()
    # site = options.site
    tile = options.tile

    print 'tile:',  options.tile
    print 'ul x: ', options.ulx
    print 'ul y: ', options.uly
    # year = options.pattern[8:12]
    print 'year:', options.year
    lrx = int(options.ulx) + int(options.cols)
    lry = int(options.uly) + int(options.rows)
    print 'lr x', lrx
    print 'lr y', lry
    # print 'year: ', year

    brf = opt_img()

    brf.do_job(options.ncd_out, options.root, options.pattern, int(options.year),\
               options.tile,\
               int(options.ulx), int(options.uly), lrx, lry,\
               int(options.cols), int(options.rows), int(options.step))

    print 'kernel_img is done!!!'
    timeAfter = time.clock()
    print 'it lasted for (h): ', (timeAfter - timeBefore)/(60.*60.)
