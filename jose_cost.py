# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 11:55:51 2016

@author: Jose Gomez-Dans, Maxim Chernetskiy

Functions for minimization of time series
with identity observation operator
"""

import numpy as np
#import datetime as dtay
#import scipy.optimize
#import scipy.special as ssp


def der_cost_function ( x, y, x_prior, sigma_obs, sigma_prior, gamma, obs_mask ):
    """This function calculates the partial derivatives of the cost function"""
    I = np.identity(x.shape[0])
    D1 = np.matrix(I - np.roll(I,1))
    
    # Second order
    D1 = D1*D1.T
    
    der_j_obs = x*0
    der_j_obs[obs_mask] = (x[obs_mask] - y)/sigma_obs**2
    der_j_prior = (x - x_prior ) /sigma_prior**2
    der_j_model = np.array(gamma*np.dot((D1).T, D1*np.matrix(x).T)).squeeze()
    return der_j_obs + der_j_prior + der_j_model



def cost_function ( x, y, x_prior, sigma_obs, sigma_prior, gamma, obs_mask ):
    """
    Variational cost function with a weak constraint

    Parameters
    -----------
    x: The current estimate of the state
    y: The observations vector
    x_prior: The prior estimate of the state
    sigma_obs: Observational uncertainty (standard deviation)
    sigma_prior: Prior uncertainty (standard deviation)
    gamma: Inverse of the model term variance
    """
    I = np.identity(x.shape[0])
    D1 = np.matrix(I - np.roll(I,1))
    
    # Second order
    D1 = D1*D1.T
    
    xa = np.matrix ( x )
    j_obs = 0.5*np.sum((x[obs_mask]-y)**2/sigma_obs**2)
    j_prior = 0.5*np.sum((x-x_prior)**2/sigma_prior**2)
    j_model = 0.5*gamma*np.dot((D1*(xa.T)).T, D1*xa.T)
    return j_obs + j_prior + j_model



def posterior_covariance (x, sigma_prior, sigma_obs, gamma, obs_mask=None):
    I = np.identity(x.shape[0])
    D1 = np.matrix(I - np.roll(I,1))
    
    hess_obs = np.zeros( x.shape[0] )
    # hess_obs[obs_mask] = (1./sigma_obs[obs_mask]**2)
    hess_obs = (1. / sigma_obs ** 2)
    hess_obs = np.diag( hess_obs )
    hess_prior = (1./sigma_prior**2)*np.diag(np.ones_like(x))

    hess_model = gamma*np.dot ( D1,np.eye(x.shape[0])).dot( D1.T)
    hess = hess_obs + hess_prior + hess_model
    try:
        post_covar = np.linalg.inv ( hess )
    except:
        print 'singular matrix?'
        print 'x, sigma_prior, sigma_obs, gamma: ', x.shape, sigma_prior.shape, sigma_obs.shape, gamma.shape
        post_covar = np.zeros(hess.shape)
    return np.sqrt(np.diagonal(post_covar))



# *************************************************************************************

def cost_function_t(x, y, x_prior, sigma_obs, sigma_prior, gamma, obs_mask, t, t_obs):
        """
        Variational cost function with a weak constraint.
        Allows several observations per time step

        Parameters
        -----------
        x: The current estimate of the state
        y: The vector of observations
        x_prior: The prior estimate of the state
        sigma_obs: Observational uncertainty (standard deviation)
        sigma_prior: Prior uncertainty (standard deviation)
        gamma: Inverse of the model term variance
        t: temporal grid
        t_obs: days of observations
        """

        x = np.squeeze(np.array(x))

        I = np.identity(x.shape[0])

        # define differential operator
        D1 = np.matrix(I - np.roll(I, 1))

        # Second order
        D1 = D1 * D1.T

        xa = np.matrix(x)
        # j_obs = 0.5 * np.sum((x[obs_mask] - y) ** 2 / sigma_obs ** 2)
        j_obs = 0
        for i in xrange(t.shape[0]-1):
            # Find observations within temporal step
            ind = np.logical_and(t_obs >= t[i], t_obs < t[i+1])
            # Calculate cost as sum of differences between observations and
            # current estimate of the state
            for j in xrange(y[ind].shape[0]):
                if y[ind][j] != 0:
                    j_obs = j_obs + (x[i] - y[ind][j]) ** 2 / sigma_obs[ind][j] ** 2

        # observational term
        j_obs = 0.5 * j_obs

        # prior term
        j_prior = 0.5 * np.sum((x - x_prior) ** 2 / sigma_prior ** 2)

        # dynamic model term
        j_model = np.array(0.5 * gamma * np.dot((D1 * (xa.T)).T, D1 * xa.T)).squeeze()

        # return sum of three terms
        return j_obs + j_prior + j_model



    #**************************************************************************



def der_cost_function_t(x, y, x_prior, sigma_obs, sigma_prior, gamma, obs_mask, t, t_obs):
        """This function calculates the partial derivatives of the cost function"""

        x = np.squeeze(np.array(x))

        I = np.identity(x.shape[0])
        D1 = np.matrix(I - np.roll(I, 1))

        # Second order
        D1 = D1 * D1.T

        der_j_obs = x * 0
        # der_j_obs = 0

        # der_j_obs[obs_mask] = (x[obs_mask] - y) / sigma_obs ** 2

        for i in xrange(t.shape[0]-1):
            ind = np.logical_and(t_obs >= t[i], t_obs < t[i+1])
            for j in xrange(y[ind].shape[0]):
                if y[ind][j] != 0:
                    # pdb.set_trace()
                    # print 'x, y', x.shape, y.shape, sigma_obs.shape
                    der_j_obs[i] = (x[i] - y[ind][j]) / sigma_obs[ind][j] ** 2

                    #der_j_obs[obs_mask] = (x[obs_mask] - y)/sigma_obs**2

                    # der_j_obs = der_j_obs + (x[i] - y[ind][j]) ** 2 / sigma_obs[ind][j] ** 2

                    # pdb.set_trace()
        # j_obs = 0.5 * j_obs


        der_j_prior = (x - x_prior) / sigma_prior ** 2
        # der_j_prior = np.sum((x - x_prior) / sigma_prior ** 2)
        der_j_model = np.array(gamma * np.dot((D1).T, D1 * np.matrix(x).T)).squeeze()

        return der_j_obs + der_j_prior + der_j_model



# ******************************************************************************