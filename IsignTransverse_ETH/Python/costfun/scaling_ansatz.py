"""
Here the differeny types of scalings ansatzes are presented
to find the proper critical point
The names should be buolt as _rescale_*, where
the * should be replaced with a shoirt name characterizing the scaling
ansatz, like KT for Kosterlitz-Thouless
"""

import numpy as np
from config import hamiltonian
from scipy.special import binom

def _rescale_spacing(x, L, crit_fun, *args):
    """Regular ansatz with power-law on L"""
    return np.sign(x - crit_fun(L, *args)) * abs(x - crit_fun(L, *args))

def _rescale_classic(x, L, crit_fun, nu, *args):
    """Regular ansatz with power-law on L"""
    return (x - crit_fun(L, *args)) * L**(nu)


def _rescale_KT(x, L, crit_fun, nu, *args):
    """ Kosterlitz Thouless type """
    return  np.sign(x - crit_fun(L, *args)) * L / np.exp(nu / np.sqrt(abs(x - crit_fun(L, *args) )) )


def _rescale_FGR(x, L, crit_fun, nu, *args):
    """Fermi Golden Rule"""
    dim = binom(L, L/2) if hamiltonian == 1 else 2**L
    return np.sign(x - crit_fun(L, *args)) * abs(x - crit_fun(L, *args)) * dim**(1. / nu)


def _rescale_RG(x, L, crit_fun, nu, *args):
    """ RG-type arguments """
    return np.sign(x - crit_fun(L, *args)) * abs(x - crit_fun(L, *args))**nu * L

def _rescale_exp(x, L, crit_fun, nu, *args):
    """Exponential scaling of parameter"""
    return np.sign(x - crit_fun(L, *args)) * np.exp(nu*abs(x - crit_fun(L, *args))) * np.exp(np.log(2) / 2 * abs(L))