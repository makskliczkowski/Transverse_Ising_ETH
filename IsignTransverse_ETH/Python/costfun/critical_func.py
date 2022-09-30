"""
This module contains predicitons of
scaling of the critical point with system size.
"""

import numpy as np
from config import hamiltonian

def _crit_free(size, *args):
    """
    Free critical point for each system size
    For Heisenberg chain the sizes are even: 12,14,16,18,..
    For Ising chain the sizes are: 11,12,13,14,...
    """
    crit = np.array(args)
    idx = int( (size - 12) / 2 if hamiltonian == 1 else size - 11 )
    return crit[idx]

def _crit_free_inv(size, *args):
    """
    Free inversed critical point for each system size if critical points appears to be very small
    For Heisenberg chain the sizes are even: 12,14,16,18,..
    For Ising chain the sizes are: 11,12,13,14,...
    """
    crit = 1. / np.array(args)
    idx = int( (size - 12) / 2 if hamiltonian == 1 else size - 11 )
    return crit[idx]

def _crit_const(size, x0):
    """
    Constant critical value
    """

    return x0


def _crit_lin(size, x0, x1):
    """
    Linear scaling of critical point with size
    """

    return x0 + size * x1

def _crit_inv(size, x0, x1, x2):
    """
    Scaling with inverse of system size
    """
    return x0 + 1.0 / float(x1 + x2 * size)

def _crit_power_law(size, x0, x1, nu):
    """
    Power law prediciton with system size
    (Only leading order taken)
    """

    return x0 + float(size)**nu * x1

def _crit_log(size, x0, x1):
    """
    Ansatz with logarithmic scaling with size
    """

    return x0 + x1 * np.log(size)


def _crit_inv_log(size, x0, x1):
    """
    Inverse of logarithmic scaling
    """
    return x0 + x1 / np.log(size)
