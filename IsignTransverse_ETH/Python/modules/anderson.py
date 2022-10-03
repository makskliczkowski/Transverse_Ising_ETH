import numpy as np
import pandas as pd
import os
from os import sep as kPSep
import config as cf
import importlib
importlib.reload(cf)
user_settings = getattr(cf.plot_settings, 'settings')

def info(L, J, W):
    return "_L=%d,J=%.2f,J0=%.2f,W=%.2f.dat"%(L, J, 0.0, W)

#-------------------------------------------------------- STATISTICS
#-------------------------------------------------------------------

def load_sff(L, W, dim = 3, settings = None):
    """
    Load spectral data along with statistical measures

    Parameters:
    -----------------
        L: int
            system size

        W: float
            disorder strength

        dim: int
            dimensionality of lattice (1-3D)
        
        settings: dict
            dictionairy for plot settings as in config.py
    """
    if settings is None:    settings = user_settings

    dir = f"..{kPSep}results{kPSep}" + f"ANDERSON{kPSep}%dD{kPSep}PBC{kPSep}SpectralFormFactor{kPSep}"%dim
    if settings['smoothed'] == 1: dir = dir + f"smoothed{kPSep}"

    filename = dir + info(L, 1.0, W)
    
    if os.path.exists(filename):
        data = pd.read_table(filename, sep="\t", header=None)
        xdata = np.array(data[0])
        ydata = np.array(data[1])
        tH = data[2][0]
        tau = data[3][0]
        gap_ratio = data[4][0]
        return True, xdata, ydata, tH, tau, gap_ratio
    else:
        return False, np.array([]), np.array([]), None, None, None


#-------------------------------------------------------- LOCALISATION
#---------------------------------------------------------------------
def load_loc(L = 4000, w = 0.2):
    """
    Loading localisation length for given system size and disorder strength

    Parameters:
    -----------
        L: int
            system size
        
        w: float
            disorder strength
    """
    name = "../results/ANDERSON/1D/PBC/LocalisationLength/Distribution/_L=%d,J=1.00,J0=0.00,w=%.2f.dat"%(L,w)
    loc_len = np.loadtxt(name, unpack=True)
    return np.array(loc_len)


def load_orbital(L = 4000, w = 0.2, n = 100):
    """
    Loading correlation function for given system size and disorder strength
    for orbital n

    Parameters:
    -----------
        L: int
            system size
        
        w: float
            disorder strength
    """
    name = "../results/ANDERSON/1D/PBC/CorrelationFunction/_L=%d_n=%d_w=%.2f.dat"%(L, n, w)
    orbital = np.loadtxt(name, unpack=True)
    return np.array(orbital)