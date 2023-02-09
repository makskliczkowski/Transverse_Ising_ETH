import numpy as np
import pandas as pd
import os
from os import sep as kPSep
from modules.sff import GOE
import config as cf
import utils.helper_functions as hfun
import importlib
importlib.reload(cf)
importlib.reload(hfun)
from scipy.signal import savgol_filter

user_settings = getattr(cf.plot_settings, 'settings')

def info(L, J, W):
    return "_L=%d,J=%.2f,J0=%.2f,W=%.2f.dat"%(L, J, 0.0, W)

#-------------------------------------------------------- STATISTICS
#-------------------------------------------------------------------

#-------------- SPECTRAL FORM FACTOR
def load_sff(L, W, dim = 3, settings = None):
    """
    Load spectral form factor and other data along with statistical measures

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
    #if settings['smoothed'] == 1: dir = dir + f"smoothed{kPSep}"

    filename = dir + info(L, 1.0, W)
    
    epsilon = 2e-1
    if os.path.exists(filename):
        data = pd.read_table(filename, sep="\t", header=None)
        times = np.array(data[0])
        sff = np.array(data[1])
        if settings['smoothed'] == 1:
            sff = hfun.remove_fluctuations(sff, 50)
            #sff = savgol_filter(sff, window_length=51, polyorder=5, mode="mirror")
        tH = data[2][0]
        tau = data[3][0]
        sff_dev = np.abs(np.log10(sff / GOE(times)))
        for i, K in reversed(list(enumerate(sff_dev))):
            if K > epsilon and times[i] < (3):
                tau = times[i-1]
                break
        gap_ratio = data[4][0]
        return True, times, sff, tH, tau, gap_ratio
    else:
        #print(filename)
        return False, np.array([]), np.array([]), None, None, None


#-------------- DENSITY OF STATES
def load_DOS(L, W, dim = 3, settings = None, unfolded = True, use_500_states = True):
    """
    Load density of states either full spectrum or 500 states at middle of spectrum

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

        unfolded: bool
            DOS for unfolded spectrum

        use_500_states: bool
            pick half of spectrum (False) or 500 at middle of spectrum (True)
    """
    if settings is None:    settings = user_settings

    dir = f"..{kPSep}results{kPSep}" + f"ANDERSON{kPSep}%dD{kPSep}PBC{kPSep}DensityOfStates{kPSep}"%dim

    prefix = "500_states" if use_500_states else ""
    if unfolded: prefix += "unfolded"
    
    filename = dir + prefix + info(L, 1.0, W)
    
    if os.path.exists(filename):
        data = pd.read_table(filename, sep="\t", header=None)
        xdata = np.array(data[0])
        ydata = np.array(data[1])
        return True, xdata, ydata
    else:
        return False, np.array([]), np.array([])


#-------------- LEVEL SPACING STATISTICS
def load_level_spacing_dist(L, W, dim = 3, settings = None, log_data = False, unfolded = True, use_500_states = True):
    """
    Load density of states either full spectrum or 500 states at middle of spectrum

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

        log_data: bool
            boolean value whether to use load x axis spaced logarithmically

        unfolded: bool
            DOS for unfolded spectrum

        use_500_states: bool
            pick half of spectrum (False) or 500 at middle of spectrum (True)
    """
    if settings is None:    settings = user_settings

    dir = f"..{kPSep}results{kPSep}" + f"ANDERSON{kPSep}%dD{kPSep}PBC{kPSep}LevelSpacingDistribution{kPSep}"%dim

    prefix = "_500_states" if use_500_states else ""
    if unfolded: prefix += "unfolded"
    if log_data: prefix += "_log"

    filename = dir + prefix + info(L, 1.0, W)
    #print(filename)
    if os.path.exists(filename):
        data = pd.read_table(filename, sep="\t", header=None)
        xdata = np.array(data[0])
        ydata = np.array(data[1])
        wH_typ_unfolded = data[2][0]
        wH = data[3][0]
        wH_typ = data[4][0]
        return True, xdata, ydata, wH_typ_unfolded, wH, wH_typ
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