#import user modules
import utils.helper_functions as hfun
import config as cf
import modules.thouless_times as thouless
import modules.spectral_functions as spec_fun
import modules.sff as sff_module
import modules.adiabatics as agp
import modules.anderson as anderson
import importlib
from utils.fit_functions import *

from utils.fit_functions import *

#--- NUMERICAL LIBS
import numpy as np
import itertools
import math
import random
from cmath import nan
import h5py   


# SCIPY LIBS
import scipy.stats as statistics
from scipy.special import binom
from scipy.special import erfinv
from scipy.special import digamma
from scipy.special import polygamma
from scipy.optimize import curve_fit as fit
from scipy.signal import savgol_filter
 
# OTHER
import warnings
warnings.filterwarnings('ignore')
from joblib import Parallel, delayed
import copy
import os
from utils import exit
from os import sep as kPSep
from os.path import exists

print(cf.base_directory)


user_settings = getattr(cf.plot_settings, 'settings')

if __name__ == '__main__':

    J=1.0
    bins = 100

    sizes = np.arange(8, 42, 8)
    W_vals = hfun.regspace(0.05, 1.05, 0.05)
    num_of_points = 40
    def get_pr(num_realis, W_vals, plane_wave_basis = False):
        E = []
        pr = []
        DOS = []
        
        E_density = []
        pr_density = []
        dirr = f"..{kPSep}results{kPSep}" + f"ANDERSON{kPSep}3D{kPSep}PBC{kPSep}MultiFractality{kPSep}" + (f"ParticipationRatio_PlaneWave{kPSep}" if plane_wave_basis else f"ParticipationRatio{kPSep}")
        for L in sizes:
            dim = L**3
            E_L = []
            pr_L = []
            DOS_L = []
            
            E_density_L = []
            pr_density_L = []
            for W in W_vals:
                E_tmp = np.zeros( (dim) )
                pr_tmp = np.zeros( (dim) )

                DOS_tmp = np.zeros ( (bins) )
                counter = 0
                E_density_tmp = np.zeros((num_of_points))
                pr_density_tmp = np.zeros((num_of_points))
                for realis in range(num_realis):
                    filename = dirr + f'realisation={realis}{kPSep}' + hfun.remove_info(anderson.info(L, J, W), 'q') + '_q=2.00.hdf5'
                    if os.path.exists(filename):
                        print(filename)
                        with h5py.File(filename, "r") as file:
                            ener = np.array(file.get('energies'))[0]
                            paratata = np.array(file.get('participation_ratio'))[0]
                            hist, _ = np.histogram(ener, bins=bins, normed=True)
                            
                            DOS_tmp += hist
                            E_tmp += ener
                            pr_tmp += paratata

                            M = len(ener) // num_of_points  # number of points in average
                            bucket_num = len(ener) // M     # number of buckets
                            for kk in range(0, num_of_points):
                                pr_tmp_tmp_tmp_tmp = 0
                                e_tmp_tmp = 0
                                cunt = 0
                                end = min( (kk + 1) * M, len(ener) - 1)
                                for jj in range(kk * M, end ):
                                    pr_tmp_tmp_tmp_tmp += paratata[jj]
                                    e_tmp_tmp += ener[jj]
                                    cunt += 1
                                if cunt > 0:
                                    pr_density_tmp[kk] += pr_tmp_tmp_tmp_tmp / cunt #np.mean(gaps_tmp[ kk * M : (kk + 1) * M ])
                                    E_density_tmp[kk]  += e_tmp_tmp / cunt #np.mean(ener    [ kk * M : (kk + 1) * M ])
                            counter += 1
                print(L, W, counter)

                E_L.append(E_tmp / counter)
                pr_L.append(pr_tmp / counter)
                DOS_L.append(DOS_tmp / counter)
                pr_density_L.append(pr_density_tmp / counter)
                E_density_L.append(E_density_tmp / counter)

            E.append(np.array(E_L))
            pr.append(np.array(pr_L))
            DOS.append(np.array(DOS_L))
            pr_density.append(np.array(pr_density_L))
            E_density.append(np.array(E_density_L))
        return E, pr, DOS, pr_density, E_density


    E, pr, DOS, pr_density, E_density = get_pr(4200, W_vals)


    with open(f'ANDERSON3D_DATA/density_of_states.npy', 'wb') as file:      np.save(file, DOS, allow_pickle=True)

    with open(f'ANDERSON3D_DATA/energy_densities_real.npy', 'wb') as file:  np.save(file, E_density, allow_pickle=True)
    with open(f'ANDERSON3D_DATA/ipr_densities_real.npy', 'wb') as file:     np.save(file, pr_density, allow_pickle=True)
    with open(f'ANDERSON3D_DATA/mean_energies_real.npy', 'wb') as file:     np.save(file, E)
    with open(f'ANDERSON3D_DATA/ipr_real.npy', 'wb') as file:               np.save(file, pr)
        
    E_kspace, pr_kspace, DOS_kspace, pr_density_kspace, E_density_kspace = get_pr(100, W_vals, plane_wave_basis=True)

    with open(f'ANDERSON3D_DATA/density_of_states_kspace.npy', 'wb') as file:      np.save(file, DOS_kspace, allow_pickle=True)

    with open(f'ANDERSON3D_DATA/energy_densities_kspace.npy', 'wb') as file:  np.save(file, E_density_kspace, allow_pickle=True)
    with open(f'ANDERSON3D_DATA/ipr_densities_kspace.npy', 'wb') as file:     np.save(file, pr_density_kspace, allow_pickle=True)
    with open(f'ANDERSON3D_DATA/mean_energies_kspace.npy', 'wb') as file:     np.save(file, E_kspace)
    with open(f'ANDERSON3D_DATA/ipr_kspace.npy', 'wb') as file:               np.save(file, pr_kspace)