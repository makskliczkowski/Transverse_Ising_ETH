import importlib
from os import sep as kPSep
import utils.plot_settings as ps
importlib.reload(ps)
#---------------------------------------------------- MODEL PARAMETERS
model = 2           # chooses model: 0-disorder / 1-symmetries / 2-local perturbation
hamiltonian = 1     # which hamiltonian?: 0-Ising / 1-Heisenberg
BC = 0              # boundaary condition: 0 - OBC / 1 - PBC

L = 18                          # system size
J = 1.00                        # spin exchange (Ising-like)
g = 0.55                       # trasnverse magnetic field (z-axis)
h = 0.0                        # longitudal magnetic field (x-axis)
#---- DISORDER PARAMETERS
w = 0.5                        # disorder on longitudonal field ( h_i \in [h-w, h+w] )
J0 = 0.0                        # disorder on spin exchange ( J_i \in [J-J0, J+J0] )
g0 = 0.0                        # disorder on longitudonal field ( h_i \in [h-w, h+w] )
#---- SYMETRY PARAMETERS
k_sym = 0                       # translational symmetry sector
p_sym = 1                       # parity symetry sector (even/odd)
x_sym = 1                       # spin-flip symmetry sector (only when h=0)
#---------------------------------------------------- PLOT SETTINGS DICTIONAIRY
"""
General settings for all plots
"""
plot_settings_dict = {
    'vs':             'L',          # set parameter on x-axis
    'scaling':        'w',          # set scaling parameter (changing in legend)

    'x_scale':      'log',       
    'y_scale':      'log',          
    
    'physical_units':   1,          # rescale by Heisenberg time?

#---- rescaling y-data
    'rescaleY':         0,      
    'func_y':       'exp',           # rescale function -> function(x, nu) (power-law = x^nu)    
    'nu_y':             -1,           # power of inversion
    
#---- rescaling x-axis
    'rescaleX':         0,          
    'func_x':       'power-law',     # rescale function -> function(x, nu) (power-law = x^nu)    
    'nu_x':             -1.0,           # power of inversion

#---- operator options
    'operator':         8,         # chosen operator according to order set in IsingModel.h
    'site':            L/2,          # chosen site for local operator
    'smoothed':        0,          # choose running-smoothed option

#---- instances set after
    'vs_idx':          -1,          # idx of vs option set after dict
    'scaling_idx':     -1,          # idx of scaling option set after dict
#- could be added later but better overwview if present here
    'func_x_name':     '',
    'func_y_name':     ''
}

parameter_critical = 0.37

#---- set array with parameters
params_arr = [L, J, g, h, w, k_sym, p_sym, x_sym, J0, g0]
#options = ['L', 'J', 'g', 'h', 'w', 'k']
names = ps.options
names.extend(['p','x','J0','x0'])

#---- DIR
model_dir = ""
if model == 0 or hamiltonian == 2:
    model_dir = f"disorder{kPSep}"
elif model == 1:
    model_dir = f"symmetries{kPSep}"
elif model == 2:
    model_dir = f"local_pert{kPSep}"
else:
    model_dir = ""

base_directory = f"..{kPSep}results{kPSep}"
if hamiltonian == 0:
    base_directory += f"ISING{kPSep}" + model_dir + (f"PBC{kPSep}" if BC else f"OBC{kPSep}") 
elif hamiltonian == 1:
    base_directory += f"HEISENBERG{kPSep}" + model_dir + (f"PBC{kPSep}" if BC else f"OBC{kPSep}") 
elif hamiltonian == 2:
    base_directory += f"QUANTUM_SUN{kPSep}" + model_dir
else:
    base_directory += ""

#---- INSTANCE OF PLOT SETTINGS CLASS --> plot_settings.py
plot_settings = ps.plot_settings_class(plot_settings_dict)


#---------------------------------------------------- OPERATOR OPTIONS
operator_names = [
    "SigmaZ_j=",
    "SigmaX_j=",
    "H_j=",
    "SigmaZ_q=",
    "SigmaX_q=",
    "H_q=",
    "TFIM_LIOM_plus_n=",
    "TFIM_LIOM_minus_n=",
    "SpinCurrent",
    "SigmaX",
    "SigmaZ",
    "SigmaX_near_neigh",
    "SigmaZ_near_neigh",
    "SigmaX_next_neigh",
    "SigmaZ_next_neigh",
    "SpinImbalance"			
]

def operator_name(operator, site):
    return operator_names[operator] + ("%d"%site if operator < 8 else "")
    
def subdir(operator, site):
    return (f"EXTENSIVE{kPSep}" if operator > 7 else ("j=%d%s" if operator < 3 else "q=%d%s")%(site, kPSep) )


operator_formuals = [
    r"$\sigma^z_j$",
    r"$\sigma^x_j$",
    (r"$H^{XXZ}_j=J_j\left(\sigma^x_j\sigma^x_{j+1} + \sigma^y_j\sigma^y_{j+1}\right) + \Delta_j\sigma^z_j\sigma^z_{j+1} + \frac{h_j}{2}\left(\sigma^z_j+\sigma^z_{j+1}\right)$" 
            if hamiltonian  else r"$H^{Ising}_j=J_j\sigma^z_j\sigma^z_{j+1} + \frac{g_j}{2}\left(\sigma^x_j+\sigma^x_{j+1}\right) + \frac{h_j}{2}\left(\sigma^z_j+\sigma^z_{j+1}\right)$"),
    r"$\sigma^z_q=\frac{1}{\sqrt{L}}\sum_\ell e^{i\frac{2\pi}{L}q\ell}\sigma^z_\ell$",
    r"$\sigma^x_q=\frac{1}{\sqrt{L}}\sum_\ell e^{i\frac{2\pi}{L}q\ell}\sigma^x_\ell$",
    r"$H_q=\frac{1}{\sqrt{L}}\sum_\ell \cos\left(\frac{2\pi}{L}q\ell\right)H_\ell$",
    "TFIM_LIOM_plus_n=...",
    "TFIM_LIOM_minus_n=...",
    r"$\frac{J}{\sqrt{L}}\sum_\ell\left(\sigma^x_j\sigma^y_{j+1}-\sigma^y_j\sigma^x_{j+1}\right)$",
    r"$\sum_\ell\sigma^x_\ell$",
    r"$\frac{1}{\sqrt{L}}\sum_\ell\sigma^z_\ell$",
    r"$\frac{1}{\sqrt{L}}\sum_\ell\sigma^x_\ell\sigma^z_{\ell+1}$",
    r"$\frac{1}{\sqrt{L}}\sum_\ell\sigma^z_\ell\sigma^z_{\ell+1}$",
    r"$\frac{1}{\sqrt{L}}\sum_\ell\sigma^x_\ell\sigma^z_{\ell+2}$",
    r"$\frac{1}{\sqrt{L}}\sum_\ell\sigma^z_\ell\sigma^z_{\ell+2}$",
    r"$\frac{1}{\sqrt{L}}\sum_\ell(-1)^\ell\sigma^z_\ell$"			
]

operator_names_latex = [
    "\sigma^z_j",
    "\sigma^x_j",
    "H_j",
    "\sigma^z_q",
    "\sigma^x_q",
    "H_q",
    "TFIM_LIOM_plus_n=...",
    "TFIM_LIOM_minus_n=...",
    "\mathcal{J}_{spin}",
    "\sigma^x_{tot}",
    "\sigma^z_{tot}",
    "U^x_n",
    "U^z_n",
    "U^x_{nn}",
    "U^z_{nn}",
    "\mathcal{I}_{spin}"
]


#---------------------------------------------------- SHORT USER FUNCTIONS
def set_params(_L = None, _J = None, _J0 = None, _g = None, _g0 = None, _h = None, _w = None, _k = None, _p = None, _x = None):
    global params_arr, L, J, J0, g, g0, h, w, k_sym, p_sym, x_sym
    if(_L != None):     L = _L 
    if(_J != None):     J = _J 
    if(_J0 != None):    J0 = _J0
    if(_g != None):     g = _g 
    if(_g0 != None):    g0 = _g0
    if(_h != None):     h = _h 
    if(_w != None):     w = _w 
    if(_k != None):     k_sym = _k 
    if(_p != None):     p_sym = _p 
    if(_x != None):     x_sym = _x
    params_arr = [L, J, g, h, w, k_sym, p_sym, x_sym, J0, g0]

def set_params2(**kwargs):
    print(kwargs)
    for key in kwargs:
        params_arr[ps.options.index(key)] = kwargs[key]

def set_plot_dict_value(option, value):
    plot_settings_dict[option] = value
