import importlib
from os import sep as kPSep
import plot_settings as ps
importlib.reload(ps)
#---------------------------------------------------- MODEL PARAMETERS
model = 0           # chooses model: 0-disorder / 1-symmetries
hamiltonian = 1     # which hamiltonian?: 0-Ising / 1-Heisenberg
BC = 1              # boundaary condition: 0 - OBC / 1 - PBC

L = 16                          # system size
J = 1.00                        # spin exchange (Ising-like)
g = 1.00                        # trasnverse magnetic field (z-axis)
h = 0.00                        # longitudal magnetic field (x-axis)
#---- DISORDER PARAMETERS
w = 1.5                        # disorder on longitudonal field ( h_i \in [h-w, h+w] )
J0 = 0.2                        # disorder on spin exchange ( J_i \in [J-J0, J+J0] )
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
    'vs':             'g',          # set parameter on x-axis
    'scaling':        'L',          # set scaling parameter (changing in legend)

    'x_scale':      'log',       
    'y_scale':      'log',          
    
    'physical_units':   1,          # rescale by Heisenberg time?

#---- rescaling y-data
    'rescaleY':         0,      
    'func_y':       'log',           # rescale function -> function(x, nu) (power-law = 1 / x^nu)    
    'nu_y':             3,           # power of inversion
    
#---- rescaling x-axis
    'rescaleX':         0,          
    'func_x':       'power-law',     # rescale function -> function(x, nu) (power-law = 1 / x^nu)    
    'nu_x':             -0.5,           # power of inversion
    
#---- instances set after
    'vs_idx':          -1,          # idx of vs option set after dict
    'scaling_idx':     -1,          # idx of scaling option set after dict
#- could be added later but better overwview if present here
    'func_x_name':     '',
    'func_y_name':     ''
}

#---- set array with parameters
params_arr = [L, J, g, h, w, k_sym, p_sym, x_sym, J0, g0]
#options = ['L', 'J', 'g', 'h', 'w', 'k']
names = ps.options
names.extend(['p','x','J0','x0'])

#---- DIR
base_directory = f"..{kPSep}results{kPSep}" + (f"Heisenberg{kPSep}" if hamiltonian else f"Ising{kPSep}")\
                                             + (f"symmetries{kPSep}" if model else f"disorder{kPSep}") \
                                              + (f"PBC{kPSep}" if BC else f"OBC{kPSep}") 

#---- INSTANCE OF PLOT SETTINGS CLASS --> plot_settings.py
plot_settings = ps.plot_settings_class(plot_settings_dict)


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