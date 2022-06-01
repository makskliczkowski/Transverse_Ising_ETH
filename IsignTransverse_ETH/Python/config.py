from os import sep as kPSep

#---------------------------------------------------- MODEL PARAMETERS
model = 0   # chooses model: 0-disorder / 1-symmetries
BC = 1     # boundaary condition: 0 - OBC / 1 - PBC

L = 12                          # system size
J = 1.00                        # spin exchange (Ising-like)
g = 0.90                        # trasnverse magnetic field (z-axis)
h = 0.80                        # longitudal magnetic field (x-axis)
#---- DISORDER PARAMETERS
w = 1.00                        # disorder on longitudonal field ( h_i \in [h-w, h+w] )
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
    'vs':             'g',          # set parameter on x-axis
    'scaling':        'L',          # set scaling parameter (changing in legend)

    'x_scale':      'linear',       
    'y_scale':      'log',          
    
    'physical_units':   1,          # rescale by Heisenberg time?

    'function':     'exp',          # rescale function
    'rescale':          0,          # inverse axis x -> 1/x^nu    
    'nu':               2,          # power of inversion

#---- instances set after
    'vs_idx':          -1,          # idx of vs option set after dict
    'scaling_idx':     -1           # idx of scaling option set after dict
}

#---- set array with parameters
params_arr = [L, J, g, h, w, k_sym, p_sym, x_sym, J0, g0]
options = ['L', 'J', 'g', 'h', 'w', 'k']
names = options
names.extend(['p','x','J0','x0'])

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
#---- DIR
base_directory = f"..{kPSep}results{kPSep}" + (f"symmetries{kPSep}" if model else f"disorder{kPSep}") + (f"PBC{kPSep}" if BC else f"OBC{kPSep}") 


#--- class for settings to control
class plot_settings_class:
    settings = {}
    def __init__(self, settings_dict_in):
        self.settings = settings_dict_in
        self.settings['vs_idx'] = options.index(self.settings['vs'])
        self.settings['scaling_idx'] = options.index(self.settings['scaling'])
        #if self.settings['rescale'] : self.settings['x_scale'] = 'log'

    def set_vs(self, vs_option):
        self.settings['vs'] = vs_option
        self.settings['vs_idx'] = options.index(self.settings['vs'])

    def set_scaling(self, scaling_option):
        self.settings['scaling'] = scaling_option
        self.settings['scaling_idx'] = options.index(self.settings['scaling'])

plot_settings = plot_settings_class(plot_settings_dict)