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
p_sym = 1                       # 
x_sym = 1                       # 

#---- set array with parameters
params_arr = [L, J, g, h, w, k_sym, p_sym, x_sym, J0, g0]

#---- DIR
base_directory = f"..{kPSep}results{kPSep}" + (f"symmetries{kPSep}" if model else f"disorder{kPSep}") + (f"PBC{kPSep}" if BC else f"OBC{kPSep}") 

#---------------------------------------------------- PLOT SETTINGS DICTIONAIRY
"""
General settings for all plots
"""
plot_settings = {
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
if plot_settings['rescale'] : plot_settings['x_scale'] = 'log'

_options = ['L', 'J', 'g', 'h', 'w', 'k']

plot_settings['vs_idx'] = _options.index(plot_settings['vs'])
plot_settings['scaling_idx'] = _options.index(plot_settings['scaling'])

