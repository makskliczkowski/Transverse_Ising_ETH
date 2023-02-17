import numpy as np
from scipy.optimize import differential_evolution
from . import critical_func as crit
from . import scaling_ansatz as rescale
import importlib
from random import seed as set_seed
from random import random
importlib.reload(crit)
importlib.reload(rescale)

import warnings

#suppress warnings
warnings.filterwarnings('ignore')

from inspect import signature, getmembers, isfunction
# --------------------------
# CRITICAL PARAMETER SCALING
# --------------------------

# which routines are available
_crit_functions_list = [fun for fun in getmembers(crit)
                        if isfunction(fun[1])]
# dictionary to be used
crit_functions_dict = {key.split('_crit_')[1]: value for
                        (key, value) in _crit_functions_list}
_crit_keys = crit_functions_dict.keys()

# -------------------
# RESCALING FUNCTIONS
# -------------------
_resc_functions_list = [fun for fun in getmembers(rescale)
                        if isfunction(fun[1])]

resc_functions_dict = {key.split('_rescale_')[1]: value for
                        (key, value) in _resc_functions_list}
_resc_keys = resc_functions_dict.keys()

# --------------------
# RESCALING FOR XLABEL
# --------------------

#------------------------------------------ XLABELS FOR SCALING ANSATZ
scale_ansatz_label = {
    'spacing':      lambda vs_str : "(" + vs_str + " - " + vs_str + "_c) \\cdot \\omega_H^{-1/\\nu}",
    'FGR':      lambda vs_str : "(" + vs_str + " - " + vs_str + "_c) \\cdot D^{1/\\nu}",
    'KT':       lambda vs_str : "L / \\xi_{KT}\ \ \ \ \\xi_{KT}^{-1}=exp(\\nu\\cdot|" + vs_str + " - " + vs_str + "_c|^{-1/2})",
    'RG':       lambda vs_str: "L / \\xi_0\ \ \ \ \\xi_0^{-1}=|" + vs_str + " - " + vs_str + "_c|^{\\nu}",
    'classic':  lambda vs_str: "(" + vs_str + " - " + vs_str + "_c) \\cdot L^{1/\\nu}",
    'exp':      lambda vs_str : "exp(\\nu\\cdot|" + vs_str + " - " + vs_str + "_c|) \\cdot e^{\\frac{ln2}{2}L}",
}


# CALCULATING THE COST FUNCTION FOR GIVEN SORTED DATASET
def calculate_cost_function(data):
    """
    Calculate cost function for given rescaled data
    """
    cost_func = 0
    for i in range(0, len(data) - 1):
        cost_func = cost_func + abs(data[i + 1] - data[i])
    return cost_func / ( max(data) - min(data) ) - 1


# FUNCTION TO BE MINIMIZED
def minimization_function(params, xvals, y, sizes, scaling_ansatz, crit_function, wH = None):
    """
    General function to used in minimization scheme
    Collects all data to a 1D array and sorts it according to non-decreasing values of the scaling ansatz
    """
    
    if crit_function not in _crit_keys:
        """
        Check if input choice is among possible critical scalings
        """
        err_message = ('Critical parameter scaling '
                       'function {} not allowed! '
                       'Allowed functions are: {}').format(crit_function,
                                                           _crit_keys)
        raise ValueError(err_message)

    if scaling_ansatz not in _resc_keys:
        """
        Check if input choice is among possible scaling solutions
        """
        err_message = ('Rescaling '
                       'function {} not allowed! '
                       'Allowed functions are: {}').format(scaling_ansatz,
                                                           _resc_keys)
        raise ValueError(err_message)
    
    if wH is None and scaling_ansatz == 'spacing':
        """
        Check if level spacing is input if spacing rescale is chosen
        """
        err_message = ( 'Mean level spacing is empty sequence'
                        'Need data of shape: {}').format(y.shape())
        raise ValueError(err_message)

    rescale_fun = resc_functions_dict[scaling_ansatz]
    crit_fun = crit_functions_dict[crit_function]

    crit_pars=np.array(params[1:]) # critical parameters to find given in critical point scaling

    xdata = []
    fulldata = []
    final_func = None
    if scaling_ansatz == 'spacing':   
        final_func = lambda params, i, j : rescale_fun(xvals[i][j], sizes[i], crit_fun, *crit_pars) * wH[i][j]**(-1. / params[0])
    else:                               
        final_func = lambda params, i, j : rescale_fun(xvals[i][j], sizes[i], crit_fun, params[0], *crit_pars)

    center_of_range = []
    for i in range(0, len(sizes)):
        center_of_range.append( final_func(params, i, int(len(xvals[i]) / 2)) )
        for j in range(0, len(xvals[i])):
            xvalue = final_func(params, i, j)
            xdata.append( xvalue )
            fulldata.append(y[i][j])
    
    fulldata = np.array(fulldata)
    permut = np.argsort(xdata)
    fulldata = fulldata[permut]

    # constraint
    """
    If rescaled datapoints for any system size exceed the rest of
    the data the constrint is not satisfied. In other words, the center of each range has to be inside the other ranges.
    """
    constraint_satisfied = True
    for center in center_of_range:
        for i in range(len(sizes)):
            if center < final_func(params, i, 0) or center > final_func(params, i, len(xvals[i])-1):
                constraint_satisfied = False
                break


    cost_fun = calculate_cost_function(fulldata)
    #if cost_fun < 1.: print(cost_fun)
    return cost_fun if constraint_satisfied else 1000 + (random() * 1e10)



def cost_func_minization(x, y, sizes, bnds, 
                            scale_func, crit_func, 
                            population_size=1e2, maxiterarions=1e3, 
                            workers=1, seed = None, wH = None):
    """
    Main function returning optimal parameters for given scaling ansatz

    Parameters:
    ----------------------
    x: 2D array 
        values of parameters (e.g. disorder) 
        at different system sizes
    
    y: 2D array 
        qunatity to find scaling solution for
        at different parameter values and sizes
    
    sizes: 1D array of system sizes present in scaling

    bnds: tuple of bounds for each fit parameter

    crit_func: string
        String to specify which critical function should be chosen
        (gives ValueError when is not amongst existing ones)

    scale_func: string
        String to specify which scaling solution should be chosen
        (gives ValueError when is not amongst existing ones)

    """

    if seed is None: seed = np.random.default_rng();
    set_seed(seed)

    result = differential_evolution(
                    minimization_function,
                    bounds=bnds,
                    args=(x, y, sizes, scale_func, crit_func, wH),
                    popsize=int(population_size), 
                    maxiter=int(maxiterarions), 
                    workers=workers, atol=1e-3,
                    seed=seed
            )
    optimal_res = np.array(result.x)
    cost_fun = result.fun
    if result.success ==  False: print('Failed convergence')
      
    return optimal_res, cost_fun, result.success


def prepare_bounds(x, crit_fun, scaling_ansatz, vals):
    """
    Generate bounds for optimization procedure
    
    For details about signature of critical functions see critical_fun.py to see
    what the boudns correspond to

    Parameters:
    ------------
        x: 2D array 
            values of parameters (e.g. disorder) 
            at different system sizes

        crit_fun:   'string'
            function of critical value, one among critical_func.py

        vals:   1D array
            scaling parameters: system sizes
    """
    x_max = -1e6
    x_min = 1e6
    for a in x: 
        for _x_ in a: 
            if x_max is None or _x_ > x_max:   x_max = _x_
            if x_min is None or _x_ < x_min:    x_min = _x_
            
    bounds = [(0.2, 5.)] if scaling_ansatz == 'FGR' or scaling_ansatz == 'spacing' else [(0.1, 10.)]
    #-- number of bounds is number of different scaling parameters
    if crit_fun == 'free':  
        for i in range(len(vals)):  
            bounds.append((x_min, x_max))
    elif crit_fun == 'free_inv':  
        for i in range(len(vals)):  
            bounds.append((1. / x_max, 1. / x_min))
    #-- constant value, only one new bounds
    elif crit_fun == 'const':   
        bounds.append((0, x_max))
    #-- power law and inverse function have 3 paramters
    elif crit_fun == 'power_law' or crit_fun == 'inv':   
        for i in range(3):
            bounds.append((-10., 10.))
    #-- logarithmic and its inverse function have 2 paramters
    elif crit_fun == 'log' or crit_fun == 'inv_log' or crit_fun == 'lin':   
        for i in range(2):
            bounds.append((-10., 10.))
    return bounds


def get_crit_points(x, y, vals, crit_fun='free', scaling_ansatz = 'classic', seed = None, wH = None):
    """
    Calculating critical points for each system size

    Parameters:
    -----------
        x, y:   2D array
            x and y data for collapse with different vals (sizes) values
        
        vals:   1D array
            scaling parameters: system sizes
        
        crit_fun:   'string'
            function of critical value, one among critical_func.py

        scaling_ansatz:   'string'
            scaling ansatz for collapse, classic, KT, FGR,...
    """
    
    params = []
    bounds = prepare_bounds(x, crit_fun, scaling_ansatz, vals)

    params, cost_fun, status = cost_func_minization(x=x, y=y, sizes=vals,
                                    scale_func=scaling_ansatz, 
                                    crit_func=crit_fun,
                                    bnds=bounds,
                                    population_size=1e2,
                                    maxiterarions=1e3, workers=10,
                                    seed = seed,
                                    wH = wH
                                )

    par = params[0]
    crit_pars = np.array(params[1:])
    return par, crit_pars, cost_fun, status