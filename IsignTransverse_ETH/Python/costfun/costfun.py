import numpy as np
from scipy.optimize import differential_evolution
from . import critical_func as crit
from . import scaling_ansatz as rescale
import importlib
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
    'FGR':      lambda vs_str : "|" + vs_str + " - " + vs_str + "_c|^{\\nu} \\cdot e^{\\frac{ln2}{2}L}",
    'KT':       lambda vs_str : "L / \\xi_{KT}\ \ \ \ \\xi_{KT}^{-1}=exp(\\nu\\cdot|" + vs_str + " - " + vs_str + "_c|^{-1/2})",
    'RG':       lambda vs_str: "L / \\xi_0\ \ \ \ \\xi_0^{-1}=|" + vs_str + " - " + vs_str + "_c|^{\\nu}",
    'classic':  lambda vs_str: "(" + vs_str + " - " + vs_str + "_c) \\cdot L^{\\nu}",
    'exp':      lambda vs_str : "exp(\\nu\\cdot|" + vs_str + " - " + vs_str + "_c|) \\cdot e^{\\frac{ln2}{2}L}",
}


# CALCULATING THE COST FUNCTION FOR GIVEN SORTED DATASET
def calculate_cost_function(data):
    """
    Calculate cost function for given rescaled data
    """
    cost_func = 0
    for i in range(0, len(data) - 1):
        cost_func = cost_func + abs(data[i+1] - data[i])
    return cost_func / ( max(data) - min(data)) - 1


# FUNCTION TO BE MINIMIZED
def minimization_function(params, xvals, y, sizes, scaling_ansatz, crit_function):
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
    
    rescale_fun = resc_functions_dict[scaling_ansatz]
    crit_fun = crit_functions_dict[crit_function]
    
    crit_pars=np.array(params[1:]) # critical parameters to find given in critical point scaling

    xdata = []
    fulldata = []
    
    for i in range(0, len(sizes)):
        for j in range(0, len(xvals[i])):
            xvalue = rescale_fun(
                    xvals[i][j], 
                    sizes[i], 
                    crit_fun, 
                    params[0], *crit_pars)
            xdata.append( xvalue )
            fulldata.append(y[i][j])
    
    fulldata = np.array(fulldata)
    permut = np.argsort(xdata)
    fulldata = fulldata[permut]

    cost_fun = calculate_cost_function(fulldata)
    #if cost_fun < 1.: print(cost_fun)
    return cost_fun



def cost_func_minization(x, y, sizes, bnds, 
                            scale_func, crit_func, 
                            population_size=1e2, maxiterarions=1e3, 
                            workers=1, realisations=1, seed = None):
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
    result = differential_evolution(
                    minimization_function,
                    bounds=bnds,
                    args=(x, y, sizes, scale_func, crit_func),
                    popsize=int(population_size), 
                    maxiter=int(maxiterarions), 
                    workers=workers, atol=1e-2,
                    seed=seed
            )
    optimal_res = np.array(result.x)
    cost_fun = result.fun
    if realisations > 1:
        for r in range(realisations - 1):
            seed_r = np.random.default_rng();
            result = differential_evolution(
                        minimization_function,
                        bounds=bnds,
                        args=(x, y, sizes, scale_func, crit_func),
                        popsize=int(population_size), 
                        maxiter=int(maxiterarions), 
                        workers=workers, atol=1e-2,
                        seed=seed_r
                )
            optimal_res = optimal_res + np.array(result.x)
            cost_fun = cost_fun + result.fun
            if result.success ==  False: print('Failed convergence')    
    return optimal_res / float(realisations), cost_fun / float(realisations)


def get_crit_points(x, y, vals, scaling_ansatz = None, seed = None):
    """
    Calculating critical points for each system size

    Parameters:
    -----------
        x, y:   1D array
            x and y data for collapse
        
        vals:   1D array
            scaling parameters: system sizes
        
        scaling_ansatz:   'string'
            scaling ansatz for collapse, classic, KT, FGR,...
    """
    if scaling_ansatz is None: scaling_ansatz='classic'

    crit_fun='free'
    params = []
    x_max=None
    for a in x: 
        for _x_ in a: 
            if x_max is None or _x_ > x_max: x_max = _x_

    bounds = [(0, x_max) for i in range(len(vals) + 1)]
    params, cost_fun = cost_func_minization(x=x, y=y, sizes=vals,
                                    scale_func=scaling_ansatz, 
                                    crit_func=crit_fun,
                                    bnds=bounds,
                                    population_size=1e2,
                                    maxiterarions=1e2, workers=10, realisations=1,
                                    seed = seed
                                )

    par = params[0]
    crit_pars = np.array(params[1:])
    return par, crit_pars, cost_fun