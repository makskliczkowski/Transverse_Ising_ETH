import importlib
from os import sep as kPSep
from os.path import exists
import numpy as np

import costfun.costfun as cost
import utils.helper_functions as hfun
import pandas as pd
import config as cf
import copy
importlib.reload(cf)
importlib.reload(cost)
importlib.reload(hfun)
#--- Global
user_settings = getattr(cf.plot_settings, 'settings')

# ------------------------------------------------------------------------------------------------------------------------


def GOE(x : np.array):
    """
    GOE shape of sff in thermodynamic limit
    
    Parameters:
    -----------------
        x : np.array
            numpy array with datapoints (times defined for unfolded data)
    """
    return np.array([2 * a - a * np.log(1 + 2 * a) if a < 1 else 2 - a * np.log( (2 * a + 1) / (2 * a - 1) ) for a in x])


# ------------------------------------------------------------------------------------------------------------------------

def load(settings = None, parameter = None):
    """
    Load spectral data along with statistical measures

    Parameters:
    -----------------
        settings : dict
            settings dictionairy as in config.py with details

        parameter : double/int
            current parameter value of 'scaling' parameter
    """
    if settings == None:
        settings = user_settings

    param_copy = copy.deepcopy(cf.params_arr)
    if parameter == None:
        raise ValueError("Input value 'parameter' unasigned. No default value.")

    cf.params_arr[settings['scaling_idx']] = parameter
    if settings['scaling_idx'] == 3 and cf.J0 == 0 and cf.g0 == 0:
        cf.params_arr[4] = int(100 * parameter / 2.) / 100.

    dir = cf.base_directory + "SpectralFormFactor" + kPSep + ("smoothed" + kPSep if settings['smoothed'] else "")
    filename = dir + hfun.info_param(cf.params_arr)
    info = hfun.info_param(cf.params_arr)
    print(info, hfun.get_params_from_info(info))
    #--- reset defaults
    cf.params_arr = param_copy

    if exists(filename):
        data = pd.read_table(filename, sep="\t", header=None)
        xdata = np.array(data[0])
        ydata = np.array(data[1])
        tH = data[2][0]
        tau = data[3][0]
        gap_ratio = data[4][0]
        dim = data[6][0]
        return True, xdata, ydata, tH, tau, gap_ratio, dim
    else:
        return False, np.array([]), np.array([]), None, None, None, None

# ------------------------------------------------------------------------------------------------------------------------

def plot(axis, settings = None, 
                    xlab = "\\tau", ylab = "K(\\tau)", xscale='log', yscale="log",
                    func_x = lambda x, a: x, func_y = lambda y, a: y, 
                    font = 12, vals = None, folded = False, axis_inset = None) :
    """
    Plot spectral form factor according to input range

    Parameters:
    -----------------
        settings : dictionairy
            general settings for plotting and rescaling given in config.py or user input

        xlab, ylab : r-string (LaTeX)
            label for x, y axis respectively

        xscale, yscale : string
            axis scale for x, y respectively

        func_x, func_y : lambda
            rescaling function for x, y axis values (input function is x, y data and scaling parameter)

        font : int
            base font size (legend and axis labels are set accordingly)
        
        folded : boolean
            choose whether to use folded spectrum not rescaled by mean level spacing

        vals : np.array
            Numpy Array with scaling parameter values to sweep through

    """

    #-- main settings
    if settings == None:
        settings = user_settings
    #-- axis scales
    if xscale != None: settings['x_scale'] = xscale
    if yscale != None: settings['y_scale'] = yscale

    #--- prepare scaling - axis
    vals = np.array(vals)
    if vals.any() == None:
        vals = hfun.get_scaling_array(settings=settings)
    taus = [];  val_at_taus = []
    gap_ratio = []; dim = 1;    times = []; y_min = 1e10;

    for x in vals:

        status, xdata, ydata, tH, tau, r, dim = load(settings=settings, parameter=x)
        
        if status:
            dim = dim
            times = xdata
            gap_ratio.append(r)
            xdata = func_x(xdata, x)
            ydata = func_y(ydata, x)
            tau = func_x(tau, x);   taus.append(tau)
            idx = min(range(len(xdata)), key=lambda i: abs(xdata[i] - tau));  val_at_taus.append(ydata[idx])
            if min(ydata) < y_min:  y_min = min(ydata)

            axis.plot(xdata, ydata, label=hfun.key_title(x, settings), linewidth=int(font / 6), markersize=font-6, zorder=int(10*x))
        else:
            taus.append(np.nan)
            val_at_taus.append(np.nan)
            gap_ratio.append(np.nan)
    
    hfun.set_plot_elements(axis = axis, xlim = (func_x(1. / (2 * np.pi * dim), min(vals)), func_x(10, max(vals))), 
                                    ylim = (None, None), ylabel = ylab, xlabel = xlab, settings=settings, font_size=font, set_legend=True)
    axis.legend(loc='lower right')
    axis.set_ylim(0.75 * y_min, None)
    title = ""
    if (settings['vs_idx'] == 3 or settings['scaling_idx'] == 3) and cf.J0 == 0 and cf.g0 == 0 and cf.h != 0:
        title = hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs'], settings['scaling'], 'w') + ',w=0.5h'
    else :
        title = hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs'], settings['scaling'])
    if settings['vs_idx'] != 2 :
        try : 
            title = list(title);    title[title.index('g')] = hfun.var_name;   title = "".join(title) # g
            #title = list(title);    title[title.index('g')] = hfun.var_name;   title = "".join(title) # g0
        except ValueError:
                print("not found")
    axis.title.set_text(r"$%s$"%title[1:])
    axis.title.set_fontsize(10)
    
    axis.scatter(taus, val_at_taus, marker='o', edgecolor='black', s=100, facecolors='none', zorder=1000)
    axis.plot(times, GOE(times), linestyle='--', color='black')
    if axis_inset is not None:
        axis_inset.plot(vals, gap_ratio, marker = 'o', linestyle='--', color='black', linewidth=int(font / 6), markersize=font-4)
        axis_inset.axhline(y=0.5307, ls='--', color='black')
        axis_inset.axhline(y=0.3867, ls='--', color='black')

# ------------------------------------------------------------------------------------------------------------------------

def plot_deviation(axis, settings = None, 
                    xlab = "\\tau", ylab = "K(\\tau)", xscale='log', yscale="log",
                    font = 12, vals = None) :
    """
    Plot spectral form factor according to input range

    Parameters:
    -----------------
        settings : dictionairy
            general settings for plotting and rescaling given in config.py or user input

        xlab, ylab : r-string (LaTeX)
            label for x, y axis respectively

        xscale, yscale : string
            axis scale for x, y respectively

        font : int
            base font size (legend and axis labels are set accordingly)

        vals : np.array
            Numpy Array with scaling parameter values to sweep through

    """

    #-- main settings
    if settings == None:
        settings = user_settings
    #-- axis scales
    if xscale != None: settings['x_scale'] = xscale
    if yscale != None: settings['y_scale'] = yscale

    #--- prepare scaling - axis
    vals = np.array(vals)
    if vals.any() == None:
        vals = hfun.get_scaling_array(settings=settings)
    dim = 1;    times = []; y_min = 1e10;

    for x in vals:

        status, xdata, ydata, tH, tau, r, dim = load(settings=settings, parameter=x)
        
        if status:
            dim = dim
            ydata = np.log10(ydata / GOE(xdata))
            if min(ydata) < y_min:  y_min = min(ydata)

            axis.plot(xdata, ydata, label=hfun.key_title(x, settings), linewidth=int(font / 6), markersize=font-6, zorder=int(10*x))
    if dim is None: dim = 1e4
    
    axis.axhline(y=1e-1, linestyle='--', color='black')
    hfun.set_plot_elements(axis = axis, xlim = (1. / (2 * np.pi * dim), 10.0), 
                                    ylim = (0.75 * y_min, 10), ylabel = ylab, xlabel = xlab, settings=settings, font_size=font, set_legend=True)
    axis.legend(loc='upper right')
    axis.set_ylim(0.75 * y_min, 10)
    title = ""
    if (settings['vs_idx'] == 3 or settings['scaling_idx'] == 3) and cf.J0 == 0 and cf.g0 == 0 and cf.h != 0:
        title = hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs'], settings['scaling'], 'w') + ',w=0.5h'
    else :
        title = hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs'], settings['scaling'])
    if settings['vs_idx'] != 2 :
        try : 
            title = list(title);    title[title.index('g')] = hfun.var_name;   title = "".join(title) # g
            #title = list(title);    title[title.index('g')] = hfun.var_name;   title = "".join(title) # g0
        except ValueError:
                print("not found")
    axis.title.set_text(r"$%s$"%title[1:])
    axis.title.set_fontsize(10)
    