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
from scipy.signal import savgol_filter

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

def load(settings = None, parameter = None, folded = False, beta = 0.0):
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
    Lx = cf.params_arr[0]


    dir = cf.base_directory + "SpectralFormFactor" + kPSep
    if beta > 0:
        n = hfun.order_of_magnitude(beta) 
        dir += str("beta={:.%df}"%n).format(round(beta, n)) + kPSep
    #dir += "smoothed" + kPSep

    filename = dir + ("folded" if folded else "") + hfun.info_param(cf.params_arr)
    info = hfun.info_param(cf.params_arr)
    
    # ------ GET STATISTICAL DATA
    filename_stats = cf.base_directory + "STATISTICS" + kPSep + "raw_data" + kPSep + hfun.info_param(cf.params_arr)
    wH = 0; gap_ratio = 0
    stats = {}
    try:
        stats = hfun.load_stats(filename_stats)
        gap_ratio = stats['gap ratio']
        wH = stats['mean level spacing']
    except FileNotFoundError:   wH = 0.0
    if np.isnan(wH):            wH = 0.0
    if wH > 0: tH = 1. / wH

    #--- reset defaults
    cf.params_arr = param_copy

    epsilon = 1e-1
    if exists(filename):
        data = pd.read_table(filename, sep="\t", header=None)
        times = np.array(data[0])
        sff = np.array(data[1])
        dim = data[6][0]

        if folded:
            wH_eh = np.sqrt(cf.params_arr[0]) / (0.341345 * dim) * np.sqrt(cf.params_arr[1]**2 + cf.params_arr[3]**2 + cf.params_arr[2]**2
												 + ( 0.0 if cf.model == 0 else (cf.params_arr[4]**2 + cf.params_arr[8]**2 + cf.params_arr[9]**2) / 3. ))
            tH_eh = 1./ wH_eh 
            time_end = int(np.ceil(np.log10(5 * tH_eh)))
            time_end = time_end + 1 if time_end / np.log10(tH_eh) < 1.5 else time_end
            times = np.logspace(-2, time_end, sff.size)

        if settings['smoothed'] == 1:    
            sff = savgol_filter(sff, window_length=int(0.1 * sff.size) + int(0.1 * sff.size) % 2 - 1, polyorder=5, mode="mirror")
            sff = hfun.remove_fluctuations(sff, 200)
            #idx = min(range(len(times)), key=lambda i: abs(times[i] - 3.0)); 
            #sff /= sff[idx]
        if wH == 0.0: 
            tH = data[2][0]
            print("Not found stats for wH")
        #tau = data[3][0]
        times_for_algorithm = times / tH if folded else times
        sff_dev = np.abs(np.log10(sff / GOE(times_for_algorithm)))
        for i, K in reversed(list(enumerate(sff_dev))):
            if K > epsilon and times[i] < (3 * tH if folded else 3):
                tau = times[i-1]
                break

        if gap_ratio == 0:
            gap_ratio = data[4][0]
            print("Not found stats for gap ratio")
        return True, times, sff, tH, tau, gap_ratio, dim
    else:
        print(filename, hfun.get_params_from_info(info))
        return False, np.array([]), np.array([]), None, None, None, None

# ------------------------------------------------------------------------------------------------------------------------

def plot(axis, settings = None, 
                    xlab = "\\tau", ylab = "K(\\tau)", xscale='log', yscale="log",
                    func_x = lambda x, a: x, func_y = lambda y, a: y, 
                    font = 20, vals = None, folded = False, axis_inset = None, zoom = False, beta = 0.0) :
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
        Lx = x if settings['scaling_idx'] == 0 else cf.L
        status, times, sff, tH, tau, r, dimensions = load(settings=settings, parameter=x, folded = folded, beta=beta)

        if status:

            dim = dimensions
            gap_ratio.append(r)
            sff = func_y(sff, x)
            tau = func_x(tau, x)
            taus.append(tau)
            idx = min(range(len(times)), key=lambda i: abs(func_x(times, x)[i] - tau));  val_at_taus.append(sff[idx])
            if min(sff) < y_min:  y_min = min(sff)

            key_tit = r"$\beta=%.2f$"%beta if len(vals) == 1 else hfun.key_title(x, settings)
            p = axis.plot(func_x(times, x), sff, label=key_tit, linewidth=int(font / 6), markersize=font-6, zorder=int(10*x))
            axis.scatter(tau, sff[idx], marker='o', edgecolor=p[0].get_color(), s=100, facecolors='none', zorder=1000)
            if folded:
                axis.plot(func_x(times, x), GOE(times / tH), linestyle='--', color='black')
        else:
            taus.append(np.nan)
            val_at_taus.append(np.nan)
            gap_ratio.append(np.nan)
    print(dim, vals)
    hfun.set_plot_elements(axis = axis, xlim = (func_x(min(times), min(vals)), func_x(0.6 * max(times), max(vals))),
                                    ylim = (None, None), ylabel = ylab, xlabel = xlab, settings=settings, font_size=font, set_legend=True)
    #axis.legend(loc='lower right')
    axis.set_ylim(0.75 * y_min, None)
    if zoom:
        axis.set_ylim(0.75 * y_min, 1.2)
        axis.set_yscale('linear')
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
    #axis.title.set_text(r"$%s$"%title[1:])
    #axis.title.set_fontsize(10)
    
    #axis.scatter(taus, val_at_taus, marker='o', edgecolor='black', s=100, facecolors='none', zorder=1000)
    if folded == False:
        axis.plot(func_x(times, x), GOE(times), linestyle='--', color='black')
        axis.axvline(x=1.0, linestyle='--', color='black', ymin=0, ymax=0.4)
    axis.tick_params(axis="both",which='major',direction="in",length=6)
    axis.tick_params(axis="both",which='minor',direction="in",length=3)
    if axis_inset is not None:
        axis_inset.plot(vals, gap_ratio, marker = 'o', linestyle='--', color='black', linewidth=int(font / 6), markersize=font-10)
        axis_inset.axhline(y=0.5307, ls='--', color='black')
        axis_inset.axhline(y=0.3867, ls='--', color='black')

# ------------------------------------------------------------------------------------------------------------------------

def plot_deviation(axis, settings = None, 
                    xlab = "\\tau", ylab = "K(\\tau)", xscale='log', yscale="log",
                    font = 12, vals = None, func_x = lambda x, a: x, folded = False) :
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

        func_x : lambda
            rescaling function for x, y axis values (input function is x, y data and scaling parameter)

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

        status, times, sff, tH, tau, r, dim = load(settings=settings, parameter=x)
        
        if status:
            dim = dim
            sff = np.abs(np.log10(sff / GOE(times)))
            if min(sff) < y_min:  y_min = min(sff)
            if folded: times *= tH
            axis.plot(func_x(times, x), sff, label=hfun.key_title(x, settings), linewidth=int(font / 6), markersize=font-6, zorder=int(10*x))
    if dim is None: dim = 1e4
    
    axis.axhline(y=3e-1, linestyle='--', color='black')
    hfun.set_plot_elements(axis = axis, xlim = (func_x(min(times), min(vals)), func_x(0.6 * max(times), max(vals))),
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
    