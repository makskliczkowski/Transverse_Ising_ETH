from turtle import width
import numpy as np
import importlib
import helper_functions as hfun
import config as cf
import copy
importlib.reload(cf)
importlib.reload(hfun)
import pandas as pd
from os import sep as kPSep
from os.path import exists

#--- Global
user_settings = getattr(cf.plot_settings, 'settings')
#--- SET SCALING RANGES AND DATA

def get_scaling_array(settings = None, x0 = 0.1, xend = 1.0, dx = 0.1):
    if settings == None:
        settings = user_settings
    vals = []
    length = int((xend-x0) / dx) + 1
    if settings['scaling_idx'] == 0:
        if cf.hamiltonian: vals = range(12, 19, 2)
        else: vals = range(11, 17, 1)
    elif settings['scaling_idx'] == 5:
        vals = range(1, int(cf.params_arr[0] / 2) + 1)
    else :
        for x in range(0, length) :
            vals.append(x0 + x * dx)
    return np.array(vals)


def load_spectral(dir = "", settings = None, parameter = None, 
                            spec = None, normalise = False, 
                            func_x = lambda x, a: x,
                            operator = -1, site = -3):
    """
    Load spectral data along with statistical measures

    Parameters:
    -----------------
        dir : f-string
            directory of file to plot

        settings : dictionairy
            general settings for plotting and rescaling given in config.py or user input

        parameter : double/int
            current parameter value of 'scaling' parameter

        normalise : boolean
            choose if function to be rescaled or not
        
        func_x : lambda
            rescaling function for x axis values (input function is x data and scaling parameter)

        operator, site : int
            Values definining which operator to choose, if None the default from config.py is chosen
    """
    importlib.reload(cf)

    if operator < 0: operator = settings['operator']
    if site < -1: site = settings['site']

    param_copy = copy.deepcopy(cf.params_arr)
    if parameter == None:
        raise ValueError("Input value 'parameter' unasigned. No default value.")
    
    cf.params_arr[settings['scaling_idx']] = parameter
    if settings['scaling_idx'] == 3 and cf.J0 == 0 and cf.g0 == 0:
        cf.params_arr[4] = int(100 * parameter / 2.) / 100.
    filename = (hfun.info_param(cf.params_arr) if cf.hamiltonian else hfun.remove_info(hfun.info_param(cf.params_arr), 'J') + ".dat")
    
    if settings['scaling_idx'] == 5 and operator < 8:
        filename = dir + ("j=%d%s" if operator < 3 else "q=%d%s")%(parameter, kPSep) + cf.smo_dir + cf.operator_names[operator] + "%d"%parameter + filename
    elif settings['scaling_idx'] == 0 and site < 0 and operator < 8:
        filename = dir + ("j=%d%s" if operator < 3 else "q=%d%s")%(parameter / 2, kPSep) + cf.smo_dir + cf.operator_names[operator] + "%d"%(parameter/2) + filename
    else :
        filename = dir + cf.subdir(operator, site, settings['smoothed']) + cf.operator_name(operator, site) + filename
    
    filename2 = cf.base_directory + "STATISTICS" + kPSep + "raw_data" + kPSep + hfun.info_param(cf.params_arr)
    #--- reset defaults
    cf.params_arr = param_copy
    #print(filename)

    if exists(filename):
        seper = "\t\t" if spec == "spec" and cf.hamiltonian == 0 else "\t";
        data = pd.read_table(filename, sep=seper, header=None)
        stats = []
        if exists(filename2): stats = pd.read_table(filename2, sep="\t", header=None)
        xdata = func_x(np.array(data[0]), parameter)
        ydata = np.array(data[1])
        
        if normalise and spec != "spec":
            norm_idx = min(range(len(xdata)), key=lambda i: abs(xdata[i] - 0.1))
            #if x > 0.4 or spec == "time": 
            norm_idx = len(ydata)-1
            if spec == "time": ydata = (ydata - data[3][0]) / np.abs(ydata[norm_idx] - ydata[0])
            else: ydata = (ydata - data[3][0]) / np.abs(ydata[norm_idx] - ydata[0])
        if exists(filename2):
            "mean"
            wH = (np.array(stats[1][4])).astype(np.float);    
            "typical"
            wHtyp = (np.array(stats[1][5])).astype(np.float); 
        else:
            wH = 1e-8; wHtyp = 1e-8
        if spec == "time": wH = 1. / wH
        if spec == "time": wHtyp = 1. / wHtyp  
        return True, xdata, ydata, wH, wHtyp
    else:
        return False, np.array([]), np.array([]), None, None



def plot_spectral(axis, settings = None, 
                    xlab = "x", ylab = "y", xscale='log', yscale=None,
                    func_x = lambda x, a: x, func_y = lambda y, a: y,
                    normalise=False, spec="time", 
                    font = 12, use_derivative = 0, 
                    vals = None,
                    operator = -1, site = -3):
    """
    Plot spectral function according to input range

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

        normalise : boolean
            choose if function to be rescaled or not

        spec : string
            choose spectral to plot: "time", "int" or "spec"

        font : int
            base font size (legend and axis labels are set accordingly)
        
        use_derivative : boolean
            choose whether to use derivative of other function or raw_data (only valid for "spec" option)

        vals : np.array
            Numpy Array with scaling parameter values to sweep through

        operator, site : int
            Values definining which operator to choose, if None the default from config.py is chosen

    """
    if spec == "time":
                            dir = cf.base_directory + "timeEvolution%s"%kPSep
    elif spec == "int":     dir = cf.base_directory + "IntegratedResponseFunction%s"%kPSep
    elif spec == "spec":    dir = cf.base_directory + ("IntegratedResponseFunction%sDERIVATIVE%s"%(kPSep,kPSep) if use_derivative else "ResponseFunction%s"%kPSep)
    else:
        raise ValueError("No spectral data possible for this option, choose among: 'time', 'int' or 'spec'")

    if operator < 0: operator = settings['operator']
    if site < -1: site = settings['site']

    #-- main settings
    if settings == None:
        settings = user_settings
    #-- axis scales
    if xscale != None: settings['x_scale'] = xscale
    if yscale != None: settings['y_scale'] = yscale

    #--- prepare scaling - axis
    vals = np.array(vals)
    if vals.any() == None:
        vals = get_scaling_array(settings=settings)

    y_min = 1.0e10;     y_max = -1.0e10;
    x_min = 1.0e10;     x_max = -1.0e10;
    #--- load data and plot one-by-one
    wH = [];    LTA = []
    wH_typ = [];    val_at_typ = [];
    for x in vals:

        status, xdata, ydata, wHnow, wHtypnow = load_spectral(dir=dir, 
                                                    settings=settings, 
                                                    parameter=x,
                                                    spec=spec,
                                                    func_x=func_x,
                                                    normalise=normalise,
                                                    operator = operator,
                                                    site = site
                                                    )

        if status:
            
            ydata = func_y(ydata, x)

            if use_derivative == 0 and spec == "spec": 
                ydata = ydata * (2**x / x if settings['scaling_idx'] == 0 else 2**cf.L / cf.L) # rescale by D

            #idx_cut = 0
            #if use_derivative == 1 and spec == "spec": idx_cut = 200
            #xdata = np.array([xdata[i] for i in range(len(xdata)) if i > idx_cut])
            #ydata = np.array([ydata[i] for i in range(len(ydata)) if i > idx_cut])
            axis.plot(xdata, ydata, label=hfun.key_title(x, settings), linewidth=int(font / 6), markersize=font-6)
            
            "mean" 
            wH.append(wHnow)
            idx = min(range(len(xdata)), key=lambda i: abs(xdata[i] - wHnow));  LTA.append(ydata[idx])
            "typical"  
            wH_typ.append(wHtypnow)
            idx = min(range(len(xdata)), key=lambda i: abs(xdata[i] - wHtypnow));  val_at_typ.append(ydata[idx])
            
            #-- xy-ranges
            mini = ydata.min();  maxi = ydata.max();
            if mini < y_min and np.isfinite(mini): y_min = mini
            if maxi > y_max and np.isfinite(maxi): y_max = maxi
            mini = xdata.min();  maxi = xdata.max();
            if mini < x_min and np.isfinite(mini): x_min = mini
            if maxi > x_max and np.isfinite(maxi): x_max = maxi

    if normalise:
        ylab = "normalised\quad" + ylab
    hfun.set_plot_elements(axis = axis, xlim = (x_min, x_max), 
                                    ylim = (0.95*y_min, 1.05*y_max), ylabel = ylab, xlabel = xlab, settings=settings, font_size=font, set_legend=False)
  
    
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
    
    axis.plot(wH, LTA, linestyle='--', marker='o', color='black', linewidth=int(font / 6), markersize=font-4)
    axis.plot(wH_typ, val_at_typ, linestyle='--', marker='o', color='black', markerfacecolor='None', linewidth=int(font / 6), markersize=font-4)
   