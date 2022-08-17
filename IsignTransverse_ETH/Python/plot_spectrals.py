from turtle import width
import numpy as np
import importlib
import costfun.costfun as cost
import helper_functions as hfun
import config as cf
import copy
importlib.reload(cf)
importlib.reload(cost)
importlib.reload(hfun)
import pandas as pd
from os import sep as kPSep
from os.path import exists

#--- Global
user_settings = getattr(cf.plot_settings, 'settings')
#--- SET SCALING RANGES AND DATA
x0 = 0.3
xend = 0.9
dx = 0.1
length = int((xend-x0) / dx) + 1

def get_scaling_array(settings = None):
    if settings == None:
        settings = user_settings
    vals = []
    if settings['scaling_idx'] == 0:
        if cf.hamiltonian: vals = range(12, 19, 2)
        else: vals = range(11, 17, 1)
    elif settings['scaling_idx'] == 5:
        vals = range(1, int(cf.params_arr[0] / 2) + 1)
    else :
        for x in range(0, length) :
            vals.append(x0 + x * dx)
    return np.array(vals)

def plot_spectral(axis, dir = "", settings = None, xlab = None, ylab = None, yscale=None, xscale=None, func_x=None, func_y=None, normalise=False, seper="\t"):
    #-- main settings
    if settings == None:
        settings = user_settings
    #-- labels
    if xlab == None: xlab = "x"
    if ylab == None: ylab = "y"
    #-- axis scales
    if xscale != None: settings['x_scale'] = xscale
    if yscale != None: settings['y_scale'] = yscale
    #-- axis rescaling
    if func_x == None: func_x = lambda x, a: x
    if func_y == None: func_y = lambda y, a: y
    param_copy = copy.deepcopy(cf.params_arr)

    #--- prepare scaling - axis
    vals = get_scaling_array(settings=settings)

    y_min = 1.0e10;     y_max = -1.0e10;
    x_min = 1.0e10;     x_max = -1.0e10;
    #--- load data and plot one-by-one
    wH = []
    LTA = []
    for x in vals:
        cf.params_arr[settings['scaling_idx']] = x
        if settings['scaling_idx'] == 3 and cf.J0 == 0 and cf.g0 == 0:
            cf.params_arr[4] = int(100 * x / 2.) / 100.
        filename = (hfun.info_param(cf.params_arr) if cf.hamiltonian else hfun.remove_info(hfun.info_param(cf.params_arr), 'J') + ".dat")

        if settings['scaling_idx'] == 5 and settings['operator'] < 8:
            filename = dir + ("j=%d%s" if settings['operator'] < 3 else "q=%d%s")%(x, kPSep) + cf.operator_names[settings['operator']] + "%d"%x + filename
        elif settings['scaling_idx'] == 0 and settings['site'] < 0:
            filename = dir + ("j=%d%s" if settings['operator'] < 3 else "q=%d%s")%(x / 2, kPSep) + cf.operator_names[settings['operator']] + "%d"%(x/2) + filename
        else :
            filename = dir + cf.subdir + cf.op_name + filename
        
        if exists(filename):
            data = pd.read_table(filename, sep=seper, header=None)
            xdata = func_x(data[0], x)
            ytmp = data[1]
            if normalise and seper == "\t":
                ytmp = (ytmp - data[3][0]) / np.abs(ytmp[len(ytmp)-1] - ytmp[0])
            ydata = func_y(ytmp, x)
            axis.plot(xdata, ydata, label=hfun.key_title(x, settings), linewidth=3, markersize=12)
           
            if seper == "\t":
                wHnow = func_x(data[2][0], x)
                wH.append(wHnow)
                idx = min(range(len(xdata)), key=lambda i: abs(xdata[i] - wHnow))
                LTA.append(ydata[idx])
            
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
                                    ylim = (0.95*y_min, 1.05*y_max), ylabel = ylab, xlabel = xlab, settings=settings, font_size=18, set_legend=False)
  
    
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
    axis.plot(wH, LTA, linestyle='--', marker='o', color='black', linewidth=3, markersize=12)
    #--- reset defaults
    cf.params_arr = param_copy