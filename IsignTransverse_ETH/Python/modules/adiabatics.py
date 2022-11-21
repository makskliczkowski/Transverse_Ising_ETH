import numpy as np
import importlib
import utils.helper_functions as hfun
import config as cf
import copy
importlib.reload(hfun)
importlib.reload(cf)
import pandas as pd
from os import sep as kPSep
from os.path import exists
from scipy.special import binom

def plot_agp(axis=None, settings_class = None, 
                which=1, operator = -1, site = -3):

    if which < 1 or which > 4: 
        print("Parameter 'which' entered with illegal value")
        return
    #-- main settings
    if settings_class == None:
        settings_class = cf.plot_settings
    settings = getattr(settings_class, 'settings')
    dir = cf.base_directory + "AGP" + kPSep
    param_copy = copy.deepcopy(cf.params_arr)

    #--- prepare scaling - axis
    vals = hfun.get_scaling_array(settings=settings)

    y_min = 1.0e10;     y_max = -1.0e10;
    x_min = 1.0e10;     x_max = -1.0e10;
    #--- load data and plot one-by-one
    print(vals)
    for x in vals:
        cf.params_arr[settings['scaling_idx']] = x
        if settings['scaling_idx'] == 3 and cf.J0 == 0 and cf.g0 == 0:
            cf.params_arr[4] = int(100 * x / 2.) / 100.
        filename = (hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs']) + ".dat" if cf.hamiltonian else hfun.remove_info(hfun.info_param(cf.params_arr), 'J', settings['vs']) + ".dat")

        if settings['scaling_idx'] == 5 and operator < 8:
            filename = dir + cf.operator_names[operator] + "%d"%x + kPSep + filename
        elif settings['scaling_idx'] == 0 and site < 0  and operator < 8:
            filename = dir + cf.operator_names[operator] + "%d"%(x/2) + kPSep + filename
        else :
            filename = dir + cf.operator_name(operator, site) + kPSep + filename
        
        filename2 = cf.base_directory + "STATISTICS" + kPSep + hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs']) + ".dat"
        if exists(filename):
            data = pd.read_table(filename, sep="\t", header=None)
            if "nan" in data[1][1]: continue
            stats = hfun.read_python_saved_dat_file(filename2)
            xdata = (np.array(data[0][1:])).astype(np.float)
            ydata = (np.array(data[which][1:])).astype(np.float)
            wH = stats[5]
            if which == 2:
                ydata = ydata * np.power(np.sqrt(x) / (binom(x, x/2)), 2.0) / x * (binom(x, x/2))
            elif which == 1:
                ydata = ydata / (binom(x, x/2))
            axis.plot(xdata, ydata, label=hfun.key_title(x, settings), marker='o')
            
            #-- xy-ranges
            mini = ydata.min();  maxi = ydata.max();
            if mini < y_min and np.isfinite(mini): y_min = mini
            if maxi > y_max and np.isfinite(maxi): y_max = maxi
            mini = xdata.min();  maxi = xdata.max();
            if mini < x_min and np.isfinite(mini): x_min = mini
            if maxi > x_max and np.isfinite(maxi): x_max = maxi
    ylab = ""
    if which == 1:
        ylab = "||\\mathcal{A}(A)||^2 / D"
    elif which == 2:
        ylab = "D\\cdot\\omega_H^2\\cdot\\chi^{typ}(A) / L"
    elif which == 3:
        ylab = "\\chi(A)"
    else :
      ylab = "||A||^2_{diag}"  
    hfun.set_plot_elements(axis = axis, xlim = (x_min, x_max), 
                                    ylim = (0.95*y_min, 1.05*y_max), ylabel = ylab, xlabel = settings['vs'], settings=settings, font_size=8)
  
    
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
    #---Thouless times
    
    #--- reset defaults
    cf.params_arr = param_copy