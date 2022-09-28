from turtle import width
import numpy as np
import importlib
import utils.helper_functions as hfun
import config as cf
import copy
import thouless_times as thouless
importlib.reload(cf)
importlib.reload(hfun)
importlib.reload(thouless)
import pandas as pd
from os import sep as kPSep
from os.path import exists
from utils.fit_functions import *
from scipy.optimize import curve_fit as fit


#--- Global
user_settings = getattr(cf.plot_settings, 'settings')

time_dir = cf.base_directory + "timeEvolution%s"%kPSep
int_dir = cf.base_directory + "IntegratedResponseFunction%s"%kPSep
spec_dir = cf.base_directory + "ResponseFunction%s"%kPSep

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
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
        filename = dir + cf.subdir(operator, parameter, settings['smoothed']) + cf.operator_names[operator] + "%d"%parameter + filename
    elif settings['scaling_idx'] == 0 and site < 0 and operator < 8:
        filename = dir + cf.subdir(operator, parameter / 2, settings['smoothed']) + cf.operator_names[operator] + "%d"%(parameter/2) + filename
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

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

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
    if spec == "time":      dir = cf.base_directory + "timeEvolution%s"%kPSep
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
        vals = hfun.get_scaling_array(settings=settings)

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
   
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
def get_thouless_time(par, set_class = None, vals = None):
    if set_class is None:   set_class = cf.plot_settings    
    if par is None:         
        print("No parameter input!")
        raise ValueError

    " Find thouless data "
    tau_data = []
    status_time = False
    try :
        tau_data = thouless.load(getattr(set_class, 'settings'))
        status_time = True
    except OSError:
        print("No Thouless data present")
        
    taus = [];  gap_ratio = []; xvalues = []
    if status_time:    
        idx = list(tau_data[0]).index(par)
        if vals is None:    vals = tau_data[1][idx]
        for x in vals:
            idx2 = hfun.find_index(tau_data[1][idx], x)
            if idx2 >= 0:   
                taus.append(tau_data[2][idx][idx2])
                gap_ratio.append(tau_data[3][idx][idx2])
            else:           
                taus.append(np.nan)
                gap_ratio.append(np.nan)
    xvalues = vals
    return status_time, np.array(xvalues), np.array(taus), np.array(gap_ratio)

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
def get_relax_times(vals = None, set_class = None, operator = -1, site = -2):
    """ 
    Find relaxation times from both integrated spectral function at I(w)=1/2 and fitting autocorrelation function
    """
    if set_class is None:   set_class = cf.plot_settings 
    settings = getattr(set_class, 'settings')

    if operator < 0: operator = settings['operator']
    if site < -1: site = settings['site']
    if vals is None:    vals = hfun.get_scaling_array()

   
    set_class_th = copy.deepcopy(set_class)
    set_class_th.set_vs(settings['scaling'])
    set_class_th.set_scaling('L')

    status_time, xvals, taus, gap_ratio = get_thouless_time(par = cf.L, set_class=set_class_th, vals=vals)
    relax_time = []
    
    relaxt_time_fit = [];   tH = [];    tH_typ = [];
    
    for i in range(0, len(vals)):
        x = vals[i]

        "Find relax time from integrated spec fun"
        status, xdata, ydata, wHnow, wHtypnow = load_spectral(dir=int_dir, 
                                                                        settings=settings, 
                                                                        parameter=x,
                                                                        spec="int",
                                                                        normalise=True,
                                                                        operator=operator,
                                                                        site=site
                                                                        )
        if status:
            idx = min(range(len(ydata)), key=lambda i: abs(ydata[i] - 0.5));
            relax_time.append(1. / xdata[idx])
            tH.append(1. / wHnow)
            tH_typ.append(1. / wHtypnow)
        else:
            relax_time.append(np.nan)
            tH.append(np.nan)
            tH_typ.append(np.nan)
    
        "Find relax time from autocorrelation function"
        status2, xdata2, ydata2, wHnow, wHtypnow = load_spectral(dir=time_dir, 
                                                                        settings=settings, 
                                                                        parameter=x,
                                                                        spec="time",
                                                                        normalise=True,
                                                                        operator=operator,
                                                                        site=site
                                                                        )
        
        if status2:
            cut = 30
            if x <= 0.2: cut = 120
            xfull = xdata2
            xdata2 = np.array([xdata2[i] for i in range(0,len(xdata2)) if (xdata2[i] < 5000 and xdata2[i] > cut)])
            ydata2 = np.array([ydata2[i] for i in range(0,len(ydata2)) if (xfull[i] < 5000 and xfull[i] > cut)])
            
            ydata2 = np.log10(np.abs(ydata2))
            idx_zero = np.argmin((ydata2))
            ydata2 = ydata2[:idx_zero - 15]
            xdata2 = xdata2[:idx_zero - 15]
            #print(pars)
            
            pars, sth = fit(f=lin_fit, 
                                xdata=xdata2, 
                                ydata=ydata2)
            relaxt_time_fit.append(pars[0])
        else:
            relaxt_time_fit.append(np.nan)

    return status_time, np.array(taus), np.array(relax_time), np.array(relaxt_time_fit), np.array(tH), np.array(tH_typ), np.array(gap_ratio)

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
def set_inset_style(axis, vals, settings, ylim = None, ylabel = None):
    """ 
    Sets style of plot with relaxation times
    """
    if ylim is None:    ylim = (1e-1, 7e3)
    if ylabel is None:  ylabel = "\\tau_{rel}"
    
    ii = settings['scaling_idx']
    xlab = "q/\\pi" if ii == 5 else (hfun.var_name if ii == 2 else settings['scaling'])

    hfun.set_plot_elements(axis = axis, xlim = (0.95*min(vals), 1.05*max(vals)), 
                                        ylim = ylim, ylabel = ylabel, xlabel = xlab, 
                                        settings=settings, font_size=16, set_legend=True)
    
    axis.set_yscale('log')
    axis.set_xscale('log')
    axis.tick_params(axis='both', which='both',length=2)

    axis.grid(b=True, which='major', color='0.75', linestyle='-')
    axis.grid(b=True, which='minor', color='0.85', linestyle='--')
    axis.legend(ncol=3, loc='lower center')

def set_inset(axis, settings, vals, taus, relax_time, relaxt_time_fit, tH, tH_typ, status_time = True):
    """ 
    Plots relaxation times
    """
    ii = settings['scaling_idx']
    xlab = "q/\\pi" if ii == 5 else (hfun.var_name if ii == 2 else settings['scaling'])
    if settings is None:
        settings = user_settings

    axis.plot(vals, relax_time, marker='o', label='int-fit')
    axis.plot(vals, relaxt_time_fit, marker='o', label='exp fit')
    axis.plot(vals, tH, linestyle='--', label=r"$t_H$", color='gray')
    #axis.plot(vals, tH_typ, linestyle=':', label=r"$t_H^{typ}$", color='gray')

    axis.plot(vals, 2e0 / vals**2, linestyle='--', color='red', label=r"$%s^{-2}$"%xlab)
    axis.plot(vals, 1e0 / vals**1., linestyle='--', color='black', label=r"$%s^{-1}$"%xlab)

    if status_time and settings['scaling_idx'] == 5: 
        axis.axhline(y=taus[0], ls='--', color='black')
        axis.annotate("$\it{Thouless}\ \it{Time}$", xy=(0.3,1.5e3), color='black', size=12)
    else: 
        if status_time == False or cf.model == 2: print('No data')
        else: axis.plot(vals, taus, linestyle='--', color='black', marker='o', label=r"$\tau_{Th}$")

    set_inset_style(axis, vals, settings)

#axis2.xaxis.set_minor_locator(plt.MultipleLocator(0.2))
#def format_func(value, tick_number):
#    return "%.1f"%value
#def format_func2(value, tick_number):
#    return "%d"%value
#axis2.xaxis.set_minor_formatter(plt.FuncFormatter(format_func))
#axis2.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
#axis2.yaxis.set_major_formatter(plt.FuncFormatter(format_func2))

