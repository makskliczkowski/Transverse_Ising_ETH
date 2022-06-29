import importlib
from os import sep as kPSep
from numpy import array
from numpy import loadtxt
from numpy import exp
from numpy import sqrt
from numpy import log
from numpy import float as npfloat
#from scipy.special import binom as binomial
import scipy
import helper_functions as hfun
import config as cf
import copy
importlib.reload(cf)

#--- Global
user_settings = getattr(cf.plot_settings, 'settings')

#--- GENERAL
def load_taus():
    """
    Function to load Thouless times from file and return as numpy array
    """

    name = f"{cf.base_directory}ThoulessTime{kPSep}" + "_all" + ( ",p=%d,x=%d.dat"%(cf.p_sym, cf.x_sym) if cf.model else ",J0=%.2f,g0=%.2f.dat"%(cf.J0, cf.g0) )
    tau_data = loadtxt(name, unpack=True)
    return array(tau_data)

#--- compare all parameters except the scaling and vs one
def compare_params(tau_data, row):
    """
    """
    bool = 1
    for i in range(0, 5) :
        if i != user_settings['vs_idx']:
            if i == 4 and user_settings['vs_idx'] == 3 and cf.J0 == 0 and cf.g0 == 0:
                bool = bool and (abs(tau_data[4][row] - tau_data[3][row] / 2.) <= 2e-2)
            else:
                bool = bool and (abs(tau_data[i][row] - cf.params_arr[i]) <= 1e-10)
    return bool

#--- get tau data according to scaling in plot_settings
def get_tau_data(tau_data) : 
        vs_column = array(tau_data[user_settings['vs_idx']])
        taus = {}
        for i in range(0, len(vs_column)): 
            if(compare_params(tau_data, i)):
                par = vs_column[i]
                taus[f"%.5f"%(par)] = (tau_data[5][i] * (tau_data[6][i] if user_settings['physical_units'] else 1.0), tau_data[7][i])
        x_float = [];   tau = [];   gap = []
        if taus:
            lists = sorted(taus.items())
            x, data = zip(*lists)
            for j in range(0, len(x)) : 
                tau.append(data[j][0]); gap.append(data[j][1]); x_float.append(float(x[j]))
        return array(x_float), array(tau), array(gap)







#--- Function to Load data from file given by plot_settings
def load() :
    """
    Function to Load data from file given by plot_settings.
    
    Parameters used in function (all in config.py):
    -----------
    params_arr: array with model parameters set as follows:

    base_directory:    directory to main results (in which ThoulessTimes folder resides)

    plot_settings:  dictionary with plot settings, see config.py
    """
    print(user_settings)
    hfun.print_vars(cf.params_arr, cf.names)
    param_copy = cf.params_arr

    #--- SET SCALING RANGES AND DATA
    x0 = 0.2
    xend = 2.0
    dx = 0.1

    length = int((xend-x0) / dx) + 1
    #--- prepare scaling - axis
    vals = []
    if user_settings['scaling_idx'] == 0:
        if cf.hamiltonian: vals = range(10, 19, 2)
        else: vals = range(10, 17, 1)
    elif cf.model and user_settings['scaling_idx'] == 4:
        vals = range(0, cf.params_arr[0])
    else :
        for x in range(0, length) :
            vals.append(x0 + x * dx)
    vals = array(vals)

    #----- find data
    tau = []
    xvals = []
    gap_ratio = []

    tau_data = load_taus()
    new_vals = []
    for x in vals:
        cf.params_arr[user_settings['scaling_idx']] = x
        if user_settings['scaling_idx'] == 3 and cf.J0 == 0 and cf.g0 == 0:
            cf.params_arr[4] = int(100 * x / 2.) / 100.
        new_x, new_tau, new_gap = get_tau_data(tau_data)
        if new_tau.size > 1 :
            xvals.append(new_x)
            tau.append(new_tau)
            gap_ratio.append(new_gap)
            new_vals.append(x)

    vals = array(new_vals)

    xvals = array(xvals)
    tau = array(tau)
    gap_ratio = array(gap_ratio)
    
    #--- reset defaults
    cf.params_arr = param_copy
    return vals, xvals, tau, gap_ratio








#--- Function to plot thouless data given by plot_settings
def plot(axis1, axis2, new_settings = None) :
    """
    Plotter of Thouless times with plot_settings defining x-axis and scaling
    """
    global user_settings
    if new_settings != None:
        user_settings = new_settings

    def key_title(x):
        scaling_str = user_settings['scaling']
        if user_settings['scaling_idx'] == 2:
            scaling_str = hfun.var_name
        return scaling_str + (f"=%d"%(vals[i]) if user_settings['scaling_idx'] == 0 else f"=%.2f"%(vals[i]))

    #--- load data 
    vals, xvals, tau, gap_ratio = load()
    num_of_plots = len(tau)
    
    #--- plot first panel with thouless times
    marker_style = [];  face_colors = [];   ec = []
    y_min = 1.0e10;     y_max = -1.0e10;
    x_min = 1.0e10;     x_max = -1.0e10;

    rescale_by_L_nu = 0
    nu = 2
    for i in range(0, num_of_plots):
        yvals = tau[i]
        if rescale_by_L_nu and user_settings['vs_idx'] > 0 : yvals = yvals / (vals[i]**nu if user_settings['scaling_idx'] == 0 else cf.L**nu)
        yvals = cf.plot_settings.rescale(yvals, 'y')
        xx = cf.plot_settings.rescale(xvals[i], 'x')
        #if cf.hamiltonian and (user_settings['vs_idx'] == 0): xx = log(scipy.special.binom(xx, xx / 2)) / log(2)
        p = axis1.plot(xx, yvals, label=key_title(vals[i]))
        m = []; fc = [];    ec.append(p[0].get_color())
        
        #-- xy-ranges
        min = yvals.min();  max = yvals.max();
        if min < y_min: y_min = min
        if max > y_max: y_max = max
        min = xx.min();  max = xx.max();
        if min < x_min: x_min = min
        if max > x_max: x_max = max
        
        #--- plot markers with additional legend according to level spacing:
        # ~0.3865   -   filled squares
        # < 0.45    -   empty sqaures
        # > 0.45    -   empty circles
        # ~0.53     -   full circles
        for r in gap_ratio[i]: 
            m.append( 's' if r <= 0.46 else 'o')
            fc.append( p[0].get_color() if ( abs(r-0.53) <= 0.01 or abs(r-0.3865) <= 0.01 ) else 'none' )
        for j in range(0, len(tau[i])) :
            axis1.scatter(xx[j], yvals[j], edgecolors=ec[i], marker=m[j], s=50, facecolor=fc[j])
        
        # save markers for gap_ratio plot
        marker_style.append(m); face_colors.append(fc)

    #-- set panel1 details
    print(x_min, x_max, y_min, y_max, user_settings['vs_idx'], user_settings['scaling_idx'])
    yrange = (0.9*y_min, 1.1*y_max)
    ylab = "\\tau/L^{%d}"%nu if rescale_by_L_nu and user_settings['vs_idx'] > 0 else "\\tau"
    vs_str = user_settings['vs']
    if user_settings['vs_idx'] == 2:
        vs_str = hfun.var_name

    hfun.set_plot_elements(axis = axis1, xlim = (0.98*x_min, 1.02*x_max), 
                                ylim = yrange, ylabel = ylab, xlabel = vs_str, settings=user_settings)
    axis1.grid()
    axis1.legend()
    title = ""
    if (user_settings['vs_idx'] == 3 or user_settings['scaling_idx'] == 3) and cf.J0 == 0 and cf.g0 == 0:
        title = hfun.remove_info(hfun.info_param(cf.params_arr), user_settings['vs'], user_settings['scaling'], 'w') + ',w=0.5h'
    else :
        title = hfun.remove_info(hfun.info_param(cf.params_arr), user_settings['vs'], user_settings['scaling'])
    if user_settings['vs_idx'] != 2 :
        try : 
            title = list(title);    title[title.index('g')] = hfun.var_name;   title = "".join(title) # g
            title = list(title);    title[title.index('g')] = hfun.var_name;   title = "".join(title) # g0
        except ValueError:
                print("not found")
    axis1.title.set_text(title)




    rescale_by_L = 1
    #--- plot second panel with gap ratios
    x_min = 1.0e10;     x_max = -1.0e10;
    for i in range(0, num_of_plots):
        D = scipy.special.binom(vals[i], vals[i] / 2.)
        #norm = float(sqrt(D / sqrt(vals[i])) ) if (rescale_by_L and user_settings['scaling_idx'] == 0) else 1.0
        #norm = float(exp(log(2)/4*vals[i]) if (rescale_by_L and user_settings['scaling_idx'] == 0) else 1.0)
        norm = float(vals[i] if (rescale_by_L and user_settings['scaling_idx'] == 0) else 1.0)
        xpoints = xvals[i] * norm

        min = xpoints.min();  max = xpoints.max()
        if min < x_min: x_min = min
        if max > x_max: x_max = max
        axis2.plot(xpoints, gap_ratio[i], label=key_title(vals[i]))
        for j in range(0, len(tau[i])) :
            axis2.scatter(xpoints[j], gap_ratio[i][j], edgecolors=ec[i], marker=marker_style[i][j], s=50, facecolor=face_colors[i][j])
    new_set_class = copy.deepcopy(cf.plot_settings)
    new_set_class.set_x_rescale(rescale=0)
    new_set = getattr(new_set_class, 'settings')
    new_set['y_scale'] = 'linear';  new_set['x_scale'] = 'linear'
    xlab = new_set['vs'] + (" \\cdot L" if rescale_by_L else "")
    #xlab = new_set['vs'] + (" \\cdot e^{\\frac{ln2}{2}L}" if rescale_by_L else "")
    #xlab = new_set['vs'] + (" \\cdot D^{1/2}L^{1/4}" if rescale_by_L else "")
    hfun.set_plot_elements(axis = axis2, xlim = (0.98*x_min, 1.02*x_max), 
                                ylim = (0.37, 0.54), xlabel = xlab, ylabel = 'r', settings=new_set)
    #--- additional lines on plot
    axis2.axhline(y=0.5307, ls='--', color='black', label='GOE')
    axis2.axhline(y=0.3863, ls='--', color='red', label='Poisson')
    axis2.legend()
    axis2.title.set_text(title)