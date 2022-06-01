import importlib
from os import sep as kPSep
from numpy import array
from numpy import loadtxt
from numpy import exp
import helper_functions as hfun
import config as cf
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
    xend = 1.95
    dx = 0.1

    length = int((xend-x0) / dx) + 1
    #--- prepare scaling - axis
    vals = []
    if user_settings['scaling_idx'] == 0:
        vals = range(12, 17)
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
        return user_settings['scaling'] + (f"=%d"%(vals[i]) if user_settings['scaling_idx'] == 0 else f"=%.2f"%(vals[i]))
    def xform(x) :
        return x if user_settings['rescaleX'] == 0 else 1. / x**user_settings['nu']

    #--- load data 
    vals, xvals, tau, gap_ratio = load()
    num_of_plots = len(tau)
    
    #--- plot first panel with thouless times
    marker_style = [];  face_colors = [];   ec = []
    for i in range(0, num_of_plots):
        yvals = tau[i]
        if(user_settings['scaling_idx'] == 0):   yvals = yvals / exp(0.01*vals[i]**2)
        p = axis1.plot(xform(xvals[i]), yvals, label=key_title(vals[i]))
        m = []; fc = [];    ec.append(p[0].get_color())
        
        #--- plot markers with additional legend according to level spacing:
        # ~0.3865   -   filled squares
        # < 0.45    -   empty sqaures
        # > 0.45    -   empty circles
        # ~0.53     -   full circles
        for r in gap_ratio[i]: 
            m.append( 's' if r <= 0.46 else 'o')
            fc.append( p[0].get_color() if ( abs(r-0.53) <= 0.01 or abs(r-0.3865) <= 0.02 ) else 'none' )
        for j in range(0, len(tau[i])) :
            axis1.scatter(xform(xvals[i][j]), yvals[j], edgecolors=ec[i], marker=m[j], s=50, facecolor=fc[j])
        
        # save markers for gap_ratio plot
        marker_style.append(m); face_colors.append(fc)

    #-- set panel1 details
    yrange = (8e-1, 1e3) if user_settings['physical_units'] else (1e-5, 1e0)
    xlab = r"$1\ /\ %s^{%.d}$"%(user_settings['vs'], user_settings['nu']) if user_settings['rescaleX'] else user_settings['vs']
    hfun.set_plot_elements(axis = axis1, xlim = [xform(xvals[0][0]), xform(xvals[0][len(xvals[0])-1])], ylim = yrange, xlabel = xlab, ylabel = 'tau', settings=user_settings)
    axis1.grid()
    axis1.legend()
    axis1.title.set_text(hfun.remove_info(hfun.info_param(cf.params_arr), user_settings['vs'], user_settings['scaling']))

    #--- plot second panel with gap ratios
    for i in range(0, num_of_plots):
        axis2.plot(xvals[i], gap_ratio[i], label=key_title(vals[i]))
        for j in range(0, len(tau[i])) :
            axis2.scatter(xvals[i][j], gap_ratio[i][j], edgecolors=ec[i], marker=marker_style[i][j], s=50, facecolor=face_colors[i][j])
    axis2.set_ylim(0.37, 0.54)
    #--- additional lines on plot
    axis2.axhline(y=0.5307, ls='--', color='black', label='GOE')
    axis2.axhline(y=0.3863, ls='--', color='red', label='Poisson')
    axis2.legend()
    axis2.title.set_text(hfun.info())