from os import sep as kPSep
from numpy import array
from numpy import loadtxt
from numpy import exp
from helper_functions import set_plot_elements
from config import *

#--- GENERAL
def load_taus():
    """
    Function to load Thouless times from file and return as numpy array
    """

    name = f"{base_directory}ThoulessTime{kPSep}" + "_all" + ( ",p=%d,x=%d.dat"%(p_sym, x_sym) if model else ",J0=%.2f,g0=%.2f.dat"%(J0, g0) )
    tau_data = loadtxt(name, unpack=True)
    return array(tau_data)

#--- compare all parameters except the scaling and vs one
def compare_params(tau_data, row):
    bool = 1
    for i in range(0, 5) :
        if i != plot_settings['vs_idx']:
            bool = bool and (abs(tau_data[i][row] - params_arr[i]) <= 1e-10)
    return bool
#--- get tau data according to scaling in plot_settings
def get_tau_data(tau_data) : 
        vs_column = array(tau_data[plot_settings['vs_idx']])
        taus = {}
        for i in range(0, len(vs_column)): 
            if(compare_params(tau_data, i)):
                par = vs_column[i]
                taus[f"%.5f"%(par)] = (tau_data[5][i] * (tau_data[6][i] if plot_settings['physical_units'] else 1.0), tau_data[7][i])
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
    Parameters:
    -----------
    params_arr: array with model parameters set as follows:
                 
            --- first input are the scaling parameters in this order
    
    base_directory:    directory to main results (in which ThoulessTimes folder resides)

    plot_settings:  dictionary with plot settings, see config.py
    """

    #--- SET SCALING RANGES AND DATA
    x0 = 0.2
    xend = 1.95
    dx = 0.1

    length = int((xend-x0) / dx) + 1
    #--- prepare scaling - axis
    vals = []
    if plot_settings['scaling_idx'] == 0:
        vals = range(12, 17)
    elif plot_settings['scaling_idx'] == 4:
        vals = range(0, params_arr[0])
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
        params_arr[plot_settings['scaling_idx']] = x
        new_x, new_tau, new_gap = get_tau_data(tau_data, params_arr, plot_settings)
        if new_tau.size > 1 :
            xvals.append(new_x)
            tau.append(new_tau)
            gap_ratio.append(new_gap)
            new_vals.append(x)

    vals = array(new_vals)

    xvals = array(xvals)
    tau = array(tau)
    gap_ratio = array(gap_ratio)
    
    return vals, xvals, tau, gap_ratio








#--- Function to plot thouless data given by plot_settings
def plot(axis1, axis2) :
    """
    Plotter of Thouless times with plot_settings defining x-axis and scaling
    """
    def key_title(x):
        return plot_settings['scaling'] + (f"=%d"%(vals[i]) if plot_settings['scaling_idx'] == 0 else f"=%.2f"%(vals[i]))
    def xform(x) :
        return x if plot_settings['rescale'] == 0 else 1. / x**plot_settings['nu']

    #--- load data 
    vals, xvals, tau, gap_ratio = load(params_arr, base_directory, plot_settings, model)
    num_of_plots = len(tau)
    
    #--- plot first panel with thouless times
    marker_style = [];  face_colors = [];   ec = []
    for i in range(0, num_of_plots):
        yvals = tau[i]
        if(plot_settings['scaling_idx'] == 0):   yvals = yvals / exp(0.01*vals[i]**2)
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
    yrange = (1e-2, 5e3) if plot_settings['physical_units'] else (1e-5, 1e0)
    set_plot_elements(axis1, [xform(xvals[0][0]), xform(xvals[0][len(xvals[0])-1])], yrange, 'tau', plot_settings)
    axis1.grid()
    axis1.legend()


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