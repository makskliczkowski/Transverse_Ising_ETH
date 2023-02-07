import importlib
from os import sep as kPSep
import numpy as np

import costfun.costfun as cost
import utils.helper_functions as hfun
import config as cf
import copy
importlib.reload(cf)
importlib.reload(cost)
importlib.reload(hfun)

#--- Global
user_settings = getattr(cf.plot_settings, 'settings')

#--- GENERAL
def load_taus():
    """
    Function to load Thouless times from file and return as numpy array
    """

    name = f"{cf.base_directory}ThoulessTime{kPSep}" + "_all" + ( ",p=%d,x=%d.dat"%(cf.p_sym, cf.x_sym) if cf.model == 1 else ",J0=%.2f,g0=%.2f.dat"%(cf.J0, cf.g0) )
    tau_data = np.loadtxt(name, unpack=True)
    return np.array(tau_data)

#--- compare all parameters except the scaling and vs one
def compare_params(tau_data, row, settings):
    """
    """
    bool = 1
    for i in range(0, 5) :
        if i != settings['vs_idx']:
            if i == 4 and settings['vs_idx'] == 3 and cf.J0 == 0 and cf.g0 == 0:
                bool = bool and (abs(tau_data[4][row] - tau_data[3][row] / 2.) <= 2e-2)
            else:
                bool = bool and (abs(tau_data[i][row] - cf.params_arr[i]) <= 1e-10)
    return bool

#--- get tau data according to scaling in plot_settings
def get_tau_data_from_all_file(tau_data, settings) : 
        vs_column = np.array(tau_data[settings['vs_idx']])
        taus = {}
        for i in range(0, len(vs_column)): 
            if(compare_params(tau_data, i, settings)):
                local_pars = copy.deepcopy(cf.params_arr)
                par = vs_column[i]
                local_pars[settings['vs_idx']] = par
                if settings['vs_idx'] != 4 or par <= 1.0:
                    filename = cf.base_directory + "STATISTICS" + kPSep + "raw_data" + kPSep + hfun.info_param(local_pars)
                    wH = 0; stats = {}
                    try:
                        stats = hfun.load_stats(filename)
                        
                        wH = stats['mean level spacing']
                        wH = hfun.heisenberg_time( system_size = tau_data[0][i], dim = tau_data[9][i]) if np.isnan(wH) else 1. / wH
                    except FileNotFoundError:   wH = hfun.heisenberg_time( system_size = tau_data[0][i], dim = tau_data[9][i])
                    if np.isnan(wH):            wH = hfun.heisenberg_time( system_size = tau_data[0][i], dim = tau_data[9][i])
                    
                    taus[f"%.5f"%(par)] = (tau_data[5][i] * (wH if settings['physical_units'] else 1.0), tau_data[7][i])
        x_float = [];   tau = [];   gap = []
        if taus:
            lists = sorted(taus.items())
            x, data = zip(*lists)
            for j in range(0, len(x)) : 
                tau.append(data[j][0]); gap.append(data[j][1]); x_float.append(float(x[j]))
        return np.array(x_float), np.array(tau), np.array(gap)

#-- get thouless times from data chosen by settings['vs']
def get_tau_data(settings):
    dir = f"{cf.base_directory}ThoulessTime{kPSep}"
    xvals = np.array([]); taus = np.array([]);  gaps = np.array([])
    try:
        data = np.loadtxt(dir + hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs']) + ".dat", unpack=True)
        if len(data) > 0:
            xvals = np.array(data[0])
            taus = np.array(data[1]) * (np.array(data[2]) if settings['physical_units'] else 1.0)
            gaps = np.array(data[3])
    except OSError:
        xvals = np.array([]); taus = np.array([]);  gaps = np.array([])
    return xvals, taus, gaps

#--- Function to Load data from file given by plot_settings
def load(settings = None, vals = None) :
    """
    Function to Load data from file given by plot_settings.
    
    Parameters used in function (all in config.py):
    -----------
    params_arr: array with model parameters set as follows:

    base_directory:    directory to main results (in which ThoulessTimes folder resides)

    plot_settings:  dictionary with plot settings, see config.py
    """
    if settings == None:
        settings = user_settings
    #print(user_settings)
    #hfun.print_vars(cf.params_arr, cf.names)
    param_copy = cf.params_arr

    #--- SET SCALING RANGES AND DATA
    if vals is None:
        vals = hfun.get_scaling_array(settings=settings)
    #----- find data
    tau = []
    xvals = []
    gap_ratio = []
    
    tau_data = load_taus()
    new_vals = []
    #print(vals)
    for x in vals:
        cf.params_arr[settings['scaling_idx']] = x
        if settings['scaling_idx'] == 3 and cf.J0 == 0 and cf.g0 == 0:
            cf.params_arr[4] = int(100 * x / 2.) / 100.
        
        new_x, new_tau, new_gap = get_tau_data_from_all_file(tau_data, settings)
        if new_tau.size > 1 :
            new_vals.append(x)
            index = []
            for i in range(len(new_gap)):
                if new_gap[i] < 0.37: index.append(i)
            xvals.append(np.delete(new_x, index))
            tau.append(np.delete(new_tau, index))
            gap_ratio.append(np.delete(new_gap, index))

    #--- reset defaults
    cf.params_arr = param_copy
    importlib.reload(cf)

    return np.array(new_vals), np.array(xvals), np.array(tau), np.array(gap_ratio)


def replot_taus(axis, vals, xvals, tau, gap_ratio, settings = None,
                linewidth=0, fontsize=14, use_grid=False, markersize=30):
    #--- plot first panel with thouless times
    ec = []
    y_min = 1.0e10;     y_max = -1.0e10;
    x_min = 1.0e10;     x_max = -1.0e10;
    print(vals)
    num_of_plots = len(vals)
    for i in range(0, num_of_plots):
        yvals = cf.plot_settings.rescale(tau[i], 'y')
        xx = cf.plot_settings.rescale(xvals[i] , 'x')
        p = None
        if linewidth == 0:  p = axis.plot(xx, yvals, linewidth=linewidth)
        else:               p = axis.plot(xx, yvals, label=hfun.key_title(vals[i], settings), linewidth=linewidth)
        m = []; fc = [];    ec.append(p[0].get_color())
        
        #-- xy-ranges
        min = yvals.min();  max = yvals.max();
        if min < y_min and np.isfinite(min): y_min = min
        if max > y_max and np.isfinite(max): y_max = max
        min = xx.min();  max = xx.max();
        if min < x_min and np.isfinite(min): x_min = min
        if max > x_max and np.isfinite(max): x_max = max
        
        for r in gap_ratio[i]: 
            m.append('o')       # leaving behind in case of change in marker
            fc.append( p[0].get_color() if abs(r-0.53) <= 0.02 else 'none' )
        for j in range(0, len(tau[i])) :
            if linewidth == 0 and j == 0:   
                axis.scatter(xx[j], yvals[j], edgecolors=ec[i], marker=m[j], s=markersize, facecolor=fc[j], label=hfun.key_title(vals[i], settings))
            else: axis.scatter(xx[j], yvals[j], edgecolors=ec[i], marker=m[j], s=markersize, facecolor=fc[j])
        

    #-- set panel1 details
    yrange = (0.9*y_min, 1.1*y_max)
    ylab = "t_{Th}"
    vs_str = settings['vs']
    if settings['vs_idx'] == 2: vs_str = hfun.var_name
    elif settings['vs_idx'] == 4 and cf.model == 2: vs_str = "\\varepsilon"

    xlab = vs_str
    hfun.set_plot_elements(axis = axis, xlim = (x_min - np.sign(x_min) * 0.02*x_min, 1.02*x_max), 
                                ylim = yrange, ylabel = ylab, xlabel = xlab, settings=settings, font_size=fontsize, set_legend=True)
    if use_grid: axis.grid()
    axis.tick_params(axis="both",which='major',direction="in",length=6)
    axis.tick_params(axis="both",which='minor',direction="in",length=3)

    return ec               
#--- Function to plot thouless data given by plot_settings

def plot_taus(axis, settings = None, vals = None, 
                linewidth=0, fontsize=14, markersize=30, return_data=False):
    """
    Plotter of Thouless times with plot_settings defining x-axis and scaling
    
    Parameters:
    --------------
        axis: plt.figure
            axis to plot data on

        new_settings: dict
            settings defining plot behavior and data extraction
    """
    if settings == None:
        settings = user_settings

    #--- load data 
    vals, xvals, tau, gap_ratio = load(settings, vals)

    ec = replot_taus(axis=axis,
                    vals=vals,
                    xvals=xvals,
                    tau=tau,
                    gap_ratio=gap_ratio,
                    settings=settings,
                    linewidth=linewidth,
                    fontsize=fontsize,
                    markersize=markersize)
    if return_data: return vals, xvals, tau, gap_ratio, ec
    #title = ""
    #if (settings['vs_idx'] == 3 or settings['scaling_idx'] == 3) and cf.J0 == 0 and cf.g0 == 0:
    #    title = hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs'], settings['scaling'], 'w') + ',w=0.5h'
    #else :
    #    title = hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs'], settings['scaling'])
    #
    #if settings['vs_idx'] != 2 :
    #    try : 
    #        title = list(title);    idx=title.index('g');   title.remove('g')
    #        dummy=list(hfun.var_name)
    #        for i in range(0,len(dummy)):
    #            title.insert(i + idx, dummy[i]) 
    #        title = "".join(title) # g
    #        #title = list(title);    title[title.index('g')] = hfun.var_name;   title = "".join(title) # g0
    #    except ValueError:
    #            print("not found")
    #axis.title.set_text(r"$%s$"%title)



#--- Function to plot thouless data given by plot_settings with gap ratio on second axis
def plot_with_gap_ratio(axis1, axis2, new_settings = None, use_scaling_ansatz = 0, scaling_ansatz=None, crit_fun=None) :
    """
    Plotter of Thouless times with plot_settings defining x-axis and scaling
    """
    if new_settings == None:
        new_settings = user_settings

    #--- load data 
    vals, xvals, tau, gap_ratio = load(new_settings)
    
    num_of_plots = len(tau)

    rescale_by_L_nu = 0

    #cost functinon
    if scaling_ansatz is None: scaling_ansatz='classic'
    if crit_fun is None: crit_fun='const'
    rescale_fun = cost.resc_functions_dict[scaling_ansatz]
    critical_fun = cost.crit_functions_dict[crit_fun]
    params = [0,0,0,0,0,0]
    par = params[0]
    crit_pars = params[1:]
    if new_settings['scaling_idx'] == 0 and use_scaling_ansatz:
        par, crit_pars, cost_fun, status = cost.get_crit_points(x=xvals,
                                                                y=gap_ratio,
                                                                vals=vals,
                                                                scaling_ansatz=scaling_ansatz,
                                                                crit_fun=crit_fun)
        
        print(scaling_ansatz + ": ", cost_fun, [critical_fun(vals[i], *crit_pars) for i in range(len(vals))])
    
    print(par)

    #--- plot first panel with thouless times
    marker_style = [];  face_colors = [];   ec = []
    y_min = 1.0e10;     y_max = -1.0e10;
    x_min = 1.0e10;     x_max = -1.0e10;

    nu = 2
    for i in range(0, num_of_plots):
        yvals = tau[i]
        if rescale_by_L_nu and new_settings['vs_idx'] > 0 : yvals = yvals / (vals[i]**nu if new_settings['scaling_idx'] == 0 else cf.L**nu)
        yvals = cf.plot_settings.rescale(yvals, 'y')
        temp = xvals[i] - critical_fun(vals[i], *crit_pars) if new_settings['scaling_idx'] == 0 and use_scaling_ansatz else xvals[i] 
        #scaling_ansatz(x=xvals[i], par_crit=crit_par, L=vals[i], par=mu) if use_scaling_ansatz else xvals[i]
        xx = cf.plot_settings.rescale(temp if new_settings['scaling_idx'] == 0 else xvals[i] , 'x')
        #if cf.hamiltonian and (new_settings['vs_idx'] == 0): xx = log(scipy.special.binom(xx, xx / 2)) / log(2)
        #yvals = np.log(yvals) / np.log(xx)
        p = axis1.plot(xx, yvals, label=hfun.key_title(vals[i], new_settings))
        m = []; fc = [];    ec.append(p[0].get_color())
        
        #-- xy-ranges
        min_val = yvals.min();  max_val = yvals.max();
        if min_val < y_min and np.isfinite(min_val): y_min = min_val
        if max_val > y_max and np.isfinite(max_val): y_max = max_val
        min_val = xx.min();  max_val = xx.max();
        if min_val < x_min and np.isfinite(min_val): x_min = min_val
        if max_val > x_max and np.isfinite(max_val): x_max = max_val
        
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
    yrange = (0.9*y_min, 1.1*y_max)
    ylab = "\\tau/L^{%d}"%nu if rescale_by_L_nu and new_settings['vs_idx'] > 0 else "\\tau"
    vs_str = new_settings['vs']
    if new_settings['vs_idx'] == 2:
        vs_str = hfun.var_name

    xlab = (vs_str + (" - " + vs_str + "_c") if use_scaling_ansatz and (new_settings['scaling_idx'] == 0) else vs_str)
    hfun.set_plot_elements(axis = axis1, xlim = (x_min - np.sign(x_min) * 0.02*x_min, 1.02*x_max), 
                                ylim = yrange, ylabel = ylab, xlabel = xlab, settings=new_settings)
    axis1.grid()
    axis1.legend()
    title = ""
    if (new_settings['vs_idx'] == 3 or new_settings['scaling_idx'] == 3) and cf.J0 == 0 and cf.g0 == 0:
        title = hfun.remove_info(hfun.info_param(cf.params_arr), new_settings['vs'], new_settings['scaling'], 'w') + ',w=0.5h'
    else :
        title = hfun.remove_info(hfun.info_param(cf.params_arr), new_settings['vs'], new_settings['scaling'])
    if new_settings['vs_idx'] != 2 :
        try : 
            title = list(title);
            idx=title.index('g')    
            title.remove('g')
            dummy=list(hfun.var_name)
            for i in range(0,len(dummy)):
                title.insert(i + idx, dummy[i]) 
            title = "".join(title) # g
            #title = list(title);    title[title.index('g')] = hfun.var_name;   title = "".join(title) # g0
        except ValueError:
                print("not found")
    axis1.title.set_text(r"$%s$"%title)



    xlab = ("(" + vs_str + (" - " + vs_str + "_c) \\cdot L^{\\nu}") if use_scaling_ansatz and (new_settings['scaling_idx'] == 0) else vs_str)
    xlab = ("(" + vs_str + (" - " + vs_str + "_c)^{\\nu} \\cdot e^{\\frac{ln2}{2}L}") if use_scaling_ansatz and (new_settings['scaling_idx'] == 0) else vs_str)
    xlab = cost.scale_ansatz_label[scaling_ansatz](vs_str)
    #--- plot second panel with gap ratios
    x_min = 1.0e10;     x_max = -1.0e10;
    for i in range(0, num_of_plots):
        xpoints = rescale_fun(xvals[i], vals[i], critical_fun, 
                    par, *crit_pars) if new_settings['scaling_idx'] == 0 and use_scaling_ansatz else xvals[i] 

        min_val = xpoints.min();  max_val = xpoints.max()
        if np.isfinite(min_val) and min_val < x_min: x_min = min_val
        if np.isfinite(max_val) and max_val > x_max: x_max = max_val
        #axis2.plot(xpoints, gap_ratio[i], label=key_title(vals[i]))
        for j in range(0, len(tau[i])) :
            axis2.scatter(xpoints[j], gap_ratio[i][j], edgecolors=ec[i], marker=marker_style[i][j], s=50, facecolor=face_colors[i][j])
    new_set_class = copy.deepcopy(cf.plot_settings)
    new_set_class.set_x_rescale(rescale=0)
    new_set = getattr(new_set_class, 'settings')
    new_set['y_scale'] = 'linear'; new_set['x_scale'] = 'linear'
    
    print(x_min, x_max, y_min, y_max, new_settings['vs_idx'], new_settings['scaling_idx'])
    hfun.set_plot_elements(axis = axis2, xlim = (0.98*x_min, 1.02*x_max), 
                                ylim = (0.37, 0.54), xlabel = xlab, ylabel = 'r', settings=new_set, set_legend=False)
    if new_set['scaling_idx'] == 0 and use_scaling_ansatz:
        axis2.annotate(r"$\nu=%.2f$"%par, xy=(150, 260), xycoords='axes points', color='black', size=20)
        for i in range(0, len(crit_pars)):
            axis2.annotate(r"$x_{%d}=%.4f$"%(i,crit_pars[i]), xy=(150, 240 - 20*i), xycoords='axes points', color='black', size=20)  
    
    #--- additional lines on plot
    axis2.axhline(y=0.5307, ls='--', color='black', label='GOE')
    axis2.axhline(y=0.3863, ls='--', color='red', label='Poisson')
    axis2.legend()
    axis2.title.set_text(r"$%s$"%title)
    axis2.set_xscale('linear')
    crit = [critical_fun(vals[i], *crit_pars) for i in range(len(vals))]
    axis2.set_xlim(min(crit) - 30, max(crit) + 30)