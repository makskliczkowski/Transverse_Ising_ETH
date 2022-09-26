import importlib
from os import sep as kPSep
import numpy as np

import costfun.costfun as cost
import helper_functions as hfun
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
def get_tau_data(tau_data, settings) : 
        vs_column = np.array(tau_data[settings['vs_idx']])
        taus = {}
        for i in range(0, len(vs_column)): 
            if(compare_params(tau_data, i, settings)):
                par = vs_column[i]
                #if par <= 1.0:
                taus[f"%.5f"%(par)] = (tau_data[5][i] * (tau_data[6][i] if settings['physical_units'] else 1.0), tau_data[7][i])
        x_float = [];   tau = [];   gap = []
        if taus:
            lists = sorted(taus.items())
            x, data = zip(*lists)
            for j in range(0, len(x)) : 
                tau.append(data[j][0]); gap.append(data[j][1]); x_float.append(float(x[j]))
        return np.array(x_float), np.array(tau), np.array(gap)



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
    print(vals)
    for x in vals:
        cf.params_arr[settings['scaling_idx']] = x
        if settings['scaling_idx'] == 3 and cf.J0 == 0 and cf.g0 == 0:
            cf.params_arr[4] = int(100 * x / 2.) / 100.
        new_x, new_tau, new_gap = get_tau_data(tau_data, settings)
        if new_tau.size > 1 :
            xvals.append(np.array(new_x))
            tau.append(np.array(new_tau))
            gap_ratio.append(np.array(new_gap))
            new_vals.append(np.array(x))

    #--- reset defaults
    cf.params_arr = param_copy
    importlib.reload(cf)
    return np.array(new_vals), np.array(xvals), np.array(tau), np.array(gap_ratio)



#--- Function to plot thouless data given by plot_settings
def plot(axis1, axis2, new_settings = None, use_scaling_ansatz = 0, scaling_ansatz=None, crit_fun=None) :
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
    params = [0,0,0,0,0,0]
    x_max=None
    for x in xvals: 
        for _x_ in x: 
            if x_max is None or _x_ > x_max: x_max = _x_
    if new_settings['scaling_idx'] == 0 and use_scaling_ansatz:
        x_min = 0 if crit_fun == 'free' else -x_max
        bounds = [ (0.0, 5.0), (x_min, x_max)]
        num_of_param = 0
        if crit_fun == 'free': num_of_param = len(vals) - 1
        elif crit_fun == 'power_law': num_of_param = 3
        elif crit_fun == 'const': num_of_param = 0
        else: num_of_param = 3
        for i in range(num_of_param): bounds.append((x_min, x_max))
        params, cost_fun = cost.cost_func_minization(x=xvals, y=gap_ratio, sizes=vals, 
                                        scale_func=scaling_ansatz, 
                                        crit_func=crit_fun,
                                        bnds=bounds,
                                        population_size=1e2,
                                        maxiterarions=1e2, workers=10, realisations=1
                                    )
        print(scaling_ansatz + ": ", cost_fun)
    par = params[0]
    crit_pars = np.array(params[1:])
    rescale_fun = cost.resc_functions_dict[scaling_ansatz]
    critical_fun = cost.crit_functions_dict[crit_fun]

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
        min = yvals.min();  max = yvals.max();
        if min < y_min and np.isfinite(min): y_min = min
        if max > y_max and np.isfinite(max): y_max = max
        min = xx.min();  max = xx.max();
        if min < x_min and np.isfinite(min): x_min = min
        if max > x_max and np.isfinite(max): x_max = max
        
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

        min = xpoints.min();  max = xpoints.max()
        if np.isfinite(min) and min < x_min: x_min = min
        if np.isfinite(max) and max > x_max: x_max = max
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
        for i in range(0, len(params)-1):
            axis2.annotate(r"$x_{%d}=%.4f$"%(i,params[i+1]), xy=(150, 240 - 20*i), xycoords='axes points', color='black', size=20)  
    
    #--- additional lines on plot
    axis2.axhline(y=0.5307, ls='--', color='black', label='GOE')
    axis2.axhline(y=0.3863, ls='--', color='red', label='Poisson')
    axis2.legend()
    axis2.title.set_text(r"$%s$"%title)