import importlib
from multiprocessing.sharedctypes import Value
import config as cf
from matplotlib.markers import MarkerStyle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import ceil
importlib.reload(cf)
import re

user_settings = getattr(cf.plot_settings, 'settings')
var_name = "\\Delta" if cf.hamiltonian == 1 else ("\\alpha" if cf.hamiltonian == 2 else "g")

def heisenberg_time(system_size, dim):
    chi = 0.341345
    par_term = np.sqrt(cf.J * cf.J + cf.h * cf.h + cf.g * cf.g + (cf.w * cf.w + cf.g0 * cf.g0 + cf.J0 * cf.J0) / 3.)
    if cf.hamiltonian == 1:
        par_term = np.sqrt( cf.J * cf.J / 8. + cf.h * cf.h + cf.g * cf.g / 16. + (cf.w * cf.w + cf.g0 * cf.g0 + cf.J0 * cf.J0) / 12.)
    return (chi * dim) / ( system_size**(0.5) * par_term)

def order_of_magnitude(a_value):
    #return 2
    if np.abs(a_value) < 1.0 and a_value != 0:
        m = np.abs(np.log10(np.abs(a_value)))
        return int(max(ceil(m) + 1., 2.))
    else: 
        return 2

def findOccurrences(s, ch):
    """
    Find all occurences of character ch in string s
    """
    return [i for i, letter in enumerate(s) if letter == ch]
    
#-------------------------- SET INFO
def info_sym(L, J, g, h, k, p, x, use_log_data =True):
    arr = [J, g, h]
    names = ['J', 'g', 'h']
    info = "_L=%d"%L
    for i, var in enumerate(arr):
        n = order_of_magnitude(var) if use_log_data else 2
        info += str(",%s={:.%df}"%(names[i], n)).format(round(var, n))
    return info + ",k=%d,p=%d,x=%d.dat"%(k,1 if p == 1 else -1,1 if x == 1 else -1)
    
def info_dis(L, J, J0, g, g0, h, w, use_log_data = True):
    arr = [J, J0, g, g0, h, w]
    names = ['J', 'J0', 'g', 'g0', 'h', 'w']
    info = "_L=%d"%L
    for i, var in enumerate(arr):
        n = order_of_magnitude(var) if use_log_data else 2
        info += str(",%s={:.%df}"%(names[i], n)).format(round(var, n))
    return info + ".dat"

def info(_L = cf.params_arr[0], _J = cf.params_arr[1], _J0 = cf.params_arr[8], 
            _g = cf.params_arr[2], _g0 = cf.params_arr[9], 
            _h = cf.params_arr[3], _w = cf.params_arr[4], 
          _k = cf.params_arr[5], _p = cf.params_arr[6], _x = cf.params_arr[7], use_log_data = True
    ):
    if(cf.model == 1) :
        return info_sym(_L, _J, _g, _h, _k, _p, _x, use_log_data)
    else :
        return info_dis(_L, _J, _J0, _g, _g0, _h, _w, use_log_data)
def info_param(params, use_log_data = True):
    return info(_L = params[0], _J = params[1], _J0 = params[8], 
            _g = params[2], _g0 = params[9], 
            _h = params[3], _w = params[4], 
          _k = params[5], _p = params[6], _x = params[7], use_log_data=use_log_data
    )

def get_params_from_info(info):
    """
    Exctract parameters from info and suplement rest from params array
    """
    found = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", info)
    found = [int(x) if x.isdigit() else float(x) for x in found]
    result = list(range(len(cf.params_arr)))
    
    if cf.model != 1:
        found = list(np.delete(found, [3, 6]))
        initial_found = found
        
        J0 = found[2];  g0 = found[4]
        result[0] = found[0]
        result[1:4] = found[1:7:2]
        result[4] = found[6]
        # symmetries
        result[5:8] = cf.params_arr[5:8]
        # rest of disordeer
        result[8:10] = [J0, g0]
        if len(initial_found) > 7: result.append(initial_found[-1])
    else:
        result[0:4] = found[:4]
        result[4] = cf.params_arr[4]
        # symmetries
        result[5:8] = found[4:7]
        # rest of disordeer
        result[8:10] = [cf.params_arr[8], cf.params_arr[9]]
        if len(initial_found) > 7: result.append(initial_found[-1])

    return result

def remove_info(info_str, *args):
    for x in args:
        try:
            idx = info_str.index(x)
            idx_com = (info_str[idx : :]).index(',') if (x != 'w' and x != 'x') else  (info_str[idx : :]).index('.dat')
            info_str = info_str[0 : idx - 1 :] + info_str[idx + idx_com : :]
        except ValueError:
            continue
    info_str = info_str[0 : info_str.index('.dat') :]
    return info_str
    
#-------------------------- PRETTY PRINT
def print_vars(arr_vals, names):
    for i in range(0,len(arr_vals)) :
        output = f"%s=%.2f"%(names[i], arr_vals[i]) if (isinstance(arr_vals[i], float)) else f"%s=%d"%(names[i], arr_vals[i])
        print(output);

#-------------------------- INDICES FUN
def find_index(data, value) :
    try:
        idx = list(data).index(value)
        return idx
    except ValueError:
        return -1

def key_title(x, settings):
    scaling_str = settings['scaling']
    if settings['scaling_idx'] == 2:
        scaling_str = var_name
    if settings['scaling_idx'] == 4 and cf.model == 2:
        scaling_str = "\\varepsilon"
    n = order_of_magnitude(x)
    return r"$" + (scaling_str + (f"=%d"%(x) if settings['scaling_idx'] == 0 or settings['scaling_idx'] == 5 else str("={:.%df}"%(n)).format(round(x, n)))) + "$"

#-------------------------- PLOT FANCY
#--------- scatter plot with differeent markertypes as list
def mscatter(x,y,ax=None, m=None, fc=None, **kw):
    if not ax: ax=plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            marker_obj = marker if isinstance(marker, MarkerStyle) else MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                        marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc

#--------- Set user settings on plot
def set_plot_elements(axis, xlim =[], ylim=[], xlabel = None, ylabel = None, settings = None, set_legend = True, font_size = 10):
    if settings == None:
        settings = user_settings

    if xlabel is not None:
        xlab = list(settings['func_x_name'])
        x_idx = xlab.index('Q')
        xlabel = list(xlabel)
        xlab.remove('Q')
        for i in range(0,len(xlabel)):
            xlab.insert(i + x_idx, xlabel[i]) 
        axis.set_xlabel("".join(xlab), rotation=0, fontsize=font_size+2, labelpad=font_size-8)
    if ylabel is not None:
        ylab = list(settings['func_y_name'])
        ylab[ylab.index('Q')] = ylabel
        axis.set_ylabel("".join(ylab), rotation=90, fontsize=font_size+2)
    
    axis.set_yscale(settings['y_scale'])
    axis.set_xscale(settings['x_scale'])
    axis.tick_params(axis='both', which='major', direction="in",length=6, labelsize=font_size)#, length=font_size-4, width=0.05*font_size)
    axis.tick_params(axis='both', which='minor', direction="in",length=6, labelsize=font_size)#, length=0.2*(font_size-4), width=0.05*font_size)
    
    if set_legend:
        axis.legend(frameon=False
                , loc='best'
                , fontsize=font_size)
    axis.set_axisbelow(True)
    x1, x2 = xlim
    y1, y2 = ylim
    if x1 != None and x2 != None:
        if x1 < x2: axis.set_xlim([x1,x2])
        else: axis.set_xlim([x2,x1])
    if y1 != None and y2 != None:
        if y1 < y2: axis.set_ylim([y1, y2])
        else: axis.set_ylim([y2, y1])


#--- SET SCALING RANGES AND DATA

def get_scaling_array(settings = None, x0 = 0.1, xend = 1.0, dx = 0.1):
    if settings == None:
        settings = user_settings
    vals = []
    length = int((xend-x0) / dx) + 1
    if settings['scaling_idx'] == 0:
        if cf.hamiltonian == 1: vals = range(12, 19, 2)
        else: vals = range(10, 17, 1)
    elif settings['scaling_idx'] == 5:
        vals = range(1, int(cf.params_arr[0] / 2) + 1)
    else :
        for x in range(0, length) :
            vals.append(x0 + x * dx)
    return np.array(vals)


def regspace(x0, xend, dx):
    return np.array(range(int(100 * x0), int(100 * xend), int(100 * dx))) / 100.

#--- SPECTRAL FORM FACTOR GOE SHAPE
def sff_GOE(x : np.array):
    return np.array([2*a-a*np.log(1+2*a) if a < 1 else 2-a * np.log((2*a+1)/(2*a-1)) for a in x])


#--- CREATE INSET FIGURE ON ANY SUBPLOT
def add_subplot_axes(ax,rect):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    #subax = fig.add_axes([x,y,width,height],facecolor=facecolor)  # matplotlib 2.0+
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.15
    y_labelsize *= rect[3]**0.15
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

#--- READ WEIRD FILE I SAVED VIA PYTHON

def read_python_saved_dat_file(filename):
    """
    Loading statistical data from scaling file (as funciton of paramter) and return 2D-array of them
    """
    status = True
    try:
        stats = pd.read_table(filename, sep="\t", header=None)
    except Exception as e:
        #print("Pandas broke down")
        status = False
    if status:
        result = []
        for i in range(len(stats) - 1):
            arr = np.array(stats[i][1:]).astype(float)
            result.append(arr)
    else:
        f = open(filename, "r")
        data = f.read()
        indices = findOccurrences(data, "\n")
        num_of_cols = len(findOccurrences(data[:indices[0]], "\t")) + 2  
        data = data[indices[0]:]
        data = data.split('\n')[1:]
        result = [[] for i in range(num_of_cols)]

        for line in data[:-1]:
            line = [int(x) if x.isdigit() else float(x) for x in line.split('\t')[:-1]]
            for i in range(0, num_of_cols):
                result[i].append(line[i])

    return np.array(result)


#--------------- LOAD STATISTICAL DATA
def load_stats(filename):
    """
    Loading statistical data to dictionairy from raw file
    """
    f = open(filename, "r")
    data = f.read()
    
    result = {}
    endline = [0] + findOccurrences(data, "\n")
    for i in range(len(endline) - 1):
        line = data[endline[i] : endline[i+1]]
        index = findOccurrences(line, "\t")
        x = line[index[0] + 1 :]
        result[line[2 if endline[i] > 0 else 1 : index[0]-1]] = int(x) if x.isdigit() else float(x)
    
    return result


#--------------- REMOVE FLUCTUATIONS FROM DATA
def remove_fluctuations(data, bucket_size=10):
    new_data = data;
    half_bucket = int(bucket_size / 2)
    for k in range(half_bucket, len(data) - half_bucket):
        average = np.sum(data[k - half_bucket : k + half_bucket])
        new_data[k - half_bucket] = average / bucket_size
    return new_data