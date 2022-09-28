import importlib
from multiprocessing.sharedctypes import Value
import config as cf
from matplotlib.markers import MarkerStyle
import matplotlib.pyplot as plt
import numpy as np
importlib.reload(cf)

model = 0   # chooses model: 0-disorder / 1-symmetries
BC = 1     # boundaary condition: 0 - OBC / 1 - PBC
user_settings = getattr(cf.plot_settings, 'settings')
var_name = "\\Delta" if cf.hamiltonian else "g"

#-------------------------- SET INFO
def info_sym(L, J, g, h, k, p, x):
    return "_L=%d,J=%.2f,g=%.2f,h=%.2f,k=%d,p=%d,x=%d.dat"%(L, J, g, h, k, p, x)
def info_dis(L, J, J0, g, g0, h, w):
    return "_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat"%(L, J, J0, g, g0, h, w)

def info(_L = cf.params_arr[0], _J = cf.params_arr[1], _J0 = cf.params_arr[8], 
            _g = cf.params_arr[2], _g0 = cf.params_arr[9], 
            _h = cf.params_arr[3], _w = cf.params_arr[4], 
          _k = cf.params_arr[5], _p = cf.params_arr[6], _x = cf.params_arr[7]
    ):
    if(model == 1) :
        return info_sym(_L, _J, _g, _h, _k, _p, _x)
    else :
        return info_dis(_L, _J, _J0, _g, _g0, _h, _w)
def info_param(params):
    return info(_L = params[0], _J = params[1], _J0 = params[8], 
            _g = params[2], _g0 = params[9], 
            _h = params[3], _w = params[4], 
          _k = params[5], _p = params[6], _x = params[7]
    )

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
    return r"$" + (scaling_str + (f"=%d"%(x) if settings['scaling_idx'] == 0 or settings['scaling_idx'] == 5 else f"=%.2f"%(x))) + "$"

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
def set_plot_elements(axis, xlim =[], ylim=[], xlabel = None, ylabel = 'y', settings = None, set_legend = True, font_size = 10):
    if settings == None:
        settings = user_settings
    
    xlab = list(settings['func_x_name'])
    x_idx = xlab.index('Q')
    if xlabel == None :
        xlab[x_idx] = settings['vs']
    else :
        xlabel = list(xlabel)
        xlab.remove('Q')
        for i in range(0,len(xlabel)):
            xlab.insert(i + x_idx, xlabel[i]) 
    ylab = list(settings['func_y_name'])
    ylab[ylab.index('Q')] = ylabel
    
    axis.set_ylabel("".join(ylab), rotation=90, fontsize=font_size+2)
    axis.set_xlabel("".join(xlab), rotation=0, fontsize=font_size+2, labelpad=font_size-8)
    axis.set_yscale(settings['y_scale'])
    axis.set_xscale(settings['x_scale'])
    axis.tick_params(axis='both', which='major', labelsize=font_size, length=font_size-4, width=0.05*font_size)
    axis.tick_params(axis='both', which='minor', labelsize=font_size, length=0.3*(font_size-4), width=0.05*font_size)
    
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
        if cf.hamiltonian: vals = range(12, 19, 2)
        else: vals = range(11, 17, 1)
    elif settings['scaling_idx'] == 5:
        vals = range(1, int(cf.params_arr[0] / 2) + 1)
    else :
        for x in range(0, length) :
            vals.append(x0 + x * dx)
    return np.array(vals)
