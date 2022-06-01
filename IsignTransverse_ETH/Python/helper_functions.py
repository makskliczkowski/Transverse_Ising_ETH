import importlib
from multiprocessing.sharedctypes import Value
import config as cf
from matplotlib.markers import MarkerStyle
import matplotlib.pyplot as plt
importlib.reload(cf)

model = 0   # chooses model: 0-disorder / 1-symmetries
BC = 1     # boundaary condition: 0 - OBC / 1 - PBC
user_settings = getattr(cf.plot_settings, 'settings')

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
    if(model) :
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
def set_plot_elements(axis, xlim =[], ylim=[], xlabel = 'x', ylabel = 'y', settings = None, set_legend = True, font_size = 10):
    if settings == None:
        settings = user_settings
    xlab = ( r"$1\ /\ %s^{%.d}$"%(settings['vs'], settings['nu']) if settings['rescaleX'] else settings['vs'] ) if xlabel == 'x' else xlabel
    axis.set(xlabel = xlab, ylabel = ylabel)
    axis.set_yscale(settings['y_scale'])
    axis.set_xscale(settings['x_scale'])
    
    if set_legend:
        axis.legend(frameon=False
                , loc='best'
                , fontsize=font_size)
                
    x1, x2 = xlim
    y1, y2 = ylim
    if x1 != None and x2 != None:
        if x1 < x2: axis.set_xlim([x1,x2])
        else: axis.set_xlim([x2,x1])
    if y1 != None and y2 != None:
        if y1 < y2: axis.set_ylim([y1, y2])
        else: axis.set_ylim([y2, y1])