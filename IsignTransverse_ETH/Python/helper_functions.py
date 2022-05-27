import os
import matplotlib
from matplotlib.markers import MarkerStyle
import matplotlib.pyplot as plt
model = 0   # chooses model: 0-disorder / 1-symmetries
BC = 1     # boundaary condition: 0 - OBC / 1 - PBC

#-------------------------- SET INFO
def info_sym(L, J, g, h, k, p, x):
    return "_L=%d,J=%.2f,g=%.2f,h=%.2f,k=%d,p=%d,x=%d.dat"%(L, J, g, h, k, p, x)
def info_dis(L, J, J0, g, g0, h, w):
    return "_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat"%(L, J, J0, g, g0, h, w)

def info(_L, _J, _J0, _g, _g0, _h, _w, _k, _p, _x) :
    if(model) :
        return info_sym(_L, _J, _g, _h, _k, _p, _x)
    else :
        return info_dis(_L, _J, _J0, _g, _g0, _h, _w)


#-------------------------- DIR
kPSep = os.sep
dir = f"..{kPSep}results{kPSep}" + (f"symmetries{kPSep}" if model else f"disorder{kPSep}") + (f"PBC{kPSep}" if BC else f"OBC{kPSep}") 

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