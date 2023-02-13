import numpy as np
import pandas as pd
import os, sys, copy, re
from fnmatch import fnmatch as fn
from fnmatch import translate
import utils.helper_functions as hfun
from config import *
import h5py    

if __name__ == '__main__':
    vs = sys.argv[1]
    print("vs = ", vs)

    set_class = copy.deepcopy(plot_settings)
    set_class.set_scaling(vs)
    
    settings = getattr(set_class, 'settings')
    print(base_directory)
    sys.stdout.flush()
    operator = settings['operator']
    site = settings['site']
    sub_dir = ""
    if site < 0 and operator < 8:
        sub_dir = operator_names[operator] + "%d"%(L/2) + kPSep
    else :
        sub_dir = operator_name(operator, site) + kPSep

    dir_out = f"{base_directory}AGP{os.sep}{sub_dir}"
    dir_in = f"{base_directory}AGP{os.sep}{sub_dir}raw_data{os.sep}"
    print(dir_in, dir_out)
    sizes = []
    collected_pars = []
    realis = []
    sys.stdout.flush()
    ii = settings['scaling_idx']
    info = ""

    # ASSUMING DISORDER MODULE FOR NOW
    arr = [J, J0, g, g0, h, w]
    names = ['J', 'J0', 'g', 'g0', 'h', 'w']
    indices = [-1, 0, 2, 4, 5]
    def create_info(use_log_data = True):
        info = "_L=*"
        for i, var in enumerate(arr):
            n = hfun.order_of_magnitude(var) if use_log_data else 2
            if i == indices[ii]:
                info += ",%s=*"%names[i]
            else: info += str(",%s={:.%df}"%(names[i], n)).format(round(var, n))
        return info + "_jobid=*.dat"
    info = create_info()
    print(info)
    #--------------------------------------- COLLECT DATA
    def append_parameter_range(filename):
        if fn(filename, info):
            print(filename)
            bare_info, _ = os.path.splitext(filename)
            
            f = os.path.join(dir_in, filename)
            if os.path.isfile(f):
                pars = hfun.get_params_from_info(bare_info)
                
                if pars[ii] not in collected_pars:  collected_pars.append(pars[ii])
                if pars[0] not in sizes:            sizes.append(int(pars[0]))
                if pars[-1] not in realis:          realis.append(int(pars[-1]))
    for filename in os.listdir(dir_in):
        info = create_info()
        append_parameter_range(filename)
        info = create_info(False)
        append_parameter_range(filename)
    pars = params_arr
    sizes = np.sort(sizes)
    collected_pars = np.sort(collected_pars)
    realis = np.sort(realis)
    print(sizes)
    print(collected_pars)
    print(realis)
    sys.stdout.flush()
    
    sizes = sizes if ii != 0 else [0];
    for L in sizes:
        new_pars = pars
        if ii != 0:
            new_pars[0] = L
        info_out = hfun.remove_info(hfun.info_param(new_pars), settings['scaling']) + ".dat"

        file_scaling = open(dir_out + info_out, "w")
        file_scaling.write("'adiabatic gauge potential'\t'typical fidelity susceptibility'\t'fidelity susceptibility'\t'diagonal normalisation factor'\n")
        
        # loop scaling parameter
        for par in collected_pars:
            file_scaling.write("%.8f\t"%par)
            new_pars[ii] = par
            info_in = hfun.info_param(new_pars)
            base_info, ext = os.path.splitext(info_in)

            info_in2 = hfun.info_param(new_pars, use_log_data=False)
            base_info2, ext2 = os.path.splitext(info_in2)
            
            stats = {
                'adiabatic gauge potential':        0.0,
                'typical fidelity susceptibility':  0.0,
                'fidelity susceptibility':          0.0,
                'diagonal normalisation factor':    0.0
            }
            counter = {
                'adiabatic gauge potential':        0,
                'typical fidelity susceptibility':  0,
                'fidelity susceptibility':          0,
                'diagonal normalisation factor':    0
            }
            for r in realis:
                base = base_info + "_jobid=%d"%r + ext
                info_in = base if os.path.exists(dir_in + base) else base_info2 + "_jobid=%d"%r + ext2
                if os.path.exists(dir_in + info_in):
                    data = pd.read_table(dir_in + info_in, sep="\t", header=None)
                    for i in range(len(data[0])):
                        indices = hfun.findOccurrences(data[0][i], "'")
                        variable_name = "".join(list(data[0][i][indices[0]+1:indices[1]]))
                    
                        if np.isnan(float(data[1][i])) == False and np.abs(float(data[1][i])) > 0.0:
                            stats[variable_name] += float(data[1][i])
                            counter[variable_name] += 1
                 
            #--------------------------------------- SAVE DATA TO SCALING FILE
            for keys in stats:
                if counter[keys] == 0:
                    file_scaling.write("nan\t")
                    stats[keys] = 1e20
                else:
                    stats[keys] /= float(counter[keys])
                    file_scaling.write("%.12f\t"%stats[keys])
            file_scaling.write("\n")
            
            #--------------------------------------- SAVE DATA AT GIVEN PARAMTER TO RAW DATA
            file = open(dir_in + hfun.info_param(new_pars), "w")
            for keys in stats:
                file.write("'%s'"%keys)
                if stats[keys] > 1e15:  file.write("\tnan\n")
                else:                   file.write("\t%.8f\n"%stats[keys])
            file.close()


        file_scaling.close()