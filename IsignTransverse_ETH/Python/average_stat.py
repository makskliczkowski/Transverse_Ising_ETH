import numpy as np
import pandas as pd
import os
from fnmatch import fnmatch as fn
from fnmatch import translate
import utils.helper_functions as hfun
from config import *

user_settings = getattr(plot_settings, 'settings')

if __name__ == '__main__':

    dir_out = f"{base_directory}STATISTICS{os.sep}"
    dir_in = f"{base_directory}STATISTICS{os.sep}raw_data{os.sep}"
    
    sizes = []
    collected_pars = []
    realis = []

    ii = user_settings['scaling_idx']
    info = ""

    # ASSUMING DISORDER MODULE FOR NOW
    if ii == 0:  info = f"_L=*,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f_jobid=*.dat"%(J, J0, g, g0, h, w)
    elif ii == 1:  info = f"_L=*,J=*,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f_jobid=*.dat"%(J0, g, g0, h, w)
    elif ii == 2:  info = f"_L=*,J=%.2f,J0=%.2f,g=*,g0=%.2f,h=%.2f,w=%.2f_jobid=*.dat"%(J, J0, g0, h, w)
    elif ii == 3:  info = f"_L=*,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=*,w=%.2f_jobid=*.dat"%(J, J0, g, g0, w)
    elif ii == 4:  info = f"_L=*,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=*_jobid=*.dat"%(J, J0, g, g0, h)
    else: info = f"_L=*,J=%.2f,J0=%.2fg=*,g0=%.2f,h=%.2f,w=%.2f_jobid=*.dat"%(J, J0, g0, h, w)

    #--------------------------------------- COLLECT DATA
    for filename in os.listdir(dir_in):
        if fn(filename, info):
            bare_info, extension = os.path.splitext(filename)
            f = os.path.join(dir_in, filename)
            if os.path.isfile(f):
                pars = hfun.get_params_from_info(bare_info)
                
                if pars[ii] not in collected_pars:  collected_pars.append(pars[ii])
                if pars[0] not in sizes:            sizes.append(int(pars[0]))
                if pars[-1] not in realis:          realis.append(int(pars[-1]))

    sizes = np.sort(sizes)
    collected_pars = np.sort(collected_pars)
    realis = np.sort(realis)
    print(sizes)
    print(collected_pars)
    print(realis)
    for L in sizes:
        new_pars = pars
        new_pars[0] = L
        info_out = hfun.remove_info(hfun.info_param(new_pars), user_settings['scaling']) + ".dat"

        file_scaling = open(dir_out + info_out, "w")
        file_scaling.write("'gap ratio'\t'ipr'\t'information entropy'\t'entropy in ~100 states at E=0\
                        '\t'mean level spacing'\t'typical level spacing'\t'entropy var in ~100 states at E=0'\t'entropy error over realisations'\n")
        
        # loop scaling parameter
        for par in collected_pars:
            file_scaling.write("%.8f\t"%par)
            new_pars[ii] = par
            info_in = hfun.info_param(new_pars)
            base_info, ext = os.path.splitext(info_in)
            
            info_in = base_info + "_jobid=%d"%realis[0] + ext
            stats = {
                'gap ratio':                            0.0,
                'ipr':                                  0.0,
                'information entropy':                  0.0,
                'entropy in ~100 states at E=0':        0.0,
                'mean level spacing':                   0.0,
                'typical level spacing':                0.0,
                'entropy var in ~100 states at E=0':    0.0,
                'entropy error over realisations':      0.0
            }
            counter = 0;
            for r in realis:
                info_in = base_info + "_jobid=%d"%r + ext
                if os.path.exists(dir_in + info_in):
                    counter += 1
                    data = pd.read_table(dir_in + info_in, sep="\t", header=None)
                    for i in range(len(data[0])):
                        indices = hfun.findOccurrences(data[0][i], "'")
                        variable_name = data[0][i][indices[0]+1:indices[1]]
                        stats[variable_name] += float(data[1][i])
            if counter == 0: 
                print(info_in)
                file_scaling.write("\n")
                continue

            for keys in stats:
                stats[keys] /= float(counter)

            #--------------------------------------- SAVE DATA AT GIVEN PARAMTER TO RAW DATA
            file = open(dir_in + info_out, "w")
            for keys in stats:
                file.write("'%s'"%keys)
                file.write("\t%.8f\n"%stats[keys])
            file.close()


            #--------------------------------------- SAVE DATA TO SCALING FILE
            for keys in stats:
                file_scaling.write("%.12f\t"%stats[keys])
            file_scaling.write("\n")
        file_scaling.close()

    a = 0.001131423423
    b = 0.32132142342
    c = 21.21235923042804298
    n = hfun.order_of_magnitude(a); print(str("{:.%df}"%n).format(round(a, n)))
    n = hfun.order_of_magnitude(b); print(str("{:.%df}"%n).format(round(b, n)))
    n = hfun.order_of_magnitude(c); print(str("{:.%df}"%n).format(round(c, n)))