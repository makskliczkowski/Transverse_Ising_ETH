#!/usr/bin/env python

"""



"""
import numpy as np
import pandas as pd
import copy
import sys
import os
import utils.helper_functions as hfun
import config as cf
import modules.thouless_times as thouless
import costfun.costfun as cost

if __name__ == '__main__':
    seed = int(sys.argv[1])
    vs = sys.argv[2]
    print("seed = ", seed)
    print("vs = ", vs)

    set_class = copy.deepcopy(cf.plot_settings)
    set_class.set_scaling('L')
    set_class.set_vs(vs)
    settings = getattr(set_class, 'settings')

    old_vals = hfun.get_scaling_array(settings=settings)
    entropy = []
    xarray = []
    gap_ratio = []
    vals = []
    param_copy = copy.deepcopy(cf.params_arr)
    for x in old_vals:
        cf.params_arr[settings['scaling_idx']] = x
        filename = cf.base_directory + "STATISTICS" + os.sep + hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs']) + ".dat" 
        if os.path.exists(filename):
            #stats = pd.read_table(filename, sep="\t", header=None)
            stats = hfun.read_python_saved_dat_file(filename)
            #print(stats)
            r_tmp = stats[1]
            
            xvals = np.array([stats[0][i] for i, r in enumerate(r_tmp) if r > 0.42])
            xarray.append(xvals[:25])
            r = np.array([stats[1][i] for i, r in enumerate(r_tmp) if r > 0.42])
            gap_ratio.append(r[:25])
            
            S = np.array([stats[4][i] for i, r in enumerate(r_tmp) if r > 0.42])
            norm = x * np.log(2) / 2. + (0.5 - np.log(2)) / 2. - 0.5
            entropy.append( S[:25] / norm )
            vals.append(x)
    vals = np.array(vals)
    
    #--- reset defaults
    cf.params_arr = param_copy

    dir = "CriticalParameters"
    try:
        os.mkdir(dir)
    except OSError as error:
        print(error)   
    def calculate_and_save(scaling_ansatz, crit_fun = 'free'):
        suffix = hfun.remove_info(hfun.info_param(cf.params_arr), settings['scaling'], settings['vs']) \
                     + "_critfun=%s_ansatz=%s_pert=%s_seed=%d"%(crit_fun, scaling_ansatz, settings['vs'], seed)

        #-- calculate gap ratio collapse
        par, crit_pars, costfun, status = cost.get_crit_points(x=xarray, y=gap_ratio, vals=vals, crit_fun=crit_fun, scaling_ansatz=scaling_ansatz, seed=seed)
        print("Gap Ratio:\t", par, crit_pars, costfun, status)
        if status:
            filename = dir + os.sep + "GapRatio" + suffix
            data = {
                "costfun": costfun,
                "crit exp'": par
            }
            for i in range(len(crit_pars)):
                data["x_%d"%i] = crit_pars[i]
            np.savez(filename, **data)

        #-- calculate entropy collapse
        #par, crit_pars, costfun, status = cost.get_crit_points(x=np.array(xarray), y=np.array(entropy), vals=vals, crit_fun=crit_fun, scaling_ansatz=scaling_ansatz, seed=seed)
        #print("Entropy:\t", par, crit_pars, costfun, status)
        #if status:
        #    filename = dir + os.sep + "Entropy" + suffix
        #    data = {
        #        "costfun": costfun,
        #        "crit exp'": par
        #    }
        #    for i in range(len(crit_pars)):
        #        data["x_%d"%i] = crit_pars[i]
        #    np.savez(filename, **data)


    calculate_and_save(scaling_ansatz='FGR', crit_fun='free')
    #calculate_and_save(scaling_ansatz='KT', crit_fun='free')
    #calculate_and_save(scaling_ansatz='RG', crit_fun='free')
    #calculate_and_save(scaling_ansatz='classic', crit_fun='free')
#
    #calculate_and_save(scaling_ansatz='FGR', crit_fun='free_inv')
    #calculate_and_save(scaling_ansatz='KT', crit_fun='free_inv')
    #calculate_and_save(scaling_ansatz='RG', crit_fun='free_inv')
    #calculate_and_save(scaling_ansatz='classic', crit_fun='free_inv')