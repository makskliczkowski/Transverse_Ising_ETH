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

    vals = hfun.get_scaling_array(settings=settings)
    entropy = []
    xarray = []
    param_copy = copy.deepcopy(cf.params_arr)
    for x in vals:
        cf.params_arr[settings['scaling_idx']] = x
        filename = cf.base_directory + "STATISTICS" + os.sep + hfun.remove_info(hfun.info_param(cf.params_arr), settings['vs']) + ".dat" 
        if os.path.exists(filename):
            stats = pd.read_table(filename, sep="\t", header=None)
            xarray.append(np.array(list(stats[0])[1:]).astype(float))

            S = np.array(list(stats[4])[1:]).astype(float)
            norm = x * np.log(2) / 2. + (0.5 - np.log(2)) / 2. - 0.5
            entropy.append(S / norm)

    #--- reset defaults
    cf.params_arr = param_copy

    vals, xvals, tau, gap_ratio = thouless.load(settings, vals)

    dir = "CriticalParameters"
    try:
        os.mkdir(dir)
    except OSError as error:
        print(error)   
    def calculate_and_save(scaling_ansatz, crit_fun = 'free'):
        suffix = hfun.remove_info(hfun.info_param(cf.params_arr), settings['scaling'], settings['vs']) \
                     + "_critfun=%s_ansatz=%s_pert=%s_seed=%d"%(crit_fun, scaling_ansatz, settings['vs'], seed)

        #-- calculate gap ratio collapse
        par, crit_pars, costfun, status = cost.get_crit_points(x=xvals, y=gap_ratio, vals=vals, crit_fun=crit_fun, scaling_ansatz=scaling_ansatz, seed=seed)
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
        par, crit_pars, costfun, status = cost.get_crit_points(x=np.array(xarray), y=np.array(entropy), vals=vals, crit_fun=crit_fun, scaling_ansatz=scaling_ansatz, seed=seed)
        print("Entropy:\t", par, crit_pars, costfun, status)
        if status:
            filename = dir + os.sep + "Entropy" + suffix
            data = {
                "costfun": costfun,
                "crit exp'": par
            }
            for i in range(len(crit_pars)):
                data["x_%d"%i] = crit_pars[i]
            np.savez(filename, **data)


    calculate_and_save('classic')
    calculate_and_save('RG')
    calculate_and_save('KT')
    calculate_and_save('FGR')
