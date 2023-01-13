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
    print("seed = ", seed)

    dir = "CriticalParameters_Bartek"
    try:
        os.mkdir(dir)
    except OSError as error:
        print(error)  
    def collapse_data_eh(scaling_ansatz, crit_fun = 'free', W=10.0):
        suffix = "_W=%.1f_critfun=%s_ansatz=%s_seed=%d"%(W, crit_fun, scaling_ansatz, seed)
        sizes = np.arange(10, 17, 2)
        gvals = []
        gap_ratio = []
        for L in sizes:
            name = f'BARTEK_DATA/L{L}_W{W}_sqrt.dat'
            if os.path.exists(name):
                data = pd.read_table(name, sep=" ", header=None)
                gvals.append(data[0])
                gap_ratio.append(data[1])

        #-- calculate gap ratio collapse
        par, crit_pars, costfun, status = cost.get_crit_points(x=gvals, y=gap_ratio, vals=sizes, crit_fun=crit_fun, scaling_ansatz=scaling_ansatz, seed=seed)
        print("Gap Ratio:\t", par, crit_pars, costfun, status)
        if costfun > 2: 
            status = Failed
            assert False, "Cost function too high: CF = %d"%costfun
        if status:
            filename = dir + os.sep + "GapRatio" + suffix
            data = {
                "costfun": costfun,
                "crit exp'": par
            }
            for i in range(len(crit_pars)):
                data["x_%d"%i] = crit_pars[i]
            np.savez(filename, **data)

    for W in [3.0, 5.0, 10.0, 20.0]:
        collapse_data_eh(scaling_ansatz='classic', crit_fun='free', W=W)
        collapse_data_eh(scaling_ansatz='FGR', crit_fun='free', W=W)
        collapse_data_eh(scaling_ansatz='KT', crit_fun='free', W=W)
        collapse_data_eh(scaling_ansatz='RG', crit_fun='free', W=W)
