#!/usr/bin/env python

"""



"""
import utils.helper_functions as hfun
import config as cf
import thouless_times as thouless
import costfun.costfun as cost
import numpy as np
import copy
import sys
import os


if __name__ == '__main__':
    seed = int(sys.argv[1])

    set_class = copy.deepcopy(cf.plot_settings)
    set_class.set_scaling('L')
    set_class.set_vs('g')
    settings = getattr(set_class, 'settings')

    #rescale_fun = cost.resc_functions_dict[scaling_ansatz]
    vals = hfun.get_scaling_array(settings=settings)
    vals, xvals, tau, gap_ratio = thouless.load(settings, vals)

    dir = "CriticalParameters"
    try:
        os.mkdir(dir)
    except OSError as error:
        print(error)   
    def calculate_and_save(scaling_ansatz):
        par, crit_pars, costfun = cost.get_crit_points(x=xvals, y=gap_ratio, vals=vals, scaling_ansatz=scaling_ansatz, seed=seed)
        filename = dir + os.sep + hfun.remove_info(hfun.info_param(cf.params_arr), settings['scaling'], settings['vs']) + ".dat" + "_critfun=free_ansatz=%s_pert=%s_seed=%d.dat"%(scaling_ansatz, settings['vs'], seed)
        data = {
            "costfun": costfun,
            "crit exp'": par
        }
        for i in range(len(crit_pars)):
            data["x_%d"%i] = crit_pars[i]
        np.savez(filename, **data)
    calculate_and_save('classic')