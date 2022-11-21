#!/usr/bin/env python

"""



"""
import numpy as np
import sys

if __name__ == '__main__':
    seed = int(sys.argv[1])
    jobid = int(sys.argv[2])
    Ndis = 10
    for r in range(Ndis):
        "diagonalize"
        "observable"
        "save data with suffix jobid=(r + 10 * jobid)"
        #eg. filename = "blabla..._jobid=%d.dat"%(r + 10*jobid) or .npy or .hdf5