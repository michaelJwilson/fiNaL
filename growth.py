import  numpy as np

from    scipy.interpolate  import  interp1d


# Generated by transfer.py
zs, Ds = np.loadtxt('dat/growth.txt', unpack=True)

Dz     = interp1d(zs, Ds, copy=True, bounds_error=True)
    