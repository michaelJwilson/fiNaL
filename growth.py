from    scipy.interpolate  import  interp1d


# Generated by transfer.py
zs, Ds = np.loadtxt('dat/growth.txt')

Dz     = scipy.interp1d(zs, Ds, copy=True, bounds_error=True)
    
