import re
import yaml
import pylab             as     pl
import numpy             as     np
import itertools

from   scipy.interpolate import interp1d
from   get_schechters    import get_schechter


def get_halos(redshift, reverse=False):
    nbars    = np.loadtxt('dat/halo_lim_nbar.txt')

    infile   = open('dat/halo_lim_nbar.txt', 'r')
    hdr      = infile.readline()

    zs       = re.findall("\d+\.\d+", hdr)
    zs       = np.array(zs).astype(np.float)

    biases   = np.loadtxt('dat/halo_eff_biases.txt')

    cut      = zs == redshift
    
    idx      = np.where(cut)[0]
    idx     += 1

    #
    Ms       = nbars[:,0]
    ns       = nbars[:,idx][:,0]
    bs       = biases[:,idx][:,0]
    
    if reverse:
        xs   = interp1d(ns.T, Ms.T, bounds_error=True)
        ys   = interp1d(ns.T, ns.T, bounds_error=True)
        
    else:
        xs   = interp1d(Ms, ns, bounds_error=True)
        ys   = interp1d(Ms, bs, bounds_error=True)

    return  xs, ys


if __name__ == '__main__':
    ns       = np.logspace(-3., 0., 100)

    xs, ys   = get_halos(2.0, reverse=True)

    pl.loglog(xs(ns), ys(ns))

    pl.show()
    
    print('\n\nDone.\n\n')
