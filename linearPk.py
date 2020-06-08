import re
import numpy  as np
import pylab  as pl

from   params import *


def Plin(z, printit=False):
    if z < 1.0:
        raise ValueError('z < 1.0 is unsupported.')
    
    close = np.around(z, decimals=1)
    dat   = np.loadtxt('dat/Transfers_z{}.txt'.format(np.int(100. * close)))

    if printit:
        print(z, close)
    
    ks    = dat[:, 0]  ##  [h/Mpc]
    Ps    = dat[:,-1]  ##  [(Mpc/h)^3]
    
    return  ks, Ps


if __name__ == '__main__':    
    for zz in [4.0, 3.0, 2.0, 1.0]:
      k, Pk = Plin(zz, printit=True)
      pl.loglog(k, Pk, label=r'$z={:.1f}$'.format(zz))

    pl.xlim(1.e-3, 1.)
    pl.ylim(10., 5.e4)

    pl.legend(frameon=False)
    
    pl.show()
