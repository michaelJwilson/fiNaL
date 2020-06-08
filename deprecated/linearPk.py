import re
import numpy  as np
import pylab  as pl

from   params import *


dat     = np.loadtxt('dat/Pks.txt')
ks      = dat[:,0]
Ps      = dat[:,1:]

infile  = open('dat/Pks.txt', 'r')
hdr     = infile.readline()

zs      = re.findall("\d+\.\d+", hdr)
zs      = np.array(zs).astype(np.float)

def Plin(z):
    idx = (np.abs(zs - z)).argmin()
    
    return  ks, Ps[:, idx]


if __name__ == '__main__':    
    for zz in [4.0, 3.0, 2.0, 1.0]:
      k, Pk = Plin(zz)
      pl.loglog(k, Pk, label=r'$z={:.1f}$'.format(zz))

      path = 'dat/Transfers_z{}.txt'.format(np.int(100. * zz))
                                            
      dat  = np.loadtxt(path)
      pl.loglog(dat[:,0], dat[:,-1] * h**3, 'k--')

    pl.xlim(1.e-3, 1.)
    pl.ylim(10., 5.e4)

    pl.legend(frameon=False)
    
    pl.show()
