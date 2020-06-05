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
    zz    = 4.00

    k, Pk = Plin(zz)

    pl.loglog(k, Pk)

    pl.show()
