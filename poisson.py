import matplotlib.pyplot  as     plt
import numpy              as     np
import pylab              as     pl

from   nbodykit.lab       import *
from   nbodykit           import style, setup_logging

setup_logging()


redshift = 0.55
cosmo    = cosmology.Planck15
Plin     = cosmology.LinearPower(cosmo, redshift, transfer='EisensteinHu')
b1       = 1.0

# Produce a pure shot noise catalogue. 
def nopower(k):
    return  0.0 * np.ones_like(k)

nbar      = 1.e-3

cat       = LogNormalCatalog(Plin=Plin, nbar=nbar, BoxSize=1000., Nmesh=256, bias=b1, seed=42, cosmo=cosmo, redshift=redshift)

# Mesh normalised as (n / <n> ) or (1 + delta).
rmesh     = cat.to_mesh(compensated=True, window='cic', position='Position', interlaced=False)
rmesh     = rmesh.compute()

# Every second element in last axis to be zero;  Not a copy.
#rmesh.value[:,:,::2] = 0.0
#rmesh.value[:,::2,:] = 0.0
#rmesh.value[::2,:,:] = 0.0

# Actual postions for each cell.
# rmesh.x

r         = FFTPower(rmesh, mode='1d', dk=0.05)

pl.loglog(r.power['k'], r.power['power'])
pl.axhline(1. / nbar, c='k')
pl.axhline(1. / nbar / 8., c='c')

pl.xlim(0.01,  0.5)
pl.ylim(1., 5.e4)

pl.show()

