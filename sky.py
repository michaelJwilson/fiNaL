import matplotlib; matplotlib.use('PDF')

import numpy as np
import pylab as pl

from   astropy.cosmology    import  Planck15, FlatLambdaCDM
from   scipy                import  signal
from   astropy.convolution  import  convolve, Box1DKernel


cosmo      = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.05)

owave, sky =  np.loadtxt('spec-sky.dat', unpack=True)

# sky has 0.1A sampling; downsample to 4A.
ssky       =  convolve(sky, Box1DKernel(40))

# Rest-frame line
lya        =  1216.

redshift   = (owave / lya) - 1.0
chis       = Planck15.comoving_distance(redshift)
chis      /= 1000.

cut        = (redshift > 2.0) & (redshift < 5.0)

pl.semilogy(chis[cut], ssky[cut], lw=0.1, c='k')
pl.xlabel(r'$\chi \ [10^{3} \rm{Mpc}/h]$')
pl.ylabel('Sky background')
pl.ylim(0.5, 500.)
pl.savefig('back.pdf')
