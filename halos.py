import pyccl
import numpy              as     np
import pylab              as     pl
import matplotlib.pyplot  as     plt

from   pyccl.halos.hmfunc import MassFunc, MassFuncTinker08
from   pyccl.halos.hbias  import HaloBiasTinker10
from   params             import *


##  https://ccl.readthedocs.io/_/downloads/en/latest/pdf/
##  https://github.com/LSSTDESC/CCL/blob/master/pyccl/halos/hmfunc.py

fig, ax1  = plt.subplots()
ax2       = ax1.twinx()

cosmo     = pyccl.Cosmology(Omega_c=omch2/h/h, Omega_b=ombh2/h/h,
                            h=h, n_s=ns, sigma8=0.8,
                            transfer_function='bbks')

colors    = plt.rcParams['axes.prop_cycle'].by_key()['color']

##  halo masses in units of M_sun.                                                                                                                      
Ms        = 10. ** np.arange(9., 13., 0.05)

nbars     = [Ms]
biases    = [Ms]

for zz, c in zip([1.0, 2.0, 3.0, 4.0], colors):
    ##  https://github.com/LSSTDESC/CCL/blob/master/pyccl/halos/hmfunc.py
    hmf   = MassFuncTinker08(cosmo, mass_def=None, mass_def_strict=True)
    hbf   = HaloBiasTinker10(cosmo, mass_def=None, mass_def_strict=True)

    ##  dn/dlog10M.  Units of Mpc^-3.
    Ns    = hmf.get_mass_function(cosmo, Ms, 1./(1+zz), mdef_other=None)
    bs    = hbf.get_halo_bias(cosmo, Ms, 1./(1+zz), mdef_other=None)

    ##  Units of (Mpc/h)^-3.  
    Ns   /= h**3.

    ##  dn/dlog10M -> dn/dlogM
    Ns   /= np.log(10.)
    
    nbars.append(Ns)
    biases.append(bs)
    
    ax1.loglog(Ms, Ns, label=r'$z={}$'.format(zz), c=c)
    ax2.semilogx(Ms, bs, c=c, alpha=0.5)

header  = 'M [Msun]'.ljust(23)
header += ''.join(['dn/dlnM({:.1f})[(Mpc/h)^-3]'.format(zz).ljust(25) for zz in [1.0, 2.0, 3.0, 4.0]])

nbars   = np.array(nbars).T
biases  = np.array(biases).T

np.savetxt('dat/halo_nbar.txt', nbars,  header=header)

header  = 'M [Msun]'.ljust(23)
header += ''.join(['b(z={:.1f})'.format(zz).ljust(25) for zz in [1.0, 2.0, 3.0, 4.0]])

np.savetxt('dat/halo_biases.txt', biases, header=header)

ax1.set_xlabel(r'$M_h \ \ [M_{\odot}]$')
ax1.set_ylabel(r'$d \bar n / d \ln(M) \ \ [({\rm Mpc}/h)^{-3}]$')

ax1.set_xlim(1.e10, 1.e13)
ax1.set_ylim(1.e-6, 1.0)

ax1.legend(loc=1, frameon=False)

ax2.set_ylabel(r'$b(M_h)$')
ax2.set_ylim(0.0, 15.0)

ax2.axhline(y=1.0, c='k', lw=0.5)

plt.tight_layout()

pl.savefig('plots/halos.pdf')

