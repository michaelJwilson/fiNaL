import pyccl
import numpy              as     np
import pylab              as     pl
import matplotlib.pyplot  as     plt

from   pyccl.halos.hmfunc import MassFunc, MassFuncTinker08
from   pyccl.halos.hbias  import HaloBiasTinker10


##  https://ccl.readthedocs.io/_/downloads/en/latest/pdf/
##  https://github.com/LSSTDESC/CCL/blob/master/pyccl/halos/hmfunc.py

cosmo = pyccl.Cosmology(Omega_c=0.25, Omega_b=0.05,
                        h=0.7, n_s=0.95, sigma8=0.8,
                        transfer_function='bbks')

zz    = 2.0

pyccl.sigma8(cosmo)
pyccl.comoving_radial_distance(cosmo, 1./(1+zz))
pyccl.nonlin_matter_power(cosmo, k=1, a=0.5)

##  https://github.com/LSSTDESC/CCL/blob/master/pyccl/halos/hmfunc.py
hmf   = MassFuncTinker08(cosmo, mass_def=None, mass_def_strict=True)
hbf   = HaloBiasTinker10(cosmo, mass_def=None, mass_def_strict=True)

##  halo masses in units of M_sun.
Ms    = 10. ** np.arange(10., 14., 0.05)

##  Units of Mpc^-3.
Ns    = hmf.get_mass_function(cosmo, Ms, 1./(1+zz), mdef_other=None)
bs    = hbf.get_halo_bias(cosmo, Ms, 1./(1+zz), mdef_other=None)

fig, ax1 = plt.subplots()

ax1.loglog(Ms, Ns, label=r'$z=2$', c='k')
ax1.set_xlabel(r'$M \ \ [M_{\odot}]$')
ax1.set_ylabel(r'$\bar n \ \ [Mpc^{-3}]$')
ax1.legend(loc=2, frameon=False)

ax2   = ax1.twinx()
ax2.semilogx(Ms, bs, c='g')
ax2.set_ylabel('b(M)')

plt.tight_layout()
pl.savefig('halos.pdf')

