import camb
import numpy             as     np
import pylab             as     pl
import matplotlib.pyplot as     plt

from   camb              import model, initialpower
from   scipy.interpolate import interp1d


# https://camb.readthedocs.io/en/latest/CAMBdemo.html
H0             = 67.6          # [km/s/Mpc].
h              = H0 / 100.

# eBOSS QSO. 
# b, zz        = 2.58, 1.54
# nbar         = 1.e-4

# Linear today. 
# b, zz        = 2.58, 1.54 
# nbar         = 1.e-2

# Goldrush bright.
mlim           = 24.0
b, zz          = 6.06, 4.00
nbar           = 1.e-4

# Goldrush faint.
# mlim         = 25.5
# b, zz        = 4.07, 4.00
# nbar         = 1.e-2

p              = 1.0           # [1., 1.6]
dc             = 1.61
fnl            = 20.0

c              = 2.9979 * 1.e5 # [km/s].
ombh2          = 0.022
omch2          = 0.122
mnu            = 0.000

Om             = (omch2 + ombh2) / h /h

zs             = np.arange(7.0, 0.0, -0.05)
zindex         = (np.abs(zs - zz)).argmin()

pars           = camb.CAMBparams()

pars.InitPower.set_params(As=2e-9, ns=0.965, r=0.0)
pars.set_matter_power(redshifts=zs, kmax=2.0)
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=0.0, tau=0.06)

pars.NonLinear = model.NonLinear_none

results        = camb.get_results(pars)
results.calc_power_spectra(pars)

kh, z, plin    = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200, var1='delta_cdm', var2='delta_cdm', params=pars)
scalar         = pars.scalar_power(kh)

trans_cdm      = plin[-zindex,:] / scalar / kh
trans_cdm     /= trans_cdm[0]
trans_cdm      = np.sqrt(trans_cdm)

Ds             = results.get_sigma8() / results.get_sigma8()[0] / (1. + zs[0])
Dz             = interp1d(zs, Ds, kind='linear', copy=True, bounds_error=True)

# pl.loglog(zs, Ds, lw=0.1)

# Matter domination. 
# pl.loglog(zs, 1. / (1. + zs), lw=0.1)
# pl.show()
'''
transfer     = camb.get_transfer_functions(pars).get_matter_transfer_data()
# back         = camb.get_background(pars)

kh             = transfer.transfer_z('k/h')               ##  [(Mpc/h)^-1].

trans_cdm      = transfer.transfer_z('delta_cdm')
trans_cdm     /= trans_cdm[0]
'''
k2T            = kh * kh * trans_cdm

dbk            = 3. * (b - p) * dc * Om * H0**2 / c /c / k2T / Dz(3.0)
Pk             = plin[-zindex,:] * (b + dbk * fnl)**2.

#k             = 10.**np.linspace(-5, 1, 50)
#scalar        = pars.scalar_power(k)
                  
pl.axhline(1. / nbar, xmin=0.0, xmax=1.0, c='k', lw=0.1)

pl.loglog(kh, Pk, label=r'$(z, m, b, \bar n) = ({}, {}, {}, {:.1e})$'.format(zz, mlim, b, nbar))
pl.loglog(kh, b * b * plin[-zindex,:]) 
pl.loglog(kh,         plin[-zindex,:])

pl.xlabel(r'$k \ [h/{\rm Mpc}]$')
pl.ylabel(r'$P(k)$')

pl.legend(frameon=False)

pl.xlim(1.e-3, 0.2)
pl.ylim(100., 1.e6)

plt.tight_layout()

pl.savefig('pk_fnl.pdf')
