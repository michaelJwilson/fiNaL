import  camb
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    camb               import  model, initialpower
from    scipy.interpolate  import  interp1d
from    params             import  *


# https://camb.readthedocs.io/en/latest/CAMBdemo.html
zs             = np.arange(zmax, 0.0, -dz)

# Call CAMB.
pars           = camb.CAMBparams()

pars.InitPower.set_params(As=As, ns=ns, r=r)
pars.set_matter_power(redshifts=zs, kmax=10.0)
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=omk, tau=tau)

pars.NonLinear = model.NonLinear_none

results        = camb.get_results(pars)
results.calc_power_spectra(pars)

# Linear dark matter power spectrum. 
kh, z, plin    = results.get_matter_power_spectrum(minkh=1e-4, maxkh=5., npoints = 500, var1='delta_cdm', var2='delta_cdm', params=pars)

# Sub-sampling in redshift. 
rate           = 5

Ps             = [kh]

for zindex, z in enumerate(zs): 
  if (zindex % rate == 0) & (z > 0.0):
    # Eqn. (5) of https://arxiv.org/pdf/0807.1770.pdf; requires linear DM power spectrum. 
    Pk           = plin[-zindex,:]
    Ps.append(Pk.tolist())

#    
Ps    = np.array(Ps).T

hdr   = 'k'.ljust(23) + ''.join(['Lin. Pk(z={:.2f})'.format(zz).ljust(25) for zz in zs[::rate]])
np.savetxt('dat/Pks.txt',         Ps,   header=hdr)

