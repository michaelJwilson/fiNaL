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
# results.calc_power_spectra(pars)

# Growth rate Normalised to (1. + z) in matter domination era.
Ds             = results.get_sigma8() / results.get_sigma8()[0] / (1. + zs[0])
Dz             = interp1d(zs, Ds, kind='linear', copy=True, bounds_error=True)

# Save growth rate to file. 
np.savetxt('dat/growth.txt', np.c_[zs, Ds], header='Linear growth rate;  Normalised to (1+z) in matter domination.')

