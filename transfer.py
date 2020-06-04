import  camb
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    camb               import  model, initialpower
from    scipy.interpolate  import  interp1d


# https://camb.readthedocs.io/en/latest/CAMBdemo.html
As             = 2.e-9
ns             = 0.965 
r              = 0.0

zmax           = 7.0
dz             = 0.05

H0             = 67.6           # [km/s/Mpc].
h              = H0 / 100.

c              = 2.9979 * 1.e5  # [km/s].
omk            = 0.0
ombh2          = 0.022
omch2          = 0.122
mnu            = 0.000
tau            = 0.06

Om             = (omch2 + ombh2) / h /h

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
scalar         = pars.scalar_power(kh)

# Normalised to (1. + z) in matter domination era.                                                                                                                                                                                        
Ds             = results.get_sigma8() / results.get_sigma8()[0] / (1. + zs[0])
Dz             = interp1d(zs, Ds, kind='linear', copy=True, bounds_error=True)

# Save growth rate to file. 
np.savetxt('growth.txt', np.c_[zs, Ds], header='Linear growth rate;  Normalised to (1+z) in matter domination.')

# [h/Mpc].
k2Ts           = [kh]
Ps             = [kh]

rate           = 20

# Transfer function calc. with k. 
for zindex, zz in enumerate(zs):
  if zindex == 0:    
    trans_cdm  = plin[-zindex,:] / scalar / kh
    trans_cdm /= trans_cdm[0]
    trans_cdm  = np.sqrt(trans_cdm)

    '''                                                                                                                                                                                                                     
    ##  Equates to transfer function above.
    transfer   = camb.get_transfer_functions(pars).get_matter_transfer_data()                                                                                                                                             
    # back     = camb.get_background(pars)                                                                                                                                                                               
                                                                                                                                                                                                                       
    kh         = transfer.transfer_z('k/h')               ##  [(Mpc/h)^-1].                                                                                                                                              
                                                                                                                                                                                                                    
    trans_cdm  = transfer.transfer_z('delta_cdm')                                                                                                                                                                       
    trans_cdm /= trans_cdm[0]                                                                                                                                                                                             
    '''

    k2T        = kh * kh * trans_cdm
    k2T        = k2T / h / h

    k2Ts.append(k2T.tolist())

  if zindex % rate == 0:
    # Eqn. (5) of https://arxiv.org/pdf/0807.1770.pdf; requires linear DM power spectrum. 
    Pk         = plin[-zindex,:]
    Ps.append(Pk.tolist())

k2Ts = np.array(k2Ts).T
Ps   = np.array(Ps).T

hdr  = 'k'.ljust(23) + ''.join(['k2T(z={:.2f})'.format(zz).ljust(25) for zz in zs[::rate]])
np.savetxt('k2Transfers.txt', k2Ts, header=hdr)

hdr  = 'k'.ljust(23) + ''.join(['Lin. Pk(z={:.2f})'.format(zz).ljust(25) for zz in zs[::rate]])
np.savetxt('Pks.txt',         Ps,   header=hdr)
    
