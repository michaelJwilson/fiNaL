import camb
import pylab    as      pl

from   params   import  *
from   camb     import  model, initialpower
from   growth   import  Dz
from   linearPk import  Plin


# Call CAMB.                                                                                                                                                                                                                              
pars           = camb.CAMBparams()

pars.InitPower.set_params(As=As, ns=ns, r=r)

# Get transfer function at z = 0
pars.set_matter_power(redshifts=[0.0], kmax=10.0)
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=omk, tau=tau)
  
pars.NonLinear = model.NonLinear_none

results        = camb.get_results(pars)
results.calc_power_spectra(pars)

trans           = results.get_matter_transfer_data()

kh              = trans.transfer_data[0,:,0]  ##  [(Mpc/h)^-1].                                                                                                                                             
k               = kh * results.Params.h
k2              = k * k

# model.Transfer_cdm = 2                                                                                                                                                                                                                  
transfer_cdm    = trans.transfer_data[2-1,:,0]
transfer_tot    = trans.transfer_data[model.Transfer_tot-1,:,0]

def alpha(z):
    return  5./3. * (H0 / 100.)**2 * kh**2 * transfer_tot * Dz(z)

def beta(b):
    return  2 * dc * (b - 1)

def dbk(b, z):
    return  beta(b) / alpha(z)

def phh(b, z, fnl):
    ks,  ps = Plin(z)

    assert  np.all(kh == ks)
    
    return  (b + fnl * dbk(b, z)) * ps


if __name__ == '__main__':
    b       = 1.2
    z       = 1.0

    ps      = phh(b, z, 1.0)
    
    pl.loglog(kh, ps)

    pl.show()
