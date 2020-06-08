import  camb
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    camb               import  model, initialpower
from    scipy.interpolate  import  interp1d
from    params             import  *


# Cell 48 of https://camb.readthedocs.io/en/latest/CAMBdemo.html
for zz in np.arange(5.0, 0.9, -0.1):
  print('Solving for z = {}'.format(zz))

  # Call CAMB.
  pars           = camb.CAMBparams()

  pars.InitPower.set_params(As=As, ns=ns, r=r)
  pars.set_matter_power(redshifts=[zz], kmax=10.0)
  pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=omk, tau=tau)

  pars.NonLinear = model.NonLinear_none

  results        = camb.get_results(pars)
  results.calc_power_spectra(pars)

  trans           = results.get_matter_transfer_data()

  kh              = trans.transfer_data[0,:,0]  ##  trans.transfer_z('k/h');  ##  [(Mpc/h)^-1].
  k               = kh * results.Params.h

  # model.Transfer_cdm = 2
  transfer_cdm    = trans.transfer_data[2-1,:,0]
  transfer_cdm   /= transfer_cdm[0]

  transfer_bar    = trans.transfer_data[3-1,:,0]
  transfer_bar   /= transfer_bar[0]
  
  transfer_phot   = trans.transfer_data[4-1,:,0]
  transfer_phot  /= transfer_phot[0]

  # Massless.
  transfer_nut    = trans.transfer_data[5-1,:,0] 
  transfer_nut   /= transfer_nut[0]

  # Massive
  transfer_nu     = trans.transfer_data[6-1,:,0]
  transfer_nu    /= transfer_nu[0]

  transfer_tot    = trans.transfer_data[model.Transfer_tot-1,:,0] ##  trans.transfer_z('delta_tot') 
  transfer_tot   /= transfer_tot[0]

  primordial      = results.Params.scalar_power(k)

  matter          = primordial * transfer_tot**2 * k**4 / (k**3/(2*np.pi**2))

  # *********** ----  Pk has units of [Mpc^3]  ---- ***********.
  _, zs, linpk    = results.get_linear_matter_power_spectrum(hubble_units=False)

  ##
  ##  DH          = 16. / Om / h  ## [Mpc/h]

  #  for x in [transfer_cdm, transfer_bar, transfer_phot, transfer_nut, transfer_nu, transfer_tot]:
  #    plt.loglog(kh, x) 

  # lim             = 1. / kh / kh
  # lim            /= lim[-1]
  # lim            *= (transfer_cdm[-1] / transfer_cdm[0])

  # plt.loglog(kh, lim)

  # plt.loglog(kh, matter, label='Transfer')
  # plt.loglog(kh, Pk[0,:], '--', label='Direct');
  # plt.xlabel(r'$k\, [h Mpc^{-1}]$');
  # pl.legend()
  # pl.show()
  
  # [h/Mpc], [k].
  result          = np.vstack((kh, k, primordial, transfer_cdm, transfer_bar, transfer_phot, transfer_nut, transfer_nu, transfer_tot, linpk[0,:]))
  result          = result.T
  
  hdr             = 'k/h [h/Mpc]'.ljust(23) + 'k [1/Mpc]'.ljust(23) + 'Primoridal'.ljust(23) + ''.join(['T_{}(k)'.format(xx).ljust(25) for xx in ['cdm', 'bar', 'phot', 'massless neutrions', 'massive neutrinos', 'total']])
  hdr            += 'Linear Pk [(Mpc)^-3]'.ljust(23)

  np.savetxt('dat/Transfers_z{}.txt'.format(np.int(100. * zz)), result, header=hdr)

