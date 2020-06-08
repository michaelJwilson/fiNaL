import  os
import  glob
import  numpy              as      np
import  pylab              as      pl

from    params             import  *
from    growth             import  Dz
from    linearPk           import  Plin
from    get_yaml           import  get_yaml
from    fnl                import  pk_fnl
from    scipy.interpolate  import  interp1d
from    scipy.integrate    import  quad
from    collections        import  OrderedDict


redshift     = 2.0

all_samples  = glob.glob('samples/{:.1f}/sample_*'.format(redshift))

# Sorting by creation date restores mass ordering. 
all_samples.sort(key=os.path.getmtime)

kmax         = 0.1
kmins        = [0.001, 0.005, 0.01]

for ii, kmin in enumerate(kmins):
 for pair in all_samples:
  samples    = get_yaml(pair)

  def addtracer(tracer, tlabel=''):
    for k, v in samples[tracer].items():
      gg = globals()

      gg.update({k+'{}'.format(tlabel): v})

  tracers  = ['lomass', 'himass']

  for i, tracer in enumerate(tracers):
    addtracer(tracer, i+1)

  ks, Ps   = Plin(z1)

  P1       = b1 * b1 * Ps
  P2       = b2 * b2 * Ps

  alpha    = b1 / b2 

  X1       = 1. / (nz1 * P2)
  X2       = 1. / (nz2 * P2)
  
  print('\n\n-----------------------------------------------------------')
  print('Masses assumed [1.e10 Mo]: \t M1= {:.2f}, \t\t M2={:.2f}'.format(M1, M2))
  print('Linear biases assumed: \t\t b1= {:.2f}, \t\t b2={:.2f} \t\t (alpha = {:.2f})'.format(b1, b2, alpha))
  print('Number densities assumed:  \t nz1= {:.6e}, \t nz2={:.6e}\n'.format(nz1, nz2))

  # pl.loglog(ks, X1, label='X1')
  # pl.loglog(ks, X2, label='X2')

  # pl.xlabel(r'$k [h/{\rm Mpc}]$')
  # pl.legend(frameon=False)
  # pl.savefig('plots/Xs.png')

  r        = b1 * b2 * Ps / np.sqrt(P1 + 1. / nz1) / np.sqrt(P2 + 1. / nz2)

  #  Perfect cross-correlation.
  r        = np.ones_like(r)
  
  # pl.clf()

  # pl.loglog(ks, r, label='')

  # pl.xlabel(r'$k [h/{\rm Mpc}]$')
  # pl.ylabel(r'$r$')
  # pl.savefig('plots/r.png')

  ##
  num        = (alpha * alpha * X2 + X1 )**2.
  
  ## 
  num        = alpha * alpha * X2 * (1. + 2. * X2) + r * r * X1 * (1. + X2) + alpha * alpha * (1. - r * r) * (2. - r * r + 3. * X2)

  den        = alpha * alpha * (1. - r * r) + alpha * alpha * X2 + X1 + X1 * X2
  den        = den**2.

  Faa        = num / den

  # Error on alpha = (b1 / b2).
  sig2       = 1. / Faa
  sig        = np.sqrt(sig2)

  db1        = pk_fnl(b1, z1, fnl)
  db2        = pk_fnl(b2, z2, fnl)

  integrand  = (db1 / b1) - (db2 / b2)
  integrand *= alpha * ks
  integrand  = integrand**2.
  integrand *= Faa

  # pl.clf()
  # pl.loglog(ks, integrand, '-', c='k')
  # pl.show()

  integrand  = interp1d(ks, integrand, bounds_error=True)

  res, err   = quad(integrand, kmin, kmax)
  res       *= V / (2. * np.pi**2.)
 
  sigf       = 1. / np.sqrt(res)
 
  print('Volume \t\t kmin \t\t kmax \t\t Error on fnl.')
  print('{:e} \t {:e} \t {:e} \t {:e}'.format(V, kmin, kmax, sigf))

  results              = {}
  results['redshift']  = redshift
  results['kmin']      = kmin
  results['kmax']      = kmax
  results['Volume']    = V
  results['alpha']     = alpha
  results['sigma_fnl'] = np.float(sigf)
  
  results.update(samples)
  
  __                   = np.int(pair.split('_')[-1].replace('.yml', ''))
  
  with open('fishers/{}/fisher_{:d}.yml'.format(redshift, __ + ii * len(all_samples)), 'w') as outfile:
      yaml.dump(results, outfile, default_flow_style=False)
  
print('\n\nDone.\n\n')
