import  numpy              as      np
import  pylab              as      pl

from    params             import  *
from    growth             import  Dz
from    linearPk           import  Plin
from    samples            import  samples
from    fnl                import  pk_fnl
from    scipy.interpolate  import  interp1d
from     scipy.integrate   import  quad


def addtracer(tracer, tlabel=''):
  for k, v in samples[tracer].items():
    gg = globals()

    gg.update({k+'{}'.format(tlabel): v})

# Update global variables with b1, b2, z1, z2, nz1, nz2 etc.
tracers    = ['GRUSH24', 'GRUSH25.5']
# tracers  = ['ebossQSO', 'ebossMatter']

for i, tracer in enumerate(tracers):
    addtracer(tracer, i+1)
    
ks, Ps     = Plin(z1)

P1         = b1 * b1 * Ps
P2         = b2 * b2 * Ps

#
alpha      = b1 / b2 

print('Linear biases assumed: b1= {}, b2={} (alpha = {})'.format(b1, b2, alpha))

X1         = 1. / (nz1 * P2)
X2         = 1. / (nz2 * P2)

print('Number densities assumed: nz1= {}, nz2={}'.format(nz1, nz2))

pl.loglog(ks, X1, label='X1')
pl.loglog(ks, X2, label='X2')

pl.xlabel(r'$k [h/{\rm Mpc}]$')
pl.legend(frameon=False)
pl.savefig('plots/Xs.png')

r          = b1 * b2 * Ps / np.sqrt(P1 + 1. / nz1) / np.sqrt(P2 + 1. / nz2)

pl.clf()

pl.loglog(ks, r)

pl.xlabel(r'$k [h/{\rm Mpc}]$')
pl.ylabel(r'$r$')

pl.legend(frameon=False)
pl.savefig('plots/r.png')

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

kmin       = 0.002  # ks.min()
kmax       = 0.1    # ks.max()

res, err   = quad(integrand, kmin, kmax)
res       *= V / (2. * np.pi**2.)

sigf       = 1. / np.sqrt(res)

print('\n\nVolume \t kmin \t kmax \t Error on fnl.')
print('{:e} \t {:e} \t {:e} \t {:e}'.format(V, kmin, kmax, sigf))

print('\n\nDone.\n\n')
