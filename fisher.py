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
for i, tracer in enumerate(['GRUSH24', 'GRUSH25.5']):
    addtracer(tracer, i+1)

    
ks, Ps     = Plin(z1)

P1         = b1 * b1 * Ps
P2         = b2 * b2 * Ps

#
alpha      = b1 / b2 

X1         = 1. / (nz1 * P2)
X2         = 1. / (nz2 * P2)

r          = b2 * b1 * Ps / np.sqrt((P1 + 1. / nz1) * (P2 + 1. / nz2))

num        = alpha * alpha * X2 * (1. + 2. * X2) + r * r * X1 * (1. + X2) + alpha * alpha * (1. - r * r) * (2. - r * r + 3. * X2)

den        = alpha * alpha * (1. - r * r) + alpha * alpha * X2 + X1 + X1 * X2
den        = den**2.

Faa        = num / den

# Error on alpha = (b2 / b1)
sig2       = 1. / Faa
sig        = np.sqrt(sig2)

# pl.loglog(ks, 100. * sig / alpha, '.', markersize=1)                                                                                                                                                                              
# pl.show() 

db1        = pk_fnl(b1, z1, fnl)
db2        = pk_fnl(b2, z2, fnl)

integrand  = (db1 / b1) - (db2 / b2)
integrand *= alpha * ks
integrand  = integrand**2.
integrand *= Faa

integrand  = interp1d(ks, integrand, bounds_error=True)

kmin       = 0.01  # ks.min()
kmax       = 0.1   # ks.max()

res, err   = quad(integrand, kmin, kmax)
res       *= V / (2. * np.pi**2.)

sigf       = 1. / np.sqrt(res)

# pl.loglog(ks, integrand, '-', c='k')                                                                                                                                                                          
# pl.show()

print(V, kmin, kmax, X1, X2, res)
