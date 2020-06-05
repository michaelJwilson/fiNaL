import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    params             import  *
from    growth             import  Dz


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

ks, k2T        = np.loadtxt('k2Transfers.txt', unpack=True)

def dbk(b, z):
    return  3. * (b - p) * dc * Om * H0**2 / c /c / k2T / Dz(zz)

def pk_fnl(b, z, fnl):
    # Eqn. (5) of https://arxiv.org/pdf/0807.1770.pdf; requires linear DM power spectrum.
    return  Plin(z)  * (b + dbk(b,z) * fnl)**2.


if __name__ == '__main__':
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

    print('\n\nDone.\n\n')
