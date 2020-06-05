import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    params             import  *
from    growth             import  Dz
from    linearPk           import  Plin
from    samples            import  samples


ks, k2T = np.loadtxt('dat/k2Transfers.txt', unpack=True)

def dbk(b, z):
    return  3. * (b - p) * dc * Om * H0**2 / c /c / k2T / Dz(z)

def pk_fnl(b, z, fnl):
    _, PP = Plin(z)

    # Eqn. (5) of https://arxiv.org/pdf/0807.1770.pdf; requires linear DM power spectrum.
    return  PP * (b + dbk(b,z) * fnl)**2.


if __name__ == '__main__':
    tracer = 'GRUSH24'

    for k, v in samples[tracer].items():
        exec(k+'=v')

    _, PP = Plin(z)
    Pk    = pk_fnl(b, z, fnl)
        
    pl.axhline(1. / nz, xmin=0.0, xmax=1.0, c='k', lw=0.1)

    pl.loglog(ks, Pk, label=r'${} (z, b, \bar n) = ({}, {}, {:.1e})$'.format(tracer, z, b, nz))
    pl.loglog(ks, b * b * PP) 
    pl.loglog(ks,         PP)

    pl.xlabel(r'$k \ [h/{\rm Mpc}]$')
    pl.ylabel(r'$P(k)$')

    pl.legend(frameon=False)

    pl.xlim(1.e-3, 0.2)
    pl.ylim(100., 1.e6)

    plt.tight_layout()

    pl.savefig('plots/pk_fnl.pdf')
    
    print('\n\nDone.\n\n')
