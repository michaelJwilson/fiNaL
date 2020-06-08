import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    params             import  *
from    growth             import  Dz
from    linearPk           import  Plin
from    get_yaml           import  get_yaml


def dbk(b, z, printit=False):
    if z < 1.0:
        raise ValueError('z < 1.0 is unsupported.')

    close   = np.around(z, decimals=1)    
    dat     = np.loadtxt('dat/Transfers_z{}.txt'.format(np.int(100. * close)))

    if printit:
        print(z, close)
    
    ks      = dat[:,1]  # [1 . / Mpc]
    Ts      = dat[:,8]  # Total transfer function, normalised to unity on large scales.  
    
    return  3. * (b - p) * dc * Om * H0**2 / c / c / ks / ks / Ts / Dz(z)

def pk_fnl(b, z, fnl, printit=False):
    # CDM power spectrum.
    _, PP   = Plin(z, printit=printit)

    # Eqn. (5) of https://arxiv.org/pdf/0807.1770.pdf; requires linear DM power spectrum.
    return  PP * (b + fnl * dbk(b, z, printit=printit))**2.


if __name__ == '__main__':
    tracer  = 'GRUSH24'

    samples = get_yaml('dat/samples.yml')
    
    for k, v in samples[tracer].items():
        exec(k+'=v')

    ks, PP  = Plin(z)
    Pk      = pk_fnl(b, z, fnl, printit=True)
         
    pl.axhline(1. / nz, xmin=0.0, xmax=1.0, c='k', lw=0.1)

    pl.loglog(ks, Pk, label=r'${} (z, b, \bar n) = ({}, {}, {:.1e})$'.format(tracer, z, b, nz))
    pl.loglog(ks, b * b * PP) 
    pl.loglog(ks,         PP)

    pl.xlabel(r'$k \ [h/{\rm Mpc}]$')
    pl.ylabel(r'$P(k) \ [({\rm Mpc}/h)^3]$')

    pl.legend(frameon=False)

    pl.xlim(1.e-3, 0.2)
    pl.ylim(100., 1.e6)

    plt.tight_layout()

    pl.show()
    # pl.savefig('plots/pk_fnl.pdf')
    
    print('\n\nDone.\n\n')
