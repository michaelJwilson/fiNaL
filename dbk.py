import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    params             import  *
from    growth             import  Dz
from    linearPk           import  Plin
from    get_yaml           import  get_yaml


def dbk(b, z, printit=False, return_ks=False):
    if z < 1.0:
        raise ValueError('z < 1.0 is unsupported.')

    close   = np.around(z, decimals=1)
    dat     = np.loadtxt('dat/Transfers_z{}.txt'.format(np.int(100. * close)))

    if printit:
        print(z, close)

    ks      = dat[:,1]  # [1 . / Mpc]                                                                                                                                                                                           
    Ts      = dat[:,8]  # Total transfer function, normalised to unity on large scales.                                                                                                                                         

    if return_ks:
       return  ks, 3. * (b - p) * dc * Om * H0**2 / c / c / ks / ks / Ts / Dz(z)
        
    else:
       return  3. * (b - p) * dc * Om * H0**2 / c / c / ks / ks / Ts / Dz(z)


if __name__ == '__main__':
    b1       = 2.0
    b2       = 1.4

    ks, dbk1 = dbk(b1, 1.0, return_ks=True)
    ks, dbk2 = dbk(b2, 1.0, return_ks=True)

    den      = dbk1 / b1 - dbk2 - b2
    
    pl.loglog(ks, dbk1 / b1, label='dbk1')
    pl.loglog(ks, dbk2 / b2, label='dbk2')
    
    # pl.loglog(ks, np.abs(den))

    pl.xlim(1.e-5, 1.e-2)
    
    pl.legend(frameon=False)
    pl.show()
