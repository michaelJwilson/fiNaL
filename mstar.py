import numpy as np


def K(z, alpha=0.0):
    return  2.5*(alpha - 1.0) * np.log10(1.+z)

def mu(z, chi, h=0.676):
    return 5. * np.log10((1. + z) * chi * 1.e6 / h / 10.)

def mstar(Mstar, z, chi, h=0.676, alpha=0.0):
    return  Mstar + mu(z, chi) + K(z, alpha) 

zs   = np.array([1.9, 2.8, 3.8])
Ms   = np.array([-19.68, -20.20, -20.71])
chis = np.array([3564.24, 4368.20, 4991.67])

for z, M, chi in zip(zs, Ms, chis):
    kcorr  = K(z)  
    mm     = mu(z, chi)
    result = mstar(M, z, chi)

    print('{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(z, M, mm, kcorr, result))
