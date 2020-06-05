import numpy             as np
import pylab             as pl
import matplotlib.pyplot as plt

from   scipy.integrate   import quad


alphas = np.arange(-2.0, -1.0, 0.2)

# x = L / L*
xmins  = np.arange(0.1, 10., 0.1) 

def integrand(x, alpha):
    return  (x**alpha) * np.exp(-x) 

output = [xmins]

for alpha in alphas:
    # y = Phi / Phi*
    ys = []

    for xmin in xmins:
      ys.append(quad(integrand, xmin,  np.inf, args=(alpha))[0])

    output.append(ys)
      
    ys = np.array(ys)

    pl.semilogy(xmins, ys, label=r'$\alpha={:.1f}$'.format(alpha))


output = np.array(output).T

header  = 'L / L*'.ljust(23)
header += ''.join(['Phi(alpha={:.1f})/Phi*'.format(a).ljust(25) for a in alphas]) 

np.savetxt('dat/schechters.txt', output, header=header)
    
pl.xlabel(r'$L_{\rm lim} / L_*$')
pl.ylabel(r'$\tilde \Phi / \Phi_*$')

pl.xlim(0.1, 5.0)
pl.ylim(1.e-3, 10.)

pl.legend(frameon=False)

plt.tight_layout()

pl.savefig('plots/schechters.png')
