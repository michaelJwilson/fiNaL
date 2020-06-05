import numpy             as np
import pylab             as pl
import matplotlib.pyplot as plt

from   scipy.integrate   import quad


alphas = np.arange(-2.0, -1.0, 0.2)

# x = L / L*
xmins  = np.arange(0.05, 10., 0.05) 

for x in np.arange(1, 6, 1):
    print('{:.4f} \t {:.4f}'.format(x, -2.5 * np.log10(x)))

def integrand(x, alpha):
    return  (x**alpha) * np.exp(-x) 

output = [xmins]

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

base   = quad(integrand, 1.0,  np.inf, args=(-1.4))[0]

for t in [1.0, 0.5, 0.25, 0.1, 0.05]:
    lvl = quad(integrand, t**-0.5,  np.inf, args=(-1.4))[0]
    
    pl.axhline(lvl, c='k', alpha=0.5, lw=0.4)

    plt.text(4.9, 1.05 * lvl, '{:.1f}'.format(np.around(base / lvl, 1)), horizontalalignment='right')
               
for alpha, color in zip(alphas, colors):
    # y = Phi / Phi*
    ys   = []
    
    for xmin in xmins:
      ys.append(quad(integrand, xmin,  np.inf, args=(alpha))[0])

    output.append(ys)
      
    ys = np.array(ys)

    pl.semilogy(xmins, ys, label=r'$\alpha={:.1f}$'.format(alpha))

output  = np.array(output).T

header  = 'L / L*'.ljust(23)
header += ''.join(['Phi(alpha={:.1f})/Phi*'.format(a).ljust(25) for a in alphas]) 

np.savetxt('dat/schechters.txt', output, header=header)

pl.xlabel(r'$L_{\rm lim} / L_*$')
pl.ylabel(r'$\tilde \Phi / \Phi_*$')

pl.xlim(0.05, 5.0)
pl.ylim(5.e-4, 1.)

pl.legend(frameon=False, loc=1, ncol=3)

ax      = pl.gca()
ax2     = ax.twiny()

ax2.set_xlim(ax.get_xlim())

ts      = np.array([0.05, 0.1, 0.25, 0.5, 1.0])
xs      = ts**-0.5

ax2.set_xticks(xs)
ax2.set_xticklabels(['{:.2f}'.format(t) for t in ts])

ax2.set_xlabel(r'$t_{\rm lim} / t_*$')

plt.tight_layout()

pl.savefig('plots/schechters.png')
