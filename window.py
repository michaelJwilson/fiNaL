import  numpy             as  np
import  pylab             as  pl
import  matplotlib.pyplot as  plt


Delta  = 20.
L      = 1000.
D      = 60.

def sinc(x):
    x /= np.pi
    
    return  np.sinc(x)

def W(k):
    W1 = sinc(k * Delta / 2.)
    W2 = sinc(k * L     / 2.)
    W3 = sinc(k * D     / 2.)
    
    return  W1 * W2 / W3  ## [Delta L / D].

def div(m):
    return  2. * m * np.pi / D

lim = 1.e3

ks  = np.arange(-lim, lim, 1.) / lim
Ws  = W(ks)

pl.plot(ks, Ws, alpha=0.3, c='k')
pl.plot(ks, sinc(ks * D / 2.), alpha=0.3, c='r', label='Sparsity')
pl.plot(ks, sinc(ks * L / 2.), alpha=0.3, c='c', label='Survey extent')
pl.plot(ks, sinc(ks * Delta / 2.), alpha=0.3, c='magenta', label='Sample extent')

for m in np.arange(1, 6, 1):
    pl.axvline(div(m), ymin=0., ymax=1., c='b', lw=0.5, alpha=0.5)

pl.axhline(0.0, xmin=0., xmax=1., c='k', lw=0.5)
    
pl.xlim(-0.01, 0.2)
pl.ylim(-0.30, 1.1)

pl.xlabel(r'$k \ \ [h {\rm Mpc}]$')
pl.ylabel(r'$\tilde{\mathcal{W}} (k) \ \ [\Delta L / D]$')

pl.legend(loc=1, frameon=False)

plt.tight_layout()

pl.savefig('wk.pdf')
