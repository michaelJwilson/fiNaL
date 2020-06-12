import glob
import pylab              as     pl
import numpy              as     np 
import matplotlib.pyplot  as     plt

from   get_yaml                              import get_yaml
from   mpl_toolkits.axes_grid1               import make_axes_locatable
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rc('font', family='serif')

plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')

zs         = [2.0, 3.0, 4.0]
kmins      = [0.01, 0.005, 0.001]

fig, axes  = plt.subplots(3, 3, figsize=(15, 15))

colors     = {1: 'k', 2: 'c'}

for ii, redshift in enumerate(zs):
 for jj, kmin in enumerate(kmins):
   for sampling in [1]:
    files      = glob.glob('fishers/{:.1f}/sampling_{:d}/kmin_{:d}/*.yml'.format(redshift, sampling, np.int(1000. * kmin)))
    nfile      = len(files)
    
    for isample in np.arange(0, nfile, 1):
      fpath    = 'fishers/{:.1f}/sampling_{:d}/kmin_{:d}/sample_{}.yml'.format(redshift, sampling, np.int(1000. * kmin), isample)     
      result   = get_yaml(fpath)
      
      Llo      = result['lomass']['L']
      sigf     = result['sigma_fnl']
      
      axes[ii,jj].scatter(Llo, sigf, c=colors[sampling], marker='.', s=1)

      axes[ii,jj].set_xlim(0.0, 5.0)
      axes[ii,jj].set_ylim(0.0, 5.0)

for jj, kmin in enumerate(kmins):
  axes[-1, jj].set_xlabel(r'$L_{\rm lim} / L_*$')

for ii, redshift in enumerate(zs):
  axes[ii, 0].set_ylabel(r'Error on $f_{\rm NL}$')
      
pl.savefig('plots/moneyplot.png')
