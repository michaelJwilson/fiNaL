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


subsample  = 1.0

zs         = [2.0, 3.0, 4.0]
kmins      = [0.01, 0.005, 0.001]

vmin       = 1.0
vmax       = 3.0

fig, axes  = plt.subplots(3, 3, figsize=(15, 15))

M1         = np.logspace(10., 13., 20) / 1.e10
M2         = np.logspace(10., 13., 20) / 1.e10
M1s, M2s   = np.meshgrid(M1, M2)

for ii, redshift in enumerate(zs):
  search   = 'fishers/{:.1f}/*_sub_{:d}.yml'.format(redshift, np.int(100. * subsample))
  results  = glob.glob(search)
  
  for jj, kmin in enumerate(kmins):
     ZZ    = np.zeros_like(M1s)

     for result in results:
          result = get_yaml(result)
    
          if result['kmin'] == kmin:
              M1     = result['himass']['M']
              M2     = result['lomass']['M']

              sigf   = result['sigma_fnl']

              ZZ[(M1s == M1) & (M2s == M2)] = sigf
              
              # print('{: 6.1f} \t {: 6.1f} \t {:.3f}'.format(M1, M2, sigf))
  
              # im   = axes[ii, jj].scatter(np.array([M1]), np.array([M2 / M1]), c=np.array([sigf]), vmin=vmin, vmax=vmax, s=10, cmap='coolwarm')
                            
              cax    = inset_axes(axes[ii, jj],
                                  width="50%",  # width = 50% of parent_bbox width
                                  height="5%",  # height : 5%
                                  loc='upper left')
              
              axes[ii, jj].set_xlim(0.9, 1001.00)
              axes[ii, jj].set_ylim(1.0, 1.01e-3)

              axes[ii, jj].set_xscale('log')
              axes[ii, jj].set_yscale('log')

              # divider = make_axes_locatable(axes[ii, jj])
              # cax     = divider.append_axes('top', pad=-0.1, size='2%')

              label     = r'$(z,k_{min})' + ' = ({:.1f}, {:.3f})$'.format(redshift, kmin)
              
              #cb        = fig.colorbar(im, cax=cax, orientation='horizontal', shrink=0.5)
              #cb.ax.tick_params(labelsize='small')
              #cb.set_label(label=label, fontsize='small', fontweight='normal')

     print(ZZ[ZZ > 0.0])

     CS  = axes[ii, jj].contourf(M1s, M2s, ZZ, 10, cmap='coolwarm', origin='lower')   
     CS2 = axes[ii, jj].contour(CS, levels=CS.levels[::2], colors='k', alpha=0.5, origin='lower')  
     
for jj, kmin in enumerate(kmins):              
  axes[-1, jj].set_xlabel(r'$M_1 \ [10^{10} \ M_{\odot}]$')

for ii, redshift in enumerate(zs):
  axes[ii, 0].set_ylabel(r'$M_2 \ [M1]$')

fig.suptitle('Sampling factor: {:.1f}'.format(subsample), y=1.02)
  
# plt.tight_layout()

# pl.show()
pl.savefig('plots/Fishers_sub_{:d}.png'.format(np.int(100. * subsample)))
