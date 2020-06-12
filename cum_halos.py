import re
import pylab             as     pl 
import numpy             as     np
import pandas            as     pd
import matplotlib.pyplot as     plt

from   scipy.integrate   import quad
from   scipy.interpolate import interp1d


infile    = open('dat/halo_nbar.txt', 'r')
hdr       = infile.readline()

zs        = re.findall("\d+\.\d+", hdr)
zs        = np.array(zs).astype(np.float)

# Solar masses, dn/dlnM(z) [(Mpc/h)^3].
names     = ['Ms'] + zs.tolist()

dat       = pd.read_csv('dat/halo_nbar.txt',   sep='\s+', comment ='#', names=names)
bat       = pd.read_csv('dat/halo_biases.txt', sep='\s+', comment ='#', names=names)

lnMs      = np.log(dat['Ms'])

fig, ax1  = plt.subplots()
ax2       = ax1.twinx()

colors    = plt.rcParams['axes.prop_cycle'].by_key()['color']

result_bs = [dat['Ms']]
result_ns = [dat['Ms']]

phi_stars = {2.0: 0.02327, 3.0: 0.01762, 4.0: 0.00666}

for i,  z in enumerate(zs):
  # dn/dln M [(Mpc/h)^3].

  integrand  = interp1d(lnMs, dat[z],          bounds_error=True)
  bintegrand = interp1d(lnMs, dat[z] * bat[z], bounds_error=True)
  
  beffs      = []
  ntildes    = []

  for lnM in lnMs:
    beffs.append(quad(bintegrand,  lnM,  np.log(1.e12), limit=5000)[0])
    ntildes.append(quad(integrand, lnM,  np.log(1.e12), limit=5000)[0])

  result_bs.append(beffs)
  result_ns.append(ntildes)
    
  beffs      = np.array(beffs)
  ntildes    = np.array(ntildes)

  beffs     /= ntildes

  if i > 0:
    # No PhiStar for z=1.0 in Parsa.
    ax1.loglog(dat['Ms'],   ntildes / phi_stars[z], label=r'$z={:.1f}$'.format(z), c=colors[i])

  ax2.semilogx(dat['Ms'], beffs, label=r'$z={:.1f}$'.format(z), alpha=0.5, c=colors[i])

result_bs = np.array(result_bs).T
result_ns = np.array(result_ns).T
  
header    = 'M [Msun]'.ljust(23)
header   += ''.join(['n({:.1f}) [(Mpc/h)^-3]'.format(zz).ljust(25) for zz in [1.0, 2.0, 3.0, 4.0]])

np.savetxt('dat/halo_lim_nbar.txt', result_ns,  header=header)

header  = 'M [Msun]'.ljust(23)
header += ''.join(['b(z={:.1f})'.format(zz).ljust(25) for zz in [1.0, 2.0, 3.0, 4.0]])

np.savetxt('dat/halo_eff_biases.txt', result_bs, header=header)

ax1.set_xlim(1.e10, 1.e13)
ax1.set_ylim(1.e-3, 100.0)
  
ax1.set_xlabel(r'${\rm Limiting} \ M_h \ [M_\odot]$')

# [({\rm Mpc}/h)^-3]
ax1.set_ylabel(r'$\tilde n \ [\Phi_*]$')  

ax1.legend(frameon=False)

ax2.set_ylabel(r'$b_{\rm eff}(M_h)$')
ax2.set_ylim(0.0, 15.0)

ax2.axhline(y=1.0, c='k', lw=0.5)

plt.tight_layout()

pl.savefig('plots/ehalos.pdf')
