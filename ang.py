import pyccl
import pylab             as     pl
import numpy             as     np
import matplotlib.pyplot as     plt

from   desimodel.io                  import load_tiles
from   astropy.table                 import Table
from   nbodykit.source.catalog.array import ArrayCatalog
from   nbodykit.transform            import SkyToCartesian
from   nbodykit.cosmology            import Planck15
from   scipy.spatial.distance        import cdist
from   pyccl.background              import comoving_angular_distance, comoving_radial_distance


zz    = 2.0

cosmo = pyccl.Cosmology(Omega_c=0.25, Omega_b=0.05,
                        h=0.7, n_s=0.95, sigma8=0.8,
                        transfer_function='bbks')

# Comoving radial distance; Mpc.
chi   = comoving_radial_distance(cosmo, 1./(1+zz))
chi  *= cosmo['h']

# Comoving angular diameter distance (https://ccl.readthedocs.io/en/latest/_modules/pyccl/background.html#comoving_angular_distance).
DA    = comoving_angular_distance(cosmo, 1./(1+zz))
DA   *= cosmo['h']

tiles = load_tiles()
tiles = Table(tiles)
tiles = tiles[(tiles['PASS'] == 1) & (tiles['IN_DESI'] == 1)]

tiles['RA'][tiles['RA'] < 100.] += 360.

pl.scatter(tiles['RA'], tiles['DEC'], c=np.log2(tiles['EXPOSEFAC']), s=4)
pl.colorbar(label=r'log$_{2}$|EXPOSURE FACTOR|')
pl.xlim(450., 90.0)
pl.xlabel('Right ascension  [deg.]')
pl.ylabel('Declination  [deg.]')
plt.tight_layout()
pl.savefig('ang.pdf')

## 
tiles['Z'] = zz

##  NGC
tiles      = tiles[tiles['RA'] < 300.]

##  Horizontal strip at high dec.
tiles      = tiles[tiles['DEC'] > 60.]
tiles.sort('RA')

tiles      = ArrayCatalog(tiles) 

##  Mpc/h
pos        = SkyToCartesian(tiles['RA'], tiles['DEC'], tiles['Z'], Planck15, observer=[0, 0, 0], degrees=True, frame='icrs').compute()
rs         = cdist(pos, pos, metric='euclidean')

print(rs.max())
