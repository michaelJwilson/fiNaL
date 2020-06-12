import re
import yaml
import numpy as np
import itertools

from   scipy.interpolate import interp1d
from   get_schechters    import get_schechter


redshift = 2.0

Ms       = np.logspace(10., 12., 20)
pairs    = list(itertools.combinations(Ms, 2))

nbars    = np.loadtxt('dat/halo_lim_nbar.txt')

infile   = open('dat/halo_lim_nbar.txt', 'r')
hdr      = infile.readline()

zs       = re.findall("\d+\.\d+", hdr)
zs       = np.array(zs).astype(np.float)

biases   = np.loadtxt('dat/halo_eff_biases.txt')

cut      = zs == redshift

idx      = np.where(cut)[0]
idx     += 1

#
Ms       = nbars[:,0]
ns       = nbars[:,idx]
bs       = biases[:,idx]

ns       = interp1d(Ms.T, ns.T, bounds_error=True)
bs       = interp1d(Ms.T, bs.T, bounds_error=True)

deets    = {2.0: {'PhiStar': 0.02327, 'alpha': -1.32}, 3.0: {'PhiStar': 0.01762, 'alpha': -1.31}, 4.0: {'PhiStar': 0.00666, 'alpha': -1.43}}

PhiStar  = deets[redshift]['PhiStar']
alpha    = deets[redshift]['alpha']

# returns L / L* given Phi / Phi_*
xs       = get_schechter(redshift, reverse=True)

print('Solving for:')

for i, pair in enumerate(pairs):
    M1     = np.float(pair[0])
    M2     = np.float(pair[1])
    
    b1     = np.float(bs(M1))
    b2     = np.float(bs(M2))

    nz1    = np.float(ns(M1))
    nz2    = np.float(ns(M2))

    phi1   = np.float(nz1 / PhiStar)
    phi2   = np.float(nz2 / PhiStar)

    good   = 1 
    
    try:
      L1   = np.float(xs(phi1))

    except:
      L1   = np.inf
      good = 0
      
    try:
      L2   = np.float(xs(phi2))

    except:
      L2   = np.inf
      good = 0
      
    print('Sample {: 3d} \t M1: {:.5f} \t M2: {:.5f} \t L1: {:.5f} \t L2: {:.5f} \t {:d}'.format(i, M1 / 1.e10, M2 / 1.e10, L1, L2, good))
      
    t1     = np.float(1. / L1**2.)
    t2     = np.float(1. / L2**2.)
    
    z1     = redshift
    z2     = redshift
    
    sample = {'lomass':  {'M': M1 / 1.e10, 'b': b1, 'z': redshift, 'nz': nz1, 'phit': phi1, 'L': L1, 't': t1}, 'himass': {'M': M2 / 1.e10, 'b': b2, 'z': redshift, 'nz': nz2, 'phit': phi2, 'L': L2, 't': t2}}

    print()
    print(sample['himass'])
    print(sample['lomass'])
    
    with open('samples/{}/sampling_1/sample_{}.yml'.format(redshift, i), 'w') as outfile:
      yaml.dump(sample, outfile, default_flow_style=False)
