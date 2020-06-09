import re
import yaml
import numpy as np
import itertools

from   scipy.interpolate import interp1d


Ms       = np.logspace(10., 13., 20)
pairs    = list(itertools.combinations(Ms, 2))

redshift = 4.0

nbars    = np.loadtxt('dat/halo_nbar.txt')

infile   = open('dat/halo_nbar.txt', 'r')
hdr      = infile.readline()

zs       = re.findall("\d+\.\d+", hdr)
zs       = np.array(zs).astype(np.float)

biases   = np.loadtxt('dat/halo_biases.txt')

cut      = zs == redshift

idx      = np.where(cut)[0]
idx     += 1

#
Ms       = nbars[:,0]
ns       = nbars[:,idx]
bs       = biases[:,idx]

ns       = interp1d(Ms.T, ns.T, bounds_error=True)
bs       = interp1d(Ms.T, bs.T, bounds_error=True)

print('Solving for:')

for i, pair in enumerate(pairs):
    M1     = np.float(pair[0])
    M2     = np.float(pair[1])

    print('M1: {} \t M2: {}'.format(M1 / 1.e10, M2 / 1.e10))
    
    b1     = np.float(bs(M1))
    b2     = np.float(bs(M2))

    nz1    = np.float(ns(M1))
    nz2    = np.float(ns(M2))

    z1     = redshift
    z2     = redshift
    
    sample = {'lomass':  {'M': M1 / 1.e10, 'b': b1, 'z': redshift, 'nz': nz1}, 'himass': {'M': M2 / 1.e10, 'b': b2, 'z': redshift, 'nz': nz2}}
    
    with open('samples/{}/sample_{}.yml'.format(redshift, i), 'w') as outfile:
      yaml.dump(sample, outfile, default_flow_style=False)
