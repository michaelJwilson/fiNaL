import pylab as pl
import fitsio


dat = fitsio.read('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randoms/randoms-inside-dr8-0.31.0-2.fits')

print(dat)
