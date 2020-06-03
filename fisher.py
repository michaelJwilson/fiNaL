import numpy as np


b1    = 6.60
b2    = 4.02

n1    = 1.1e-5
n2    = 6.0e-4

Pm    = 1.e4

P1    = b1 * b1 * Pm
P2    = b2 * b2 * Pm

#
alpha = b1 / b2 

X1    = 1. / (n1 * P2)
X2    = 1. / (n2 * P2)

r     = b2 * b1 * Pm / np.sqrt((P1 + 1. / n1) * (P2 + 1. / n2))

num   = alpha * alpha * X2 * (1. + 2. * X2) + r * r * X1 * (1. + X2) + alpha * alpha * (1. - r * r) * (2. - r * r + 3. * X2)

den   = alpha * alpha * (1. - r * r) + alpha * alpha * X2 + X1 + X1 * X2
den   = den**2.

Faa   = num / den

# Error on alpha = (b2 / b1)
sig2  = 1. / Faa
sig   = np.sqrt(sig2)

print(sig)
