from coffee.swsh import w3j
from coffee.swsh import clebschgordan as cg
from coffee.swsh import spinsfastpy as sfpy
import numpy as np
import math
from numpy import typeDict
import cmath


# Set up
w3j_calc = w3j.W3jBen("test.hdf")
sfpy.set_clebsch_gordan_default(cg.CGW3j(w3j_calc))

Ntheta = 11
Nphi = 11
spins = be.array([0])
lmax = 1
dtheta = math.pi / (Ntheta - 1)
dphi = 2 * math.pi / Nphi
delta_method = "RISBO_PRECOMPUTE"
Y000 = be.empty(Ntheta * Nphi * spins.shape[0], dtype=typeDict["complex128"])
Y01m1 = be.empty(Ntheta * Nphi * spins.shape[0], dtype=typeDict["complex128"])
Y010 = be.empty(Ntheta * Nphi * spins.shape[0], dtype=typeDict["complex128"])
Y011 = be.empty(Ntheta * Nphi * spins.shape[0], dtype=typeDict["complex128"])

for i in range(Ntheta):
    for j in range(Nphi):
        index = (i * Nphi) + j
        Y000[index] = 0.5 * math.sqrt(1 / math.pi)
        Y01m1[index] = 0.5 * (
            math.sqrt(3 / (2 * math.pi))
            * math.sin(i * dtheta)
            * cmath.exp(-complex(0, 1) * j * dphi)
        )
        Y010[index] = 0.5 * (math.sqrt(3 / (math.pi)) * math.cos(i * dtheta))
        Y011[index] = -0.5 * (
            math.sqrt(3 / (2 * math.pi))
            * math.sin(i * dtheta)
            * cmath.exp(complex(0, 1) * j * dphi)
        )

# Calculation

f = be.empty((spins.shape[0], Ntheta, Nphi), dtype=typeDict["complex128"])
f[0] = (Y000 + Y01m1).reshape(Ntheta, Nphi)
# print "f is %s"%repr(a)
salm = sfpy.forward(f, spins, lmax)
# salm.bandlimit_mult = True
# print "f in a spectral decomposition is"
# print salm

# print "f*f in a spectral decomposition is"
# print salm * salm

# print "f + f in a spectral decomposition is"
# g = salm + salm
# print g

# print "f - f in a spectral decomposition is"
# print salm  - salm

f2 = sfpy.backward(salm, Ntheta, Nphi)
# print "f2 is %s"%repr(f2)
# print "f + f %s"%repr(f + f)
print("f-f2 is %s" % repr(f - f2))
print("error is %s" % repr(be.sum(be.absolute(f - f2))))
