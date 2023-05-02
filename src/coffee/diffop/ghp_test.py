from builtins import range
import numpy as np
from numpy import typeDict
import sys
import math
import cmath

from . import ghp
from coffee.swsh import spinsfastpy as sfpy


def compute(power):
    Ntheta = 2 ** (power)
    Nphi = 2**power
    spins = be.array([0])
    lmax = 3
    dtheta = math.pi / (Ntheta - 1)
    dphi = 2 * math.pi / Nphi

    Y000 = be.empty((Ntheta, Nphi), dtype=typeDict["complex"])
    Y01m1 = be.empty_like(Y000)
    Y010 = be.empty_like(Y000)
    Y011 = be.empty_like(Y000)
    Y110 = be.empty_like(Y000)
    Y11m1 = be.empty_like(Y000)
    Y111 = be.empty_like(Y000)

    eth = ghp.eth()

    for i in range(Ntheta):
        for j in range(Nphi):
            Y000[i, j] = 0.5 * math.sqrt(1 / math.pi)

            Y01m1[i, j] = 0.5 * (
                math.sqrt(3 / (2 * math.pi))
                * math.sin(i * dtheta)
                * cmath.exp(-complex(0, 1) * j * dphi)
            )
            Y010[i, j] = 0.5 * (math.sqrt(3 / (math.pi)) * math.cos(i * dtheta))
            Y011[i, j] = -0.5 * (
                math.sqrt(3 / (2 * math.pi))
                * math.sin(i * dtheta)
                * cmath.exp(complex(0, 1) * j * dphi)
            )

            Y110[i, j] = math.sqrt(3 / 8 * math.pi) * math.sin(i * dtheta)
            Y11m1[i, j] = (
                -0.25
                * math.sqrt(3 / math.pi)
                * (1 + math.cos(i * dtheta))
                * cmath.exp(-complex(0, 1) * j * dphi)
            )
            Y111[i, j] = (
                -math.sqrt(3 / 16 * math.pi)
                * (1 - math.cos(i * dtheta))
                * cmath.exp(complex(0, 1) * j * dphi)
            )

    salm_Y01m1 = sfpy.forward(Y01m1, be.array([0]), lmax)
    # print "salm_Y01m1 = %s"%salm_Y01m1
    eth_Y01m1 = eth(Y01m1, spins, lmax)
    # print eth_Y01m1
    eth(eth_Y01m1, spins + 1, lmax)


#    salm_eth_Y01m1 = sfpy.forward(eth_Y01m1, be.array([1]), lmax)
#    #print "salm_eth_Y01m1 = %s"%sfpy.forward(eth_Y01m1, be.array([1]), lmax)
#    salm_Y11m1 = sfpy.forward(Y11m1, be.array([1]), lmax)
#    #print "math.sqrt(2)*salm_Y11m1 = %s"%(math.sqrt(2)*salm_Y11m1)
#    #print Y01m1.shape
#    #print eth_Y01m1.shape
#    fact = math.sqrt(2)
#    raw_error = be.abs( eth_Y01m1 - fact*Y11m1 )
#    print power
#    print be.max(raw_error)
#    #print be.sum(raw_error) * 2 * math.pi**2 / (Ntheta * Nphi)
#    salm_raw_error = be.abs(salm_eth_Y01m1 - fact * salm_Y11m1)
#    print be.max(salm_raw_error.view(be.ndarray))
#    print ""
#    #print be.sum(salm_raw_error.view(be.ndarray))

for i in range(4, 5):
    compute(i)
