#!/usr/bin/env python
# encoding: utf-8
"""
This module calculates the quadrature weights.
"""

from coffee.settings import be
import scipy.fft as fft

def indx2p(ip, wsize):
    return ip - wsize if ip > wsize / 2 else ip

def Compute_quadrature_weights(W, wsize):
    w = be.zeros(wsize, dtype=be.complex128)
    
    for ip in range(wsize):
        p = indx2p(ip, wsize)
        eo = abs(p % 2)
        
        if p == -1:
            w[ip] = complex(0, be.pi/2)
        elif p == 1:
            w[ip] = complex(0, -be.pi/2)
        elif eo == 0:
            w[ip] = 2.0 / (1.0 - p*p)
        else:
            w[ip] = 0
    
    w = fft.ifft(w, norm='forward')
    W[:] = w[:]

def Precompute_quadrature_weights(Ntheta):
    Sampling_on_torus = 2 * (Ntheta - 1)
    W = be.zeros(Sampling_on_torus, dtype=be.complex128)
    Compute_quadrature_weights(W, Sampling_on_torus)
    return W