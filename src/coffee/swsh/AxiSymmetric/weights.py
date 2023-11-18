#!/usr/bin/env python
# encoding: utf-8
"""
This module calculates the quadrature weights.
"""

from coffee.settings import be
import scipy.fft as fft

def indx2p(ip, wsize):
    return ip - wsize if ip > wsize/2 else ip

def compute_quadrature_weights(wsize):
    w = be.zeros(wsize, dtype=be.complex128)
    for ip in range(wsize):
        p = indx2p(ip, wsize)
        eo = abs(p % 2)
        if p == -1:
            w[ip] = 1j * be.pi/2.
        elif p == 1:
            w[ip] = -1j * be.pi/2.
        elif eo == 0:
            w[ip] = 2./(1.-p*p)
        else:
            w[ip] = 0
    W = fft.ifft(w, norm='forward')
    return W

def precompute_quadrature_weights(Ntheta):
    Sampling_on_torus = 2*(Ntheta-1)
    return compute_quadrature_weights(Sampling_on_torus)