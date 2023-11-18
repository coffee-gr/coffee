#!/usr/bin/env python
# encoding: utf-8
"""
This module calculates the axi-symmetric forward transform.
"""

from coffee.settings import be
import scipy.fft as fft

def Wn(ip, size):
    index = abs(ip)
    if ip == -1:
        w = 1j * be.pi/2.
    elif ip == 1:
        w = -1j * be.pi/2.
    else:
        p = index - size if index > size/2 else index
        if p % 2 == 0:
            w = 2./(1.-p*p)
        else:
            w = 0.
    return w

def compute_Imo_by_extending_f_to_torus(f, W, s, Ntheta, \
                                        Sampling_on_torus):
    F = be.zeros(Sampling_on_torus, dtype = be.complex128)
    for itheta in range(Ntheta):
        F[itheta] = W[itheta] * f[itheta] * (2 * be.pi / Sampling_on_torus)
        if 0 < itheta < Ntheta:
            F[Sampling_on_torus - itheta] = (-1.)**s * \
                                            W[Sampling_on_torus - itheta] * \
                                            f[itheta] * \
                                            (2 * be.pi / Sampling_on_torus)
    Imo = fft.fft(F, norm='backward')
    return Imo

def compute_Jmo(f, Jmo, W, s, lmax, Ntheta, Sampling_on_torus):
    Imo = compute_Imo_by_extending_f_to_torus(f, W, s, Ntheta, Sampling_on_torus)
    for m in range(lmax+1):
        if m == 0:
            Jmo[m] = Imo[m]
        else:
            Jmo[m] = Imo[m] + (-1.)**s * Imo[Sampling_on_torus - m]

def axial_spin_forward_transform(f, al, s, lmax, Ntheta, D, W):
    Sampling_on_torus = 2 * (Ntheta - 1)
    Jmo = be.zeros(lmax + 1, dtype = be.complex128)
    compute_Jmo(f, Jmo, W, s, lmax, Ntheta, Sampling_on_torus)
    for l in range(abs(s), lmax + 1):
        for m in range(l % 2, l + 1, 2):
            if s > 0:
                Dlms = (-1) ** (l + m) * D[l][m][abs(s)]
            else:
                Dlms = D[l][m][abs(s)]
            al[l] += be.power(-1.j, s) * be.sqrt((2 * l + 1) / (4 * be.pi)) * \
                        D[l][m][0] * Dlms * Jmo[m]
    del Jmo