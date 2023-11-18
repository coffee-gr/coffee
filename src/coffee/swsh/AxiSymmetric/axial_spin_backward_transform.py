#!/usr/bin/env python
# encoding: utf-8
"""
This module calculates the axi-symmetric backward transform.
"""

from coffee.settings import be
import scipy.fft as fft

def compute_Gmo(al, Gmo, D, s, lmax, Sampling_on_torus):
    sizeGmo_naive = 2*lmax + 1
    Gmo_naive = be.zeros(sizeGmo_naive, dtype = be.complex128)
    for m in range(lmax+1):
        for l in range(m % 2, lmax+1, 2):
            if s > 0:
                Dlms = (-1)**(l+m) * D[l][m][abs(s)]
            else:
                Dlms = D[l][m][abs(s)]
            Gmo_naive[m] += (1j)**s * be.sqrt((2.*l+1.) / (4.*be.pi)) * \
                                D[l][m][0] * Dlms * al[l]
    for m in range(1, lmax+1):
        if m != 0:
            Gmo_naive[sizeGmo_naive - m] = (-1.)**s * Gmo_naive[m]
    Gmo[0] = Gmo_naive[0]
    Gmo[1:lmax+1] = Gmo_naive[1:lmax+1]
    Gmo[Sampling_on_torus-lmax:] = Gmo_naive[sizeGmo_naive-lmax:]

def axial_spin_backward_transform(f, al, s, lmax, Ntheta, D):
    Sampling_on_torus = 2*(Ntheta-1)
    Gmo = be.zeros(Sampling_on_torus, dtype = be.complex128)
    compute_Gmo(al, Gmo, D, s, lmax, Sampling_on_torus)
    f_torus = fft.ifft(Gmo) * Sampling_on_torus
    f[:Ntheta] = f_torus[:Ntheta]