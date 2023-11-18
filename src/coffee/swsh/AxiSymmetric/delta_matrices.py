#!/usr/bin/env python
# encoding: utf-8
"""
This module calculates the Wigner-D matrices.
"""

from coffee.settings import be

def initialize_delta_matrices(s, lmax):
    D = be.empty((lmax+1, lmax+1, s+1))
    for l in range(lmax+1):
        D[l] = be.zeros((lmax+1, s+1))
    return D

def delete_delta_matrices(D):
    del D

def compute_delta_matrices(D, s, lmax):
    D[0, 0, 0] = 1.0
    for l in range(1, lmax+1):
        for m2 in range(min(l, s)+1):
            if m2 == 0:
                D[l, l, 0] = -be.sqrt((2.*l-1)/(2.*l)) * D[l-1, l-1, 0]
            else:
                D[l, l, m2] = be.sqrt((l/2.)*(2.*l-1) / ((l+m2)*(l+m2-1.))) * \
                                 D[l-1, l-1, m2-1]
            for m1 in range(l-1, -1, -1):
                if m1 == l-1:
                    D[l, l-1, m2] = (2.*m2)/be.sqrt((l-m1)*(l+m1+1.)) * \
                                        D[l, l, m2]
                else:
                    D[l, m1, m2] = (2.*m2)/be.sqrt((l-m1)*(l+m1+1.)) * \
                                    D[l, m1+1, m2] - \
                                    be.sqrt(((l-m1-1.)*\
                                    (l+m1+2.))/((l-m1)*(l+m1+1.))) * \
                                    D[l, m1+2, m2]

def show_delta_matrices(D, s, lmax):
    print("\nChecking components of D(l,m1,m2)\n")
    for i in range(lmax+1):
        for j in range(i+1):
            for k in range(min(i, s)+1):
                print("D({}, {}, {}) = {}".format(i, j, k, D[i, j, k]))

def precompute_delta_matrices(s, lmax):
    D = initialize_delta_matrices(abs(s) + 1, lmax)
    compute_delta_matrices(D, abs(s) + 1, lmax)
    return D

def delete_delta_matrices_and_weights(D, s, lmax, W):
    del W
    delete_delta_matrices(D, s, lmax)