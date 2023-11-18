#!/usr/bin/env python
# encoding: utf-8
"""
Here we pre-compute the Delta term and the quadratures weights.

For best accuracy

a. The FFT transform in theta fastest for Ntheta = 2^n + 1.
b. The FFT transform in theta fastest for Nphi   = 2^n.
"""

# Standard imports
from coffee.settings import be
import coffee.swsh.General.delta_matrices as dm
import coffee.swsh.General.weights as weights

################################################################################
# Function to precompute the Delta terms
################################################################################

def delta_delta_matrices(lmax, smax):    

    number_of_Dlmn = 0
    for l in range(lmax + 1):                
        number_of_Dlmn += ((2*l+1))**(2.)    

    delta = dm.Precompute_Delta_Delta_Matrices(lmax, smax)

    return delta

################################################################################
# Function to delete the delta matrices and weights
################################################################################

def delete_delta_delta_matrices_and_weights(DD, lmax, smax, W): 
    dm.Delete_Delta_Delta_Matrices(DD, lmax, smax, W) 

################################################################################
# Function to precompute the weights
################################################################################

def quadrature_weights(Ntheta):

    pointer_to_weights = weights.Precompute_quadrature_weights(Ntheta)  

    return pointer_to_weights

################################################################################
# Create the mesh according to the transform
################################################################################

def create_mesh(Ntheta, Nphi):

    dtheta = be.pi/(Ntheta-1) 
    dphi   = 2.*be.pi/(Nphi)
    NTHETA = be.arange(0.0, be.pi + dtheta, dtheta)
    NPHI   = be.arange(0.0, 2*be.pi, dphi)

    PHI, THETA = be.meshgrid(NPHI, NTHETA)

    return  PHI, THETA








