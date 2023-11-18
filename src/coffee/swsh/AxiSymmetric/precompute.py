#!/usr/bin/env python
# encoding: utf-8
"""
Here we pre-compute the Delta term and the quadratures weights.
"""

from builtins import range
from coffee.settings import be
import coffee.swsh.AxiSymmetric.delta_matrices as dm
import coffee.swsh.AxiSymmetric.weights as weights

###################################################
# Function to precompute the Delta terms
###################################################

def delta_matrices(smax, lmax):    

    Pointer_to_Deltas = dm.precompute_delta_matrices(smax + 2, lmax) 
    return Pointer_to_Deltas

###################################################
# Function to precompute the Weights
###################################################

def quadrature_weights(Ntheta):

    Pointer_to_weights = weights.precompute_quadrature_weights(Ntheta)   
    return Pointer_to_weights

def delete_precomputed_matrices(Deltas, s, lmax, Weights):

    dm.delete_delta_matrices_and_weights(Deltas, s, lmax, Weights)  

###################################################
# To create the mesh according to the transform
###################################################

def create_mesh(Ntheta):

    dtheta =  be.pi / (Ntheta-1)
    mesh = be.arange(0.0, be.pi + dtheta, dtheta)
    if mesh.shape != Ntheta:
        mesh_fixed = be.zeros(Ntheta)
        for i in range(Ntheta):
            mesh_fixed[i] = mesh[i]
        return mesh_fixed
    else:    
        return mesh






