#!/usr/bin/env python
# encoding: utf-8
"""
This module contains the forward transform. The input of the function is an 
object of the class "function_in_s2". From the attributes of this object the 
band limit and the grid dimension are calculated in order to pass them on as 
arguments. The output of the function is a object of the class "salm". 
"""

from builtins import range
from coffee.settings import be

# Library from this module
from coffee.swsh.AxiSymmetric import function_class
from coffee.swsh.AxiSymmetric import axial_spin_backward_transform
from coffee.swsh.AxiSymmetric import axial_spin_forward_transform

#####################################################
# Defining the forward transform
#####################################################

def forward(f, lmax, Ntheta, deltas, weights):

    als = be.zeros((f.spin.size, lmax+1), dtype = be.complex128)
    for i in range(len(f.spin)):
        axial_spin_forward_transform.axial_spin_forward_transform(\
            f[i].map, als[i], f.spin[i], lmax, Ntheta, \
            deltas, weights) 
    return function_class.salm(als, f.spin, lmax, Ntheta, f.shape, \
                               f.spin.shape)

#####################################################
# Defining the backward transform
#####################################################

def backward(al, lmax, Ntheta, deltas):  
    
    spins = al.spins
    if spins.size == 1:
        f = be.zeros(Ntheta, dtype = be.complex128)
        axial_spin_backward_transform.axial_spin_backward_transform(\
             f, al, spins[0], lmax, Ntheta, deltas)      
        return function_class.function_in_s2(f, spins)
    else:
        f = be.zeros((spins.size, Ntheta), dtype = be.complex128)  
        for i in range(spins.size):
            axial_spin_backward_transform.axial_spin_backward_transform(\
            f[i], al[i,:], spins[i], lmax, Ntheta, \
            deltas)
        return function_class.function_in_s2(f, spins) 