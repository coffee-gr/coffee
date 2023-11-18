#!/usr/bin/env python
# encoding: utf-8
"""
This module contains the forward transform. The input of the function is an 
object of the class "function_in_s2". From the attributes of this object the 
band limit and the grid dimension are calculated in order to pass them on as 
arguments. The output of the function is a object of the class "salm". 
"""
  
# standard imports
from coffee.settings import be

# library from this module
from coffee.swsh.General import function_class
from coffee.swsh.General import spin_backward_transform
from coffee.swsh.General import spin_forward_transform

#####################################################
# Defining the forward transform
#####################################################

def forward(f, lmax, smax, Ntheta, Nphi, deltas, weights): 

    f_flat     = f.map.reshape(f.spin.size, Ntheta*Nphi) 
    number_of_alm_per_lmax = lmax*(lmax+2) + 1   
    alm_flat   = be.zeros((f.spin.size, number_of_alm_per_lmax), \
                          dtype = be.complex128)

    for i in range(f.shape[0]):
        spin_forward_transform.spin_forward_transform(f_flat[i], alm_flat[i], f.spin[i], \
                                  lmax, smax, Ntheta, Nphi, deltas, weights)
        
    return function_class.salm(alm_flat, f.spin, lmax, Ntheta, Nphi, \
                               f.shape, f.spin.shape)
 
#####################################################
# Defining the backward transform
#####################################################

def backward(alm, lmax, smax, Ntheta,  Nphi, deltas):

    spins = alm.spins   
    if len(spins) == 1:
        f = be.zeros((Ntheta, Nphi), dtype = be.complex128)
        alm_flat = alm.vector.flatten()
       
        spin_backward_transform.spin_backward_transform(f, alm_flat, spins[0], lmax, smax, Ntheta, Nphi, \
                                    deltas)    

        return function_class.function_in_s2(f, spins[0])

    else:             
        spins_flat = spins.reshape(spins.size)  
        f_flat     = be.zeros((spins.size, Ntheta, Nphi), \
                                dtype = be.complex128)
        alm_flat   = filling_alm_flat(alm, lmax)

        for i in range(len(spins_flat)):     
            spin_backward_transform.spin_backward_transform(f_flat[i], alm_flat[i], spins_flat[i], \
                                    lmax, smax, Ntheta, Nphi, deltas)    

        return function_class.function_in_s2(f_flat, spins)


def filling_alm_flat(alm, lmax ):
    spins  = alm.spins
    number_of_alm_per_lmax = lmax*(lmax+2) + 1 
    alm_flat = be.empty((spins.size, number_of_alm_per_lmax), \
                        dtype = be.complex128) 
   
    for s in range(spins.size):
        for l in range(alm.lmax + 1): 
            m = -l
            while (m <= l):    
                alm_flat[s, l*(l+1) + m] = alm[s, l, m] 
                m += 1
    return alm_flat
