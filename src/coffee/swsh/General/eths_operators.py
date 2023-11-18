#!/usr/bin/env python
# encoding: utf-8
"""
This module contains the eth-operators. "ethU" is the name for the eth-operator 
that raises the spin and "ethD" for the one who lowers it. The input of these 
functions is an object of the class "function_in_s2". From the attributes of 
this object the band limit and the grid dimension are calculated in order to 
pass them on as arguments. The output of the functions is a object of the class 
"function_in_s2" with the spin raised or lowered depending the spin weight of 
the original function.
"""


# standard imports
from coffee.settings import be

# Library from this module
from coffee.swsh.General import function_class
from coffee.swsh.General import eths

#####################################################
# Defining eth that raise spin
#####################################################

def ethU(f, lmax, smax, Ntheta, Nphi, deltas, weights): 

    spins = be.copy(f.spin)
    f     = be.copy(f.map)

    spins_flat = spins.reshape(spins.size)  
    f_flat = f.reshape(spins.size, Ntheta*Nphi)   
    number_of_alm_per_lmax = lmax*(lmax+2) + 1   
   
    for i in range( len(spins_flat) ): 
        alm = be.zeros(number_of_alm_per_lmax, dtype = be.complex128) 
        eths.spin_ethU(f_flat[i], alm, spins_flat[i], lmax, smax, Ntheta, \
                       Nphi, deltas, weights)

    fnew = f_flat.reshape(f.shape)

    return function_class.function_in_s2(fnew, spins + be.ones_like(spins))

#####################################################
# Defining eths that lower spin
#####################################################

def ethD(f_, lmax, smax, Ntheta, Nphi, deltas, weights): 
   
    spins = be.copy(f_.spin)  
    f     = be.copy(f_.map)

    #flatting the tensors             
    spins_flat = spins.reshape(spins.size)  
    f_flat     = f.reshape(spins.size, Ntheta * Nphi)    
    number_of_alm_per_lmax = lmax*(lmax+2) + 1   
   
    for i in range( len(spins_flat) ):         
        alm = be.zeros(number_of_alm_per_lmax, dtype = be.complex128) 
        eths.spin_ethD(f_flat[i], alm, spins_flat[i], lmax, smax, Ntheta, \
                        Nphi, deltas, weights)
   
    fnew = f_flat.reshape(f.shape)

    return function_class.function_in_s2(fnew, spins - be.ones_like(spins))


    