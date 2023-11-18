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

from coffee.settings import be

# Library from this module
from coffee.swsh.AxiSymmetric import function_class
from coffee.swsh.AxiSymmetric import axial_spin_backward_transform
from coffee.swsh.AxiSymmetric import axial_spin_forward_transform

#####################################################
# Defining eth that raise spin
#####################################################

def ethU(fns, lmax, Ntheta, deltas, weights):

    Deth_fn = be.zeros(fns.shape, dtype = be.complex128)
    for i in range(0, len(fns)):
        al = be.zeros(lmax + 1, dtype = be.complex128)
        axial_spin_forward_transform.axial_spin_forward_transform(\
                                    fns[i].map, al, fns.spin[i], lmax, Ntheta, \
                                    deltas, weights)
        for l in range(lmax + 1):
            if l < abs(fns.spin[i]):
                al[l] = 0.0
            else:
                al[l] = -be.sqrt((l - fns.spin[i]) * (l + fns.spin[i] + 1.)) * al[l]
        axial_spin_backward_transform.axial_spin_backward_transform(\
                                        Deth_fn[i], al, fns.spin[i] + 1, lmax, \
                                        Ntheta, deltas)
    
    return function_class.function_in_s2(Deth_fn, fns.spin+1) 

#####################################################
# Defining eth that lower spin
#####################################################

def ethD(fns, lmax, Ntheta, deltas, weights):  

    Dethp_fn = be.zeros(fns.shape, dtype = be.complex128)
    for i in range(0, len(fns)):
        al = be.zeros(lmax + 1, dtype = be.complex128)
        axial_spin_forward_transform.axial_spin_forward_transform(\
                                    fns[i].map, al, fns.spin[i], lmax, Ntheta, \
                                    deltas, weights)
        for l in range(lmax + 1):
            if l < abs(fns.spin[i]):
                al[l] = 0.0
            else:
                al[l] = be.sqrt((l + fns.spin[i]) * (l - fns.spin[i] + 1.)) * al[l]
        axial_spin_backward_transform.axial_spin_backward_transform(\
                                        Dethp_fn[i], al, fns.spin[i] - 1, \
                                        lmax, Ntheta, deltas)

    return function_class.function_in_s2(Dethp_fn, fns.spin-1)