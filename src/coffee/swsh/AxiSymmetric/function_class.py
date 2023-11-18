#!/usr/bin/env python
# encoding: utf-8
"""
This module contains the class of "function_in_s2" and "salm". The first one 
(the most important class of this module) defines functions in S2 with 
spin weight and the rules to conduct basic operations

(addition, substraction, multiplication, division and potentiation). 

To build a object of this class two arguments have to be provided:

    1. The definition of the function (or arrays of functions) with dimension of 
        Ntheta. i.e., the function in the collocation space. 

    2. The spin weight (or arrays of spins weights) associated to the function
        (or list of functions). 

This class then defines a function, say f, with two attributes: 

f.map (functions over the grid), f.spin(spin weight)

For the general case when you are dealing with arrays of functions each 
component has these attributes.

The second class defines spectral coefficients and the way to print them. 
The objects of this class are built by applying the "forward transform" to
an object of the previous class. i.e., a well defined function with spin weight. 
Each object of the class salm contains information (the attributes) from the 
function they come from that is used for the backward transform. The spectral 
coefficients are shown as follows: 

salm is a collection of sets of spectral coefficients. Each set correspond in 
the same order to one component of the list of spins that results from 
flatting the spin array. i.e., if the spin array is of the form

[[s1, s2], [s3, s4]], 

the list of spins has 4 elements listed by [s1, s2, s3 ,s4]. You can get these 
elements of the object salm in several ways:

            1. salm[0]    : returns from the first set of elements salm.
            
            2. salm[0,1]  : returns from the first set salm, the coefficients 
                            with the first spin given by l=1

Finally, we have to point out that unlike the first class (function_in_s2), 
this class DOES NOT defines operations among its objects. So, standard rules 
among arrays apply here. 
"""

from builtins import str
from builtins import range
from coffee.settings import be

###############################################
# This is the class for spin-weighted functions
###############################################

class function_in_s2(be.ndarray): 

    # Defining a spin weight function
    def __new__(self, array, spins):        
        obj      = be.array(array).view(self)
        obj.map  = be.array(array)       
        obj.spin = be.array(spins)  
        return obj

    def __str__(self):
        return str(self.map)    

    def __getitem__(self, key):  
        key = be.index_exp[key] 
        if  len( self.spin) < len(key):
           raise IndexError("Index of the array functions is not valid!") 
        return function_in_s2(self.map[key], self.spin[key])

    ######## Defining the operations between spin weight functions #########  
    
    def __add__(self, other):   
        if isinstance(other, (int, float, complex)):
                if self.spin == 0:
                    return function_in_s2(self.map + other, self.spin)  
                else:
                    raise IndexError("Spin weight of the functions do not match!")   

        elif self.spin.all() == other.spin.all():
                return function_in_s2(self.map + other.map, self.spin) 
        else: 
            raise IndexError("Spin weight of the functions do not match!")
  
    def __radd__(other, self):
        if isinstance(self, (int, float, complex)):
                if other.spin == 0 :
                    return function_in_s2( other.map + self , other.spin ) 
                else:
                    raise IndexError("Spin weight of the functions do not match!")
       
        elif self.spin.all() == other.spin.all():
                return function_in_s2(self.map + other.map, self.spin) 
        else: 
            raise IndexError("Spin weight of the functions do not match!")   

    def __sub__(self, other):
        if isinstance(other, (int, float, complex)):
                if self.spin == 0:
                    return function_in_s2(self.map - other, self.spin)  
                else:
                    raise IndexError("Spin weight of the functions do not match!")                           
        elif self.spin.all() == other.spin.all():
                return function_in_s2(self.map - other.map, self.spin) 
        else: 
            raise IndexError("Spin weight of the functions do not match!") 

    def __rsub__(other, self):
        if isinstance(self, (int, float, complex)):
                if other.spin == 0:
                    return function_in_s2( - other.map + self, other.spin)  
                else:
                    raise IndexError("Spin weight of the functions do not match!") 

    def __mul__(other, self): 
            if isinstance(other, (int, float, complex)):
                return function_in_s2(self.map * other, self.spin) 
            elif isinstance(self, (int, float, complex)):
                return function_in_s2(self * other.map, other.spin) 
            else:
                return function_in_s2(self.map * other.map, self.spin + \
                                                            other.spin)   

    def __rmul__(self, other):
            return function_in_s2(self.map * other, self.spin)

    def __div__(self, other):
            if isinstance(other, (int, float, complex)):
                return function_in_s2(self.map / other, self.spin) 
            elif isinstance(self, (int, float, complex)):
                return function_in_s2(self /  other.map, other.spin) 
            else:
                return function_in_s2(self.map / other.map, self.spin - \
                                                            other.spin)                    
            
    def __pow__(self, n ):
        return function_in_s2(self.map**n, int(be.rint(self.spin * n))) 

    def __neg__(self):
        return function_in_s2(-self.map, self.spin)

    def info(self):
        return ['map','spin']

###########################################
# This is the class for spectral elements
###########################################

class salm(be.ndarray):
    
    def __new__(self, array, spins, lmax, Ntheta, f_shape, spin_shape):             
        obj = be.asarray(array).view(self)      
        obj.spins = be.asarray(spins)
        obj.spin_shape = spin_shape
        obj.lmax = lmax      
        obj.Ntheta = Ntheta
        obj.f_shape = f_shape        
        obj.vector = array

        return obj  
   
    def __str__(self):
        """
        Return: Returns all info of the object sal
        """ 

        s = "------------------------------------------------------------------\n"             
        s += "Total number of sets of a(l) : %s \n"% self.spins.size
        s += "-----------------------------------------------------------------\n" 
        if self.spins.size == 1: 
            s += " s = %d  \n"%self.spins                      
            for l in range(self.lmax + 1):
               s += "a(" + str(l) + ") = " + str(self.view(be.ndarray)[0][l]) + "\n"                      
            s += "------------------------------------------------------------------\n"                   
       
        else:
            for j in range(self.spins.size):                 
                s += "(set #: %s), spin = %d \n"%(str(j), self.spins[j])   
                for l in range(self.lmax + 1):
                    s += "a(" + str(l) + ") = " + str(self.view(be.ndarray)[j][l]) + "\n"         
                s += "------------------------------------------------------------------\n"   
        return s 
    
    def __getitem__(self, key):      
        key = be.index_exp[key]

        if self.spins.size == 1: 
            if len(key) == 1 :                                                          
                    return be.around(self.view(be.ndarray)[0, key[0]], 15) 

            elif len(key) >= 2: 
                    raise IndexError("\nToo many indices! Object'ss elements are given in this case as al[l]\n") 

        else:  
            if len(key) == 1:  
                s = "(set #: %s), spin = %d \n" % (key[0], self.spins[key[0]])                
                for l in range(self.lmax + 1):
                     s += "a(%i) = %e \n" % \
                        (l, be.around(self.view(be.ndarray)[key[0], l], 15))
                return s                 

            if len(key) == 2:                    
                    return be.around(self.view(be.ndarray)[key[0], key[1]], 15)

            elif len(key) >= 3: 
                    raise IndexError("\nToo many indices! Object's elements are given in this case as al[ #set ,l]\n") 


