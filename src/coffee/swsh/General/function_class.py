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

from coffee.settings import be

###############################################
# This is the class for spin weight  functions
###############################################

class function_in_s2(be.ndarray): 

    # Defining a spin weight function
    def __new__(self, array, spins ):        
        obj      = be.asarray( array ).view(self)
        obj.map  = be.asarray(array) 
        obj.spin = be.asarray(spins) 
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
    def __new__(self, array, spins, lmax, Ntheta, Nphi, f_shape , spin_shape ):  
        obj = be.asarray(array).view(self)      
        obj.spins = be.asarray(spins)
        obj.spin_shape=spin_shape
        obj.lmax = lmax
        obj.Nphi = Nphi
        obj.Ntheta = Ntheta
        obj.f_shape = f_shape   
        obj.vector= array

        return obj  
  
    def __str__(self):
        """
        Return: Returns all info of the object salm
        """ 
        s = "------------------------------------------------------------------\n"             
        s += "Total number of sets of salm : %s \n"% self.spins.size
        s += "------------------------------------------------------------------\n" 
        if self.spins.shape == (0,) or self.spins.shape == ():
            s += " s = %d  [ -m,... 0,.. m ] \n"%self.spins   
            m = 0
            i = 0         
            for l in range(self.lmax + 1):
                line = []  
                m = (2*l+1)+m 
                while (i < m):  
                    line.append(be.around(self.view(be.ndarray)[0,i],13))
                    i += 1
                s +=" l = %d: %s \n"%(l, line)                
            s += "------------------------------------------------------------------\n"                   
       
        else:
            for j in range(self.spins.size):        
                s += "(set #: %s)\n" % str(j+1)
                s += " s = %d  [ -m,... 0,... m ] \n" % self.spins[j]
                m = 0
                i = 0
                for l in range(self.lmax + 1):
                    line = []  
                    m = (2*l+1)+m 
                    while (i < m):  
                        line.append(be.around(self.view(be.ndarray)[j,i],13))
                        i += 1
                    s += " l = %d: %s \n"%(l, line)                
                s += "------------------------------------------------------------------\n"   
        return s
  
    def __getitem__(self, key):      
        key = be.index_exp[key] 

        if self.spins.shape == (0,) or self.spins.shape == ():
            if len(key) == 1 :
                    line = []  
                    i = key[0]*(key[0]+1)
                    m = -key[0]  
                    while (m < key[0]+1): 
                        line.append(be.around(self.view(be.ndarray)[0,i+m], 15))
                        m = m + 1                                            
                    return be.asarray(line)                    

            if len(key) == 2:
                    i = key[0]*(key[0]+1)  
                    return be.around(self.view(be.ndarray)[0, i+key[1]], 15)

            elif len(key) >= 3: 
                    raise IndexError("\nToo many indices! Object's elements are given in this case as  alm[l,m]\n")
        else:  
            if len(key) == 1:  
                s = " s = %d  [ -m,... 0,... m ] \n" % self.spins[key[0]]
                m = 0
                i = 0
                for l in range(self.lmax + 1):
                    line = []  
                    m = (2*l+1)+m 
                    while (i < m):  
                        line.append(be.around(self.view(be.ndarray)[key[0],\
                                                                     i], 15))
                        i += 1
                    s += " l = %d: %s \n" % (l, line) 
                return s  

            if len(key) == 2:
                line = []  
                i = key[1]*(key[1]+1)
                m = -key[1]  
                while (m < key[1] + 1):  
                    line.append(be.around(self.view(be.ndarray)[key[0],\
                                                                i+m], 15))
                    m = m + 1
                return be.asarray(line)                    

            if len(key) == 3:
                    i = key[1]*(key[1]+1)  
                    return be.around(self.view(be.ndarray)[key[0],\
                                                           i+key[2]], 15)

            elif len(key) >= 3: 
                    raise IndexError("\nToo many indices! Object's elements are given in this case as alm[ s_n ,l,m]\n")
