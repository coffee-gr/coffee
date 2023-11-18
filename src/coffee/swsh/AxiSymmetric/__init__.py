#!/usr/bin/env python
# encoding: utf-8
"""
This module contains operators to work numerically with functions over S2 
with axi-symmetry. Even thought the code is based on 
"Fast and exact spin-s Spherical Harmonic Transformations" in Astrophysical 
Journal Supplement Series (2010) 189:255-260." it is structurally different due 
to simplifications arising in axi-symmetry.

The first simplification is respect to the 2-dimensional FFT done over the 
extended function in the torus. Here, due to the axi-symmetry, we use a 
1-dimension FFT. the second simplification lies in the computation of the 
Wigner matrices.
"""

from builtins import object
from coffee.swsh.AxiSymmetric import transforms
from coffee.swsh.AxiSymmetric import eths_operators 
from coffee.swsh.AxiSymmetric import precompute 
from coffee.swsh.AxiSymmetric import function_class 

class axi_symmetric_swsh(object):
	def __init__(self, Ntheta, smax):	
		
		# Define the band limit
		self.Ntheta = Ntheta
		self.lmax = int((2./3.)*self.Ntheta)

		# Initialise the Wigner matrices
		self.deltas  = precompute.delta_matrices(smax, self.lmax)
		self.weights = precompute.quadrature_weights(self.Ntheta)

	################################### 
	def function_in_s2(self, f, spins):
		"""Returns an instance of the function_in_s2 class

        Parameters
        ----------

        f : be.array
            An array of 1D arrays representing the discretised functions.

        spins : be.array
            An array containing the spins of the functions.

        Returns
        -------

        function_class.function_in_s2
            The overloaded be.array that represents an array of descretised 
	    	functions.

        """
		return function_class.function_in_s2(f, spins)

	###################################
	def salm(self, al_flat, spins_flat, lmax, Ntheta, f_shape , spins_shape):
		return function_class.salm(al_flat, spins_flat, lmax, Ntheta, \
			     					f_shape, spins_shape)

	################################### 
	def create_mesh(self, Ntheta):
		return precompute.create_mesh(Ntheta)

	###################################  
	def forward(self, f):	
		return transforms.forward(f, self.lmax, self.Ntheta, \
			    					self.deltas, self.weights) 

	###################################  
	def backward(self, al):
		return transforms.backward(al, self.lmax, self.Ntheta, \
			     					self.deltas) 

	###################################  
	def ethU(self, f):
		return eths_operators.ethU(f, self.lmax, self.Ntheta, \
			     					self.deltas, self.weights) 		

	###################################  
	def ethD(self, f):	
		return eths_operators.ethD(f, self.lmax, self.Ntheta, \
			     					self.deltas, self.weights) 
	
	###################################  
	def C(self, f):  
	    return function_class.function_in_s2(f.map.conjugate(), (-1) * f.spin)


