#!/usr/bin/env python
# encoding: utf-8
"""
This module contains operators to work numerically with functions over S2.
The code is based on 
"Fast and exact spin-s Spherical Harmonic Transformations" in Astrophysical 
Journal Supplement Series (2010) 189:255-260."
"""

from coffee.swsh.General import transforms
from coffee.swsh.General import eths_operators 
from coffee.swsh.General import precompute 
from coffee.swsh.General import function_class 
from coffee.swsh.General import generator_of_swsh
from coffee.settings import be
import sys


class general_swsh:
	
	def __init__(self, Ntheta, Nphi, smax ):			 

		# Define the band limit
		self.Ntheta = Ntheta
		self.Nphi   = Nphi 
		self.smax   = smax
		self.lmax = int((2./3.)*self.Ntheta)

		# Initialise the Wigner matrices
		self.deltas  = precompute.delta_delta_matrices(self.lmax , self.smax)
		self.weights = precompute.quadrature_weights(self.Ntheta)	

		# Initialise mesh
		self.PHI, self.THETA = precompute.create_mesh(self.Ntheta, self.Nphi)
	
	################################### 
	def function_in_s2(self, f, spins):
		"""Returns an instance of the function_in_s2 class

        Parameters
        ----------

        f : be.array
            An array of 2D arrays representing the discretised functions.

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
	def salm(self, al_flat, spins_flat, lmax, Ntheta, nphi, f_shape, spins_shape):
		return function_class.salm(al_flat, spins_flat, lmax, Ntheta, nphi, \
			     					f_shape, spins_shape)

	################################### 
	def coefficient_aslm(self, s, l, m, values_of_coefficients):

		spins        = be.asarray(s)
		spins_shape  = spins.shape
		coefficients = be.asarray(values_of_coefficients)
		
		L  = be.asarray(l)
		M  = be.asarray(m)

		if (spins_shape != L.shape) or (spins_shape != M.shape) \
									or (spins_shape != coefficients.shape):  	
			print("Error in the input of the aslm")
			sys.exit()

		spins_shape = spins.shape
		f_shape     = self.THETA.shape 

		number_of_alm_per_lmax = self.lmax*(self.lmax+2) + 1 
			
		if spins.shape == (0,) or spins.shape == ():
			spins_flat = spins.reshape(spins.size) 
			alm_flat   = be.zeros((1, number_of_alm_per_lmax), \
									dtype = be.complex128) 
			alm_flat[0, l*(l+1) + m] = coefficients 

		else:
			spins_flat = spins.reshape(spins.size) 
			alm_flat   = be.zeros((spins.size, number_of_alm_per_lmax), \
									dtype = be.complex128) 
			for i in range(spins.size):	
				alm_flat[spins_flat[i], L[i]*(L[i]+1) + M[i]] = coefficients[i] 
		
		return function_class.salm(alm_flat, spins_flat, self.lmax, \
								self.Ntheta, self.Nphi, f_shape, spins_shape)

	###################################  
	def forward(self, f):	
		return transforms.forward(f, self.lmax, self.smax, self.Ntheta, \
						self.Nphi, self.deltas, self.weights) 

	###################################  
	def backward(self, al):
		return transforms.backward(al, self.lmax, self.smax, self.Ntheta, \
						self.Nphi, self.deltas) 

	###################################  
	def ethU(self, f):
		return eths_operators.ethU(f, self.lmax, self.smax, self.Ntheta, \
						self.Nphi, self.deltas, self.weights) 	

	###################################  
	def ethD(self, f):	
		return eths_operators.ethD(f, self.lmax, self.smax, self.Ntheta, \
						self.Nphi, self.deltas, self.weights) 

	###################################  
	def C(self, f):   
		return function_class.function_in_s2(f.map.conjugate(), (-1) * f.spin)

	###################################  
	def Yslm(self, s, l, m):
		return generator_of_swsh.Yslm(self.THETA, self.PHI, s, l, m, \
										self.lmax, self.smax)