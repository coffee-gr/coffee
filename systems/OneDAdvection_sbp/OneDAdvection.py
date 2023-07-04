#!/usr/bin/env python
# encoding: utf-8


"""
OneDAdvection.py

Created by Chris Stevens 2023
"""

# Import Python libraries
from coffee.settings import be
from coffee.tslices import tslices
from coffee.system import System

# A simple bump function for the boundary condition
def bump(t):
	return be.sin(8.*t)**7.

# Class that describes a simple advection equation
class OneDAdvection(System):

	############################################################################
	# Constructor
	############################################################################

	def __init__(self, D, tau, CFL, char_speed):
		super(OneDAdvection, self).__init__()
		self.D          = D
		self.tau        = tau
		self.CFL        = CFL
		self.char_speed = char_speed
		self.name       = "1D Advection equation"
		self.numvar     = 1

	############################################################################
	# Boundary functions
	############################################################################
		
	# Left boundary condition
	def left(self, t):
		return bump(t)
	
	############################################################################
	# Evolution Routine
	############################################################################

	def evaluate(self, t, U, intStep = None):

		# Define useful variables
		Psi = U.data[0]
		dx  = U.domain.step_sizes[0]

		# Calculate spatial derivative
		DxPsi = self.D(Psi, dx)

		# Evaluate right hand side
		DtPsi = -self.char_speed*DxPsi

		# Impose boundary condition with SAT method
		pt_l     = self.D.penalty_boundary(dx, "left")
		C_left   = self.tau * pt_l * (-self.char_speed)
		l_Psi    = self.left(t)
		DtPsi[0] -= C_left * (Psi[0] - l_Psi)

		# Output time if desired
		# print('t =', t)

		# Return RHS of evolution equation with SAT boundary condition imposed
		return tslices.TimeSlice([DtPsi], U.domain, time = t)

	############################################################################
	# Initial Data Routine
	############################################################################
	
	def initial_data(self, t0, grid):

		# Set initial data
		x   = grid.meshes[0]
		Psi = be.zeros_like(x)

		# Return Timeslice object
		return tslices.TimeSlice([Psi], grid, time = t0)

	############################################################################
	# Timestep Routine
	############################################################################

	def timestep(self, U):
		return self.CFL * U.domain.step_sizes[0] / self.char_speed