#!/usr/bin/env python
# encoding: utf-8


"""
OneDAdvection.py

Created by Chris Stevens 2023
"""

# Import Python libraries
from coffee.settings import be
from coffee.system import System
from coffee.tslices import tslices
from coffee.diffop.fd import ghost_point_processor
from functools import partial

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
	# Boundary functions (called from ibvp.)
	############################################################################

	# Impose Dirichlet boundary data on the left boundary
	def dirichlet_boundary(self, U):

		# Left boundary process
		if U.domain.mpi.comm.rank == 0:
			U.data[0][0] = self.left(U.time)
		return U
		
	# The specific boundary condition for the left boundary
	def left(self, t):
		return bump(t)
	
	############################################################################
	# Evolution Routine
	############################################################################

	def evaluate(self, t, U):

		# Define useful variables
		Psi = U.data[0]
		dx  = U.domain.step_sizes[0]

		# Calculate spatial derivative
		DxPsi = self.D(Psi, dx)

		# Communicate spatial derivative through neighbouring MPI process(es)
		new_derivatives, _ = U.communicate(
			partial(ghost_point_processor),
			data=be.array([
				DxPsi
			])
		)
		DxPsi = new_derivatives[0]

		# Evaluate right hand side
		DtPsi = -self.char_speed*DxPsi

		# Output time if desired
		# print('t =', t)

		# Return RHS of evolution equation
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