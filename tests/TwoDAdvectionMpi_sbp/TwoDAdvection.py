#!/usr/bin/env python
# encoding: utf-8


"""
TwoDAdvection.py

Created by Chris Stevens 2023
"""

# Import Python libraries
from functools import partial

from coffee.settings import be
from coffee.tslices import tslices
from coffee.system import System
from coffee.diffop.fd import ghost_point_processor

# A simple bump function for the boundary condition
def bump(t, grid, grid_global):
	start = grid_global[0]
	stop  = grid_global[-1]
	return be.sin(be.pi*start / (start - stop) - \
	       be.pi / (start - stop)*grid)**8.* \
		   be.sin(8.*t)**7.

# Class that describes a simple advection equation
class TwoDAdvection(System):

	################################
	# Constructor
	################################

	def __init__(self, Dx, Dy, tau, CFL, char_speed_x, char_speed_y, \
	      		x_global, y_global):
		super(TwoDAdvection, self).__init__()
		self.Dx           = Dx
		self.Dy           = Dy
		self.tau          = tau
		self.CFL          = CFL
		self.char_speed_x = char_speed_x
		self.char_speed_y = char_speed_y
		self.x_global     = x_global
		self.y_global     = y_global
		self.name         = "2D Advection equation"
		self.numvar       = 1

	############################################################################
	# Boundary functions
	############################################################################
		
	# Left x boundary condition
	def left_x(self, t, y, y_global):
		return bump(t, y, y_global)
	
	# Left y boundary condition
	def left_y(self, t, x, x_global):
		return bump(t, x, x_global)

	############################################################################
	# Evolution Routine
	############################################################################

	def evaluate(self, t, U, intStep = None):

		# Define useful variables
		Psi = U.data[0]
		x   = U.domain.axes[0]
		y   = U.domain.axes[1]
		dx  = U.domain.step_sizes[0]
		dy  = U.domain.step_sizes[1]

		# Calculate spatial derivatives by applying the 1D FD operators
		# over the mesh
		DxPsi = be.apply_along_axis(
			lambda f:self.Dx(f, dx),
			0, 
			Psi
			)
		DyPsi = be.apply_along_axis(
			lambda f:self.Dy(f, dy),
			1, 
			Psi
			)

		# Communicate spatial derivative through neighbouring MPI process(es)
		new_derivatives, _ = U.communicate(
			partial(ghost_point_processor),
			data=be.array([
				DxPsi, DyPsi
			])
		)
		DxPsi, DyPsi = new_derivatives

		# Impose boundary conditions with SAT method
		pt_x_l = self.Dx.penalty_boundary(dx, "left")
		pt_y_l = self.Dy.penalty_boundary(dy, "left")
	
		# Only two of the four boundaries require boundary conditions
		C_x_left     = self.tau * pt_x_l * (-self.char_speed_x)
		l_x_Psi      = self.left_x(t, y, self.y_global)
		C_y_left     = self.tau * pt_y_l * (-self.char_speed_y)
		l_y_Psi      = self.left_y(t, x, self.x_global)

		# Evaluate right hand side
		DtPsi = -self.char_speed_x*DxPsi - self.char_speed_y*DyPsi

		# Do the external boundaries
		b_data = U.external_slices()
		for dim, direction, _ in b_data:
			if dim == 0 and direction == -1: # Left x boundary)
				DtPsi[0,:] -= C_x_left * \
					(Psi[0,:] - l_x_Psi)
			elif dim == 1 and direction == -1: # Left y boundary
				DtPsi[:,0] -= C_y_left * \
					(Psi[:,0] - l_y_Psi)

		# # Output time if desired
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

		return self.CFL * min(U.domain.step_sizes[0] / self.char_speed_x, \
	  						  U.domain.step_sizes[1] / self.char_speed_y)