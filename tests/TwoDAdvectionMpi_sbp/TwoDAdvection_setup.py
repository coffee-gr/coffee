#!/usr/bin/env python
# encoding: utf-8


"""
TwoDAdvection_setup.py

Created by Chris Stevens 2023
"""

# Import python libraries
from builtins import range
import sys
import os
import numpy as np
import h5py
from mpi4py import MPI

filename     = sys.argv[1]
backend_type = sys.argv[2]

# Import standard code base and initialize NumPy backend
from coffee.settings import init
from coffee.backend import backend as be
my_backend = be.set_backend(backend_type)
init(my_backend)

from coffee import ibvp, actions, solvers, grid
from coffee.diffop.sbp import sbp

# Import system to use
import TwoDAdvection

np.set_printoptions(threshold=np.inf, precision=16)

################################################################################  
# These are the commonly altered settings
################################################################################

# How many systems?
num_of_grids = 1

# How many grid points?
Nx = 400
Ny = 400

# What grid to use?
xstart = -1
xstop  = 1
ystart = -1
ystop  = 1

# Times to run between
tstart = 0.0
tstop  = 1.

# Configuration of System
CFL = 0.5
tau = 1.0

# Differential operator for x-direction
diffop_x = sbp.D43_Strand(sbp.BOUNDARY_TYPE_GHOST_POINTS)
# diffop_x = sbp.D42(sbp.BOUNDARY_TYPE_GHOST_POINTS)
# diffop_x = sbp.D43_Tiglioetal(sbp.BOUNDARY_TYPE_GHOST_POINTS)
# diffop_x = sbp.D43_CNG(sbp.BOUNDARY_TYPE_GHOST_POINTS)

# Differential operator for y-direction
diffop_y = sbp.D43_Strand(sbp.BOUNDARY_TYPE_GHOST_POINTS)
# diffop_y = sbp.D42(sbp.BOUNDARY_TYPE_GHOST_POINTS)
# diffop_y = sbp.D43_Tiglioetal(sbp.BOUNDARY_TYPE_GHOST_POINTS)
# diffop_y = sbp.D43_CNG(sbp.BOUNDARY_TYPE_GHOST_POINTS)

################################################################################      
## MPI set up                                                                         
################################################################################ 

dims     = MPI.Compute_dims(MPI.COMM_WORLD.size, [0,0])                                   
periods  = [0,0]                                                                        
reorder  = True                                                                       
mpi_comm = MPI.COMM_WORLD.Create_cart(dims, periods=periods, reorder=reorder)        

################################################################################
# Grid construction
################################################################################

# Grid point data      
axes = [(Nx*2**i,Ny*2**i) for i in range(num_of_grids)]

# Determine the boundary data
ghost_points = (diffop_x.ghost_points(), diffop_y.ghost_points())
internal_points = (diffop_x.internal_points(), diffop_x.internal_points())
b_data = grid.MPIBoundary(
    ghost_points, 
    internal_points, 
    mpi_comm=mpi_comm, 
    number_of_dimensions=2
)

# Build grids
grids = [
    grid.UniformCart(
        axes[i], 
        [(xstart, xstop), (ystart, ystop)],
        comparison = i,
        mpi_comm = mpi_comm,
        boundary_data=b_data
    ) 
    for i in range(num_of_grids)
]

# Define global grids
x_global = np.linspace(xstart, xstop, Nx+1)
y_global = np.linspace(ystart, ystop, Ny+1)

################################################################################
# Initialise systems
################################################################################

char_speed_x = 1. # Must be positive
char_speed_y = 1. # Must be positive
systems = []
for i in range(num_of_grids):
    systems += [TwoDAdvection.TwoDAdvection(\
        diffop_x, diffop_y, tau, CFL, char_speed_x, char_speed_y, \
        x_global, y_global
        )]
    
# Configuration of IBVP
solvers = [solvers.RungeKutta4(sys) for sys in systems]
maxIteration = 10000000

################################################################################
# Set up hdf file to store output
################################################################################

hdf_file = h5py.File(filename, "w")

################################################################################
# Set up action types for data storage in hdf file
################################################################################

output_actions = [
    actions.SimOutput.Data(),
    actions.SimOutput.Times(),
    actions.SimOutput.System(),
    actions.SimOutput.Domains()
    ]

################################################################################
# Perform computation
################################################################################

for i, system in enumerate(systems):
        
        #Construct Actions
        actionList = []
        if mpi_comm.rank == 0:
            actionList += [actions.SimOutput(\
                hdf_file,\
                solvers[i], \
                system, \
                grids[i], \
                output_actions,\
                overwrite = True,\
                name = grids[i].name,\
                cmp_ = grids[i].comparison\
                )]
            
        # Construct and run problem
        problem = ibvp.IBVP(solvers[i], system, grid = grids[i],\
                maxIteration = 1000000, action = actionList,\
                minTimestep = 1e-12)
        problem.run(tstart, tstop)