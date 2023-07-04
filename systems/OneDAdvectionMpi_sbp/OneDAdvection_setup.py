#!/usr/bin/env python
# encoding: utf-8


"""
OneDAdvection_setup.py

Created by Chris Stevens 2023
"""

# Import python libraries
from builtins import range
import sys
import os
import numpy as np
import h5py
import argparse
from mpi4py import MPI

# Import standard code base and initialize NumPy backend
from coffee.settings import init
from coffee.backend import backend as be
my_backend = be.set_backend("numpy")
init(my_backend)

from coffee import ibvp, actions, solvers, grid
from coffee.diffop.sbp import sbp

# Import system to use
import OneDAdvection
import OneDAdvection_plotter as plotter

np.set_printoptions(threshold=np.inf, precision=16)

################################################################################
# Parser settings 
################################################################################

# Initialise parser
parser = argparse.ArgumentParser(description=\
"""This program numerically solves an IBVP for a simple advection equation.""")

# Parse files
parser.add_argument('-f','-file', help=\
"""The name of the hdf file to be produced.""")
parser.add_argument('-d','-display', default = False, 
    action='store_true', help=\
"""A flag to indicate if visual display is required.""")
args = parser.parse_args()

################################################################################  
# These are the commonly altered settings
################################################################################

# Output settings
store_output = args.f is not None
display_output = args.d
if store_output and args.f is None:
    print("OneDAdvection_setup.py: error: argument -f/-file is required")
    sys.exit(1)

# File settings
if store_output:
    args.f = os.path.abspath(args.f)
    try:
        if not os.path.exists(os.path.split(args.f)[0]):
            os.makedirs(os.path.split(args.f)[0])
    except OSError as oserror:
        if oserror.errno is not errno.EEXIST:
            raise oserror

# How many systems?
num_of_grids = 1

# How many grid points?
N = 200

# What grid to use?
xstart = -1
xstop  = 1

# Times to run between
tstart = 0.0
tstop  = 10.

# Configuration of System
CFL = 0.5
tau = 1.0

# Differential operator
diffop = sbp.D43_Strand(sbp.BOUNDARY_TYPE_GHOST_POINTS)
# diffop = sbp.D42(sbp.BOUNDARY_TYPE_GHOST_POINTS)
# diffop = sbp.D43_Tiglioetal(sbp.BOUNDARY_TYPE_GHOST_POINTS)
# diffop = sbp.D43_CNG(sbp.BOUNDARY_TYPE_GHOST_POINTS)

################################################################################      
## MPI set up                                                                         
################################################################################ 

dims     = MPI.Compute_dims(MPI.COMM_WORLD.size, [0])                                    
periods  = [0]                                                                        
reorder  = True                                                                       
mpi_comm = MPI.COMM_WORLD.Create_cart(dims, periods=periods, reorder=reorder)        

################################################################################
# Grid construction
################################################################################

# Grid point data      
res = [N*2**i for i in range(num_of_grids)]

# Determine the boundary data
ghost_points = (diffop.ghost_points(),)
internal_points = (diffop.internal_points(),)
b_data = grid.MPIBoundary(
    ghost_points, 
    internal_points, 
    mpi_comm=mpi_comm, 
    number_of_dimensions=1
)

# Build grids
grids = [
    grid.UniformCart(
        (res[i],), 
        [[xstart,xstop]],
        comparison = res[i],
        mpi_comm = mpi_comm,
        boundary_data=b_data
    ) 
    for i in range(num_of_grids)
]

################################################################################
# Initialise systems
################################################################################

char_speed = 1. # Must be positive
systems = []
for i in range(num_of_grids):
    systems += [OneDAdvection.OneDAdvection(\
        diffop, tau, CFL, char_speed
        )]
    
# Configuration of IBVP
solvers = [solvers.RungeKutta4(sys) for sys in systems]
maxIteration = 10000000

################################################################################
# Set up hdf file to store output
################################################################################

if store_output and mpi_comm.rank==0:
    hdf_file = h5py.File(args.f,"w")

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
        if display_output and mpi_comm.rank == 0:
            actionList += [
                plotter.Plotter(
                    system,
                    frequency=1,
                    xlim=(xstart, xstop),
                    ylim=(-1, 1),
                    delay=0.0001
                )
            ]
        if store_output and mpi_comm.rank == 0:
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