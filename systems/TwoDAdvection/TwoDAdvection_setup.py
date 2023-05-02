# import python libraries
from builtins import range
import os
import time
import h5py
from mpi4py import MPI
import argparse

# Import standard code base
from coffee import ibvp, actions, solvers, grid
from coffee.diffop.sbp import sbp

# import system to use
import TwoDAdvection
import TwoDAdvection_plotter

from coffee.backend import backend as be

be.set_backend("numpy")

import jax.numpy as jnp
from jax import grad, jit, vmap
from jax import random

################################################################################
# Argument Parsing and File Nameing
################################################################################
# Initialise parser
parser = argparse.ArgumentParser(
    description="""This script contains the settings for running SpinWave simulations."""
)

# Parse files
parser.add_argument("-f", "-file", help="""The name of the hdf file to be produced.""")
parser.add_argument(
    "-d",
    "-display",
    default=False,
    action="store_true",
    help="""A flag to indicate if visual display is required.""",
)
parser.add_argument(
    "-i",
    "-info",
    default=False,
    action="store_true",
    help="""A flag to indicate if information about progress of simulation is
required.""",
)
parser.add_argument(
    "-tstop",
    type=float,
    default=1.0,
    help="""The time that the
        simulation will run to. Defaults to 1.0.""",
)
parser.add_argument(
    "-tstart",
    type=float,
    default=0.0,
    help="""The initial time
        of the simulation. Defaults to 0.0.""",
)
parser.add_argument(
    "-nsim",
    type=int,
    default=1,
    help="""The number of
simulations that you want to run. Defaults to 1.""",
)

# Actually parse the arguments
args = parser.parse_args()

###############################################################################
# File and logging settings
###############################################################################
# output settings
store_output = args.f is not None
display_output = args.d

# file settings
if store_output:
    args.f = os.path.abspath(args.f)
    try:
        if not os.path.exists(os.path.split(args.f)[0]):
            os.makedirs(os.path.split(args.f)[0])
    except OSError as oserror:
        if oserror.errno is not errno.EEXIST:
            raise oserror

###############################################################################
# MPI set up
###############################################################################
dims_list = [0, 0]
dims = MPI.Compute_dims(MPI.COMM_WORLD.size, dims_list)
periods = [1, 1]
reorder = True
mpi_comm = MPI.COMM_WORLD.Create_cart(dims, periods=periods, reorder=reorder)

################################################################################
# System, IBVP and Grid settings
################################################################################

# How many systems?
num_of_grids = args.nsim

# How many grid points?
Nx = 100
Ny = 100

# What grid to use?
xstart = 1
xstop = 3
ystart = 0
ystop = 2

# Times to run between
tstart = 0.0
tstop = 50.0

# Configuration of System
CFLs = [0.5 for i in range(num_of_grids)]
tau = 1.0

# 1st derivative x
xaxis_1D_diffop = sbp.D43_Strand()
# raxis_1D_diffop = fd.FD12()
# raxis_1D_diffop = fd.FD14()
# raxis_1D_diffop = sbp.D42(log)
# raxis_1D_diffop = fft.FFT_diff_scipy(1,rstop-rstart)
# raxis_1D_diffop = fft.FFT(1,xstop-xstart)
# raxis_1D_diffop = fft.RFFT(1)
# raxis_1D_diffop = fft.FFT_scipy(1,xstop-xstart)
# raxis_1D_diffop = fft.FFT_lagrange1(N,xstop-xstart)

# 1nd derivative y
yaxis_1D_diffop = sbp.D43_Strand()

# Configuration of IBVP
maxIteration = 1000000

################################################################################
# Grid construction
################################################################################

# Grid point data
axes = [(Nx * 2**i, Ny * 2**i) for i in range(num_of_grids)]

# Determine the boundary data
ghost_points = (xaxis_1D_diffop.ghost_points(), yaxis_1D_diffop.ghost_points())
internal_points = (
    xaxis_1D_diffop.internal_points(),
    yaxis_1D_diffop.internal_points(),
)
b_data = grid.MPIBoundary(
    ghost_points, internal_points, mpi_comm=mpi_comm, number_of_dimensions=2
)

# Build grids
grids = [
    grid.UniformCart(
        axes[i],
        [(xstart, xstop), (ystart, ystop)],
        comparison=i,
        mpi_comm=mpi_comm,
        boundary_data=b_data,
    )
    for i in range(num_of_grids)
]

################################################################################
# Initialise systems
################################################################################
systems = []
for i in range(num_of_grids):
    systems += [
        TwoDAdvection.TwoDadvection(
            1,
            1,
            xaxis_1D_diffop,
            yaxis_1D_diffop,
            CFL=CFLs[i],
            equation_coords="Cartesian",
            tau=tau,
        )
    ]

################################################################################
# Initialise Solvers
################################################################################
RKsolvers = [solvers.RungeKutta4(system) for system in systems]

################################################################################
# Set up hdf file to store output
################################################################################
if store_output and mpi_comm.Get_rank() == 0:
    hdf_file = h5py.File(args.f)

################################################################################
# Set up action types for data storage in hdf file
################################################################################
if store_output and mpi_comm.Get_rank() == 0:
    output_actions = [
        actions.SimOutput.Data(),
        actions.SimOutput.Times(),
        actions.SimOutput.System(),
        actions.SimOutput.Domains(),
    ]

################################################################################
# Perform computation
################################################################################
for i, system in enumerate(systems):
    # Construct Actions
    actionList = []
    if display_output and mpi_comm.rank == 0:
        actionList += [
            TwoDAdvection_plotter.Plotter(
                system,
                frequency=1,
                xlim=(xstart, xstop),
                ylim=(ystart, ystop),
                delay=0.0001,
                fignum=2,
            )
        ]
    if store_output and mpi_comm.Get_rank() == 0:
        actionList += [
            actions.SimOutput(
                hdf_file,
                RKsolvers[i],
                system,
                grids[i],
                output_actions,
                overwrite=True,
                name=grids[i].name,
                cmp_=grids[i].comparison,
            )
        ]
    problem = ibvp.IBVP(
        RKsolvers[i],
        system,
        grid=grids[i],
        maxIteration=maxIteration,
        action=actionList,
    )
    start_time = time.time()
    problem.run(tstart, tstop)
    stop_time = time.time()
    print(("Simulation", i, "took", stop_time - start_time, "seconds."))
