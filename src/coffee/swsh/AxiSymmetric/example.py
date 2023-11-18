#!/usr/bin/env python
# encoding: utf-8
"""
An example of the axi-symmetric SWSH functionalities.
"""

# Import standard code base and initialize NumPy backend
from coffee.settings import init
from coffee.backend import backend as be
my_backend = be.set_backend("numpy")
init(my_backend)

# Import SWSH class
from coffee.swsh.AxiSymmetric import axi_symmetric_swsh
import numpy as np

# Number of \theta grid points
Ntheta = 5

# Maximum spin-weight
smax = 1

# Instantiate the class
sf = axi_symmetric_swsh(Ntheta, smax)

# Create the grid
theta = sf.create_mesh(Ntheta)

# Define a simple test function, say 0Y10 and something else
fn1 = 0.5*np.sqrt(3. / np.pi)*np.cos(theta)
fn2 = np.sin(theta)

# Convert this to a spin-weighted function (sf)
fns_sw = sf.function_in_s2([fn1,fn2],[0,-1])

# Forward transform to calculate the spectral coefficients
fns_salm = sf.forward(fns_sw)

# Print them out
print(fns_salm)

# Compute the backward transform of the spectral coefficients
fns_sw_alt = sf.backward(fns_salm)

# Print out the difference to original function
print("|fns - backward(forward(fns))|")
print(fns_sw_alt - fns_sw)

# Take an eth-derivative
Deth_fns = sf.ethU(fns_sw)

# Print it out
print("Deth_fns")
print(Deth_fns)

# Take an eth'-derivative
Dethp_fns = sf.ethD(fns_sw)

# Print it out
print("Dethp_fns")
print(Dethp_fns)