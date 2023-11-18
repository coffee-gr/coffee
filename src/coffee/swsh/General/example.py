#!/usr/bin/env python
# encoding: utf-8
"""
An example of the SWSH functionalities.
"""

# Import standard code base and initialize NumPy backend
from coffee.settings import init
from coffee.backend import backend as be
my_backend = be.set_backend("numpy")
init(my_backend)

# Import SWSH class
from coffee.swsh.General import general_swsh

# Number of \theta and \phi grid points
Ntheta = 10
Nphi   = 20

# Instantiate the class
sf = general_swsh(Ntheta, Nphi, 5)

# Create the grid
theta = sf.THETA
phi   = sf.PHI

# Define a simple test function, say 0Y10 and something else
fn1 = -0.5*be.sqrt(1.5 / be.pi)*be.sin(theta)*be.exp(1.j*phi)
fn2 = -0.5*be.sqrt(3. / be.pi)*(be.cos(0.5*theta)**2.)*be.exp(1.j*phi)
fn3 = 0.5*be.sqrt(5. / be.pi)*(be.sin(0.5*theta)**2.)*be.sin(theta)*be.exp(2.j*phi)

# Convert this to a spin-weighted function (sf)
# fns_sw = sf.function_in_s2([fn1,fn2,fn3], [0,-1,1])
fns_sw = sf.function_in_s2([fn1,fn2,fn3], [0,-1,1])

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