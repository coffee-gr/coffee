"""
This is the package spinsfastpy. It contains python bindings to Huffenberger's &
Wandelt's spinsfast C code. We recommend using the import statement: 
import spinsfast as sfpy.

Importing this module provides two methods, forward and backward, and one 
dictionary. Please see the individual documentation for how to use
forward and backward. 

Currently no control is provided over the generation of wigner d functions
as the code here relies on Huffenberger and Wandelt's python module api which
does not expose the required methods and arguments.

Some example code, both C and Python can be found in ./examples.
"""

from coffee.swsh.spinsfastpy.forward_transform import forward
from coffee.swsh.spinsfastpy.backward_transform import backward

del forward_transform
del backward_transform


def set_clebsch_gordan_default(cg):
    """This function sets the default method for computation of
    Clebsch Gordan coefficients in the module coffee.diffop.swsh.spinsfastpy."""
    salm.sfpy_salm.cg_default = cg
