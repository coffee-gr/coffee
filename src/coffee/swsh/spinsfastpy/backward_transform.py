""" 
This module wraps the Huffenberger and Wandelt spinsfastpy python module. It
does this to provide an additional level of indirection to make it easier to
update the C spinsfast code and associated python module if needed.
It currently is based on revision 104
of the forward transformations
implemented in Huffenberger's & Wandelt's spinsfast code. See "Fast and exact 
spin-s Spherical Harmonic Transformations" in Astrophysical Journal Supplement 
Series (2010) 189:255-260. All references to sections, appendices or equations 
in the documentation below refer to equations given in this paper.
The spinsfast code website is
http://astrophysics.physics.fsu.edu/~huffenbe/research/spinsfast/index.html

The revision 104 python module does not provide access to delta method
selection and the Gmm routines.
Support for this could be added via ctypes or editting the 
python/spinsfast_module.c code.
"""

import numpy as np
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(
    0, os.path.join(dir_path, "huffenberger_wandelt", "spinsfast_rev104", "lib")
)
import spinsfast


def backward(salm, Ntheta, Nphi):
    """
    Returns a function on S^2 given it's spin coefficients in the form of an
    salm object.

    Parameters
    ----------
    salm : swsh.salm.salm
        a subclass of numpy.ndarray that stores the spin coefficients
    Ntheta : int
        an int giving the number of points discritising the theta variable of S^2
    Nphi : an int giving the number of points discritising the phi variable of S^2

    Returns
    -------
    numpy.ndarray :
         a numpy.ndarray of shape (salm.spins.shape[[0], Ntheta, Nphi)
         containing the values of the
         function on S^2, parameterised via the ecp discretisation, see
         section 2.3
    """
    spins = be.atleast_1d(salm.spins)
    Ntransform = spins.shape[0]
    lmax = salm.lmax
    if len(spins.shape) != 1:
        raise ValueError("spins must be an int or a one dimensional array of ints")
    data = be.asarray(salm)

    f = be.array(
        [
            spinsfast.salm2map(data[i], spins[i], lmax, Ntheta, Nphi)
            for i in range(Ntransform)
        ]
    )
    return f
