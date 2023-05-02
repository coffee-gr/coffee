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
selection, the Imm and Jmm routines, the even odd method for real data
or the multi versions of these functions.
Support for this could be added via ctypes or editting the 
python/spinsfast_module.c code.

The failure to wrap the multi functions requires this module to
perform some array manipulations before calling spinsfast.
"""

import numpy as np
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(
    0, os.path.join(dir_path, "huffenberger_wandelt", "spinsfast_rev104", "lib")
)
import spinsfast
from coffee.swsh.spinsfastpy import salm


def forward(f, spins, lmax):
    """
    Returns the spin spherical harmonic coefficients for f.

    Parameters
    ----------
    f : numpy.ndarray of shape (Nmap, Ntheta, Nphi)
         the values of a function on S^2
         according to the ecp discretisation, see
         Section 2.3
    spins : float, numpy.ndarray of shape (Nmap,) of floats
        Specifies
        the required spins for the calculation. Note that allowed
        floats are signed half integers. This is not checked.
    lmax : float
        A float representation of an integer or half-integer. It specifies
        the maximum bandwidth.

    Returns
    -------
    salm : swsh.salm.salm
        Represents the spin weighted spherical harmonic expansion of f.
    """
    if len(f.shape) == 2:
        f = be.array([f])
    Nmaps, Ntheta, Nphi = f.shape
    spins = be.atleast_1d(spins)
    if len(spins.shape) != 1:
        raise ValueError("spins must be an int or a one dimensional array of ints")
    Ntransform = spins.shape[0]
    if Nmaps != Ntransform:
        raise ValueError("number of spins must equal the number of functions")

    alm = be.array([spinsfast.map2salm(f[i], spins[i], lmax) for i in range(Nmaps)])
    return salm.sfpy_salm(alm, spins, lmax)


# a = be.zeros((2,8,8), dtype="complex")
# print forward(a, [1,3], 3)
