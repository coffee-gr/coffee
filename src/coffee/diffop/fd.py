#!/usr/bin/env python
# encoding: utf-8
"""
This module contains code for the use of finite difference stencils as
differential operators.

Stencil values are calculated using Fornberg's paper. Code for the
calculations is available in
`finite_difference_weights_generator.py` located in Code/Utils.
See Fornberg's paper for the notation.
"""

from builtins import object
from coffee.settings import be

################################################################################
# Finite difference default ghost point processor
################################################################################
def ghost_point_processor(data, b_values):
    """A utility function that understands how data defined over ghost points
    for finite difference stencils should be handled when communicated via MPI."""
    for b_slice, b_data in b_values:
        data[b_slice] = b_data

################################################################################
# Finite difference stencil base class.
################################################################################


class FD_stencil(object):
    """This class takes an array of numbers and uses them to act on a
    vector of numbers.

    That is an FD_stencil represents
    a 'standard' finite difference stencil, e.g.:

        [-1, 2, -1]

    as a linear operator.
    """

    def __init__(self, stencil_coefs, loffset, roffset):
        """The initialiser for an FD_stencil.

        Parameters
        ==========
        stencil_coefs : list of floats
            Gives a finite difference stencil over a fixed step size.
        loffset : int
            How many elements of the `stencil_coefs` occurs before the
            point at which the derivative is being calculated.
        roffset: int
            Similar to loffset, but counted from the right.
        """
        self.stencil_coefs = stencil_coefs
        self.loffset = loffset
        self.roffset = roffset

    def __repr__(self):
        return "<FD_stencil %s, [%d,%d]>" % (
            repr(self.stencil_coefs),
            self.loffset,
            self.roffset,
        )

    def __call__(self, u, apply_at=None):
        """Applies the finite difference stencil to calculate a vector of
        derivatives of the function represented by the given data.

        Parameters
        ==========
        u : tslice.TimeSlice
            The data to be differentiated.
        apply_at : int, Optional
            The index of the data at which the derivative is to be calculated.
            If not given then all derivatives are calculated.

        Returns
        =======
        numpy.ndarray:
            A list of derivatives, or a list with one derivative.
        """
        if apply_at is not None:
            lbound = apply_at - self.loffset
            rbound = apply_at + self.roffset + 1
            if rbound == 0:
                rbound = None
            return be.dot(self.stencil_coefs, u[lbound:rbound])
        else:
            """The convolve method does not quite do the right thing, it
            flips the direction of iteration. Hence the ::-1 below.

            Note also that the convolve method is based on the
            multiarray.correlate a C routine."""
            return be.convolve(u, self.stencil_coefs[::-1], mode="same")


def flip(stencil):
    """This method takes, for example, a forward difference operator and returns
    a backwards difference operator using 'the same' stencil.

    For example:
        flip(F12_stencil) = B12_stencil
        flip(C12_stencil) = C12_stencil

    Parameters
    ==========
    stencil : list of floats
        A 'standard' finte difference stencil over a fixed grid step.
    """
    return FD_stencil(-stencil.stencil_coefs[::-1], stencil.roffset, stencil.loffset)


################################################################################
# First derivative stencils
################################################################################
"""Further stencil can be calculated using 
`finite_difference_weights_generator.py' located in Code/Utils.
I suggest reading Fornberg's paper for the notation."""


class C12_stencil(FD_stencil):
    """The standard central second order first derivative finite difference."""

    def __init__(self):
        stencil = be.array([-1.0 / 2.0, 0.0, 1.0 / 2.0])
        loff = 1
        roff = 1
        super(C12_stencil, self).__init__(stencil, loff, roff)


class C14_stencil(FD_stencil):
    """The standard central fourth order first derivative finite difference."""

    def __init__(self):
        stencil = be.array([1.0 / 12, -2.0 / 3, 0, 2.0 / 3, -1.0 / 12])
        loff = 2
        roff = 2
        super(C14_stencil, self).__init__(stencil, loff, roff)


class F12_stencil(FD_stencil):
    """The standard forward second order first derivative finite difference."""

    def __init__(self):
        stencil = be.array([-3.0 / 2.0, 2, -1.0 / 2.0])
        loff = 0
        roff = 2
        super(F12_stencil, self).__init__(stencil, loff, roff)


class F13_stencil(FD_stencil):
    """The standard forward third order first derivative finite difference."""

    def __init__(self):
        stencil = be.array([-11.0 / 6.0, 3.0, -3.0 / 2.0, 1.0 / 3])
        loff = 0
        roff = 3
        super(F13_stencil, self).__init__(stencil, loff, roff)


class F14_stencil(FD_stencil):
    """The standard forward fourth order first derivative finite difference."""

    def __init__(self):
        stencil = be.array([-25.0 / 12, 4.0, -3.0, 4.0 / 3, -1.0 / 4])
        loff = 0
        roff = 4
        super(F14_stencil, self).__init__(stencil, loff, roff)


class MF14_stencil(FD_stencil):
    """An offset forward fourth order first derivative finite difference
    that calculates the derivative at the second grid point."""

    def __init__(self):
        stencil = be.array([-1.0 / 4, -15.0 / 18, 3.0 / 2, -1.0 / 2, 15.0 / 180])
        loff = 1
        roff = 3
        super(MF14_stencil, self).__init__(stencil, loff, roff)


class B12_stencil(FD_stencil):
    """The standard backward second order first derivative finite difference."""

    def __init__(self):
        stencil = be.array([1.0 / 2.0, -2.0, 3.0 / 2.0])
        loff = 2
        roff = 0
        super(B12_stencil, self).__init__(stencil, loff, roff)


################################################################################
# Second derivative stencils
################################################################################


class C22_stencil(FD_stencil):
    """The standard central second order second derivative finite difference."""

    def __init__(self):
        stencil = be.array([1.0, -2.0, 1.0])
        loff = 1
        roff = 1
        super(C22_stencil, self).__init__(stencil, loff, roff)


class F22_stencil(FD_stencil):
    """The standard forward second order second derivative finite difference."""

    def __init__(self):
        stencil = be.array([2.0, -5.0, 4.0, -1.0])
        loff = 0
        roff = 3
        super(F22_stencil, self).__init__(stencil, loff, roff)


class B22_stencil(FD_stencil):
    """The standard backward second order second derivative finite difference."""

    def __init__(self):
        stencil = be.array([1.0, -4.0, 5.0, -2.0])
        loff = 3
        roff = 0
        super(B22_stencil, self).__init__(stencil, loff, roff)


################################################################################
# Finite Difference base class
################################################################################


class FD_diffop(object):
    """The FD_diffop class represents the action of a collection of
    FD_stencils on a vector of values."""

    name = "FD_diffop"

    def __init__(self):
        pass

    def __call__(self, u, dx):
        """Gives the result of the linear operator represented by this
        class on a vector of data.

        Parameters
        ==========
        u : tslice.TimeSlice
            The data to take a derivative of.
        dx: float
            The spatial step size.

        Returns
        =======
        numpy.ndarray:
            The vector of derivatives.
        """
        ru = self.central(u)
        for i, b in self.boundaries:
            ru[i] = b(u, apply_at=i)
        return ru / (dx**self.order)

    def __str__(self):
        return "Differential operator " % self.name

    def save(self):
        """This is a convenience method that saves a textual representation,
        using be.savetxt, of the data used in the finite difference operator
        to a file in the users home directory."""
        filename = os.path.expanduser("~/" + self.name)
        print(filename)
        be.savetxt(filename + "_left.txt", self.central)
        be.savetxt(filename + "_right.txt", self.boundaries)

    def ghost_points(self):
        """Returns the number of ghost points for this operator.

        See the mpi module for documentation of ghost points.

        Returns
        =======
        two tuple of ints:
            The number of ghost points for the left and right application of
            this operator.
        """
        return self.central.loffset, self.central.roffset

    def internal_points(self):
        """Returns the number of internal points for this operator.

        See the mpi module for documentation of internal points.

        Returns
        =======
        two tuple of ints:
            The number of interal points for the left and right application of
            this operator.
        """
        return self.ghost_points()


################################################################################
# First Derivative operators
################################################################################
class FD12(FD_diffop):
    """Implements a second order FD routine for first order derivatives.

    On the boundaries uses forward/backward differences. Sufficiently
    away from the boundaries uses a central scheme. Assumes that the
    FD grid is evenly spaced."""

    name = "FD12"
    central = C12_stencil()
    boundaries = be.array([(0, F12_stencil()), (-1, B12_stencil())])
    order = 1


class FD14(FD_diffop):
    """Implements a fourth order FD routine for first order derivatives.

    On the boundaries uses forward/backward differences. On
    interor points where the central difference scheme cannot be applied,
    skewed forward/backward difference stencils are used. Sufficiently
    away from the boundaries uses a central scheme. Assumes that the
    FD grid is evenly spaced."""

    name = "FD14"
    central = C14_stencil()
    boundaries = be.array(
        [
            (0, F14_stencil()),
            (1, MF14_stencil()),
            (-2, flip(MF14_stencil())),
            (-1, flip(F14_stencil())),
        ]
    )
    order = 1


################################################################################
# Second Derivative operators
################################################################################


class FD22(FD_diffop):
    """Implements a second order FD routine for second order derivatives.

    On the boundaries uses forward/backward differences. Sufficiently
    away from the boundaries uses a central scheme. Assumes that the
    FD grid is evenly spaced."""

    name = "FD22"
    central = C22_stencil()
    boundaries = be.array([(0, F22_stencil()), (-1, B22_stencil())])
    order = 2
