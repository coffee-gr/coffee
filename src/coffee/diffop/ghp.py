"""Numerical implementations of eth and eth'.

These are the differential operators, of the same name, given in the GHP spinor
formalism. To use these operators coffee.swsh is used to convert function data
into a spin weighted spherical harmonic expansion - if the data is not already
represented in that form.
"""

from builtins import range
from builtins import object
import math
import os
import numpy as np

from coffee.swsh import spinsfastpy as sfpy

# This module could be reimplemented using meta-classes as eth and ethp
# essentially just differ by a boolean value. I have already tried to avoid
# having copies of essentially the same code, but this would completely
# remove this remainin duplicates


class _ghp_operator(object):
    """A utilty class that implements both the eth and eth' operator."""

    def __init__(self, prime):
        """The initialiser for the class.

        Parameters
        ==========
        prime : bool
            If true then the operator is eth prime, otherwise the operator
            is eth.
        """
        self.prime = prime

    def __call__(self, u, spins=None, lmax=None):
        """Calculate the derivative.

        Parameters
        ==========
        u : tslice.TimeSlice
        spins : list of ints, Optional
            Required if u is not expressed as a spectral decomposiiton of
            spin weighted spherical harmonic functions.
        lmax : list of ints, Optional
            Required if u is not expressed as a spectral decomposiiton of
            spin weighted spherical harmonic functions.
        """
        if isinstance(u, sfpy.salm.sfpy_sralm):
            return self._eval_sralm(u)
        else:
            return self._eval_salm(u, spins, lmax)

    def _eval_salm(self, u, spins=None, lmax=None):
        Ntheta = None
        Nphi = None
        if spins is not None and lmax is not None:
            Ntheta, Nphi = u.shape
        salm, spins, lmax = _transform_to_harmonic_space(u, spins, lmax)
        if self.prime:
            r_salm = sfpy.salm.sfpy_salm(
                be.empty_like(salm, dtype=be.typeDict["complex"]), spins - 1, lmax
            )
            _do_derivative(salm, spins, lmax, _ethp_eigen, r_salm)
        else:
            r_salm = sfpy.salm.sfpy_salm(
                be.empty_like(salm, dtype=be.typeDict["complex"]), spins + 1, lmax
            )
            _do_derivative(salm, spins, lmax, _eth_eigen, r_salm)
        return _transform_from_harmonic_space(r_salm, Ntheta, Nphi)

    def _eval_sralm(self, u):
        if u.spins.shape is ():
            rv_temp = self._eval_salm(u[0])
            rv = be.empty((u.shape[0],) + (rv_temp.shape[0],), dtype=u.dtype)
            rv[0, :] = be.asarray(rv_temp)
            for i in range(1, u.shape[0]):
                rv[i, :] = self._eval_salm(u[i])
            return sfpy.salm.sfpy_sralm(
                rv,
                rv_temp.spins,
                rv_temp.lmax,
                cg=rv_temp.cg,
                bandlimit_multiplication=rv_temp.bl_mult,
            )
        else:
            rv_temp = self._eval_salm(u[:, 0])
            rv = be.empty(
                (rv_temp.shape[0],) + (u.shape[1],) + (rv_temp.shape[1],), dtype=u.dtype
            )
            rv[:, 0, :] = be.asarray(rv_temp)
            for i in range(1, u.shape[1]):
                rv[:, i, :] = self._eval_salm(u[:, i])
            return sfpy.salm.sfpy_sralm(
                rv,
                rv_temp.spins,
                rv_temp.lmax,
                cg=rv_temp.cg,
                bandlimit_multiplication=rv_temp.bl_mult,
            )


def _eth_eigen(s, j):
    if abs(s + 1) > j:
        return 0
    else:
        return -math.sqrt(j * (j + 1) - s * (s + 1))


def _ethp_eigen(s, j):
    if abs(s - 1) > j:
        return 0
    else:
        return math.sqrt(j * (j + 1) - s * (s - 1))


def _transform_to_harmonic_space(u, spins, lmax):
    if spins is not None and lmax is not None:
        salm = sfpy.forward(u, spins, lmax)
        spins = be.asarray(spins)
    else:
        salm = u
        spins = u.spins
        lmax = u.lmax
    return salm, spins, lmax


def _do_derivative(salm, spins, lmax, eigen_calc, r_salm):
    for j in range(lmax + 1):
        if spins.shape is ():
            sm_values = salm[j]
            eigenvalue = eigen_calc(spins, j)
            r_salm[j] = eigenvalue * sm_values
        else:
            sm_values = salm[:, j]
            eigenvalues = be.vectorize(eigen_calc)(spins, j)
            r_salm[:, j] = eigenvalues[:, be.newaxis] * sm_values


def _transform_from_harmonic_space(r_salm, Ntheta, Nphi):
    if Ntheta is not None:
        return sfpy.backward(r_salm, Ntheta, Nphi)
    return r_salm


class eth(_ghp_operator):
    """
    Implements the differential operator given by
    \eth{}_sY_{lm} = -\sqrt{l(l+1)-s(s+1)}{}_{s+1}Y_{lm}
    """

    def __init__(self):
        super(eth, self).__init__(False)

    def __repr__(self):
        return "eth"


class ethp(_ghp_operator):
    """
    Implements the differential operator given by
    \eth'{}_sY_{lm} = \sqrt{l(l+1)-s(s-1)}{}_{s-1}Y_{lm}
    """

    def __init__(self):
        super(ethp, self).__init__(True)

    def __repr__(self):
        return "ethp"
