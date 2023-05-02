#!/usr/bin/env python
# encoding: utf-8
"""
This module implements a number of differential operators based on fourier
transforms.
"""


from builtins import range
from builtins import object
import cmath, math
import numpy as np
import scipy
import pyfftw as fftw
from scipy import fftpack
from coffee.backend import backend as be


class FFT_lagrange(object):
    """Manually calculates the derivative using a cached matrix whose
    values are computed assuming application of a fourier transform
    to calculate a derivative."""

    name = "FFT_lagrange"

    def __init__(self, num_grid_points, period):
        """Initialisation for this FFT differential operator.

        Parameters
        ---------
        num_grid_points: int
            The number of grid points that the function will have.
        period: float
        """
        N = num_grid_points
        assert N % 2 == 0
        """The derivative matrix is calculated using the 'Negative Sum Trick'."""
        M = be.empty((N, N))
        for i in range(N):
            M[i, i] = 0
            for j in range(N):
                if j != i:
                    M[i, j] = (
                        (1 / 2)
                        * (-1) ** (i + j)
                        * (1 / be.tan(((i - j) * be.pi) / (N)))
                    )
        for i in range(N):
            M[i, i] = be.sum(sorted(M[i, :]))
        self.M = M * (2 * be.pi) / period

    def __call__(self, u, dx):
        """Returns the derivative.

        Parameters
        ---------
        u : tslice.TimeSlice
            The data to be differentiated.
        dx : float
            The spatial step size.
        """
        ru = be.empty_like(u)
        ru[:-1] = be.dot(self.M, u[:-1])
        ru[-1] = ru[0]
        return ru


class FFT(object):
    """Uses `numpy.fft.fft`, `numpy.fft.fftfreq` and `numpy.fft.ifft` to implement
    the differential operator."""

    name = "FFT"

    def __init__(self, order, period):
        """Initialisation for this FFT differential operator.

        Parameters
        ---------
        order : int
            The order of the required derivative.
        period: float
        """

        self.order = order
        self.period = period

    def __call__(self, u, dx):
        """Returns the derivative.

        Parameters
        ---------
        u : tslice.TimeSlice
            The data to be differentiated.
        dx : float
            The spatial step size.
        """
        # transform into fourier space
        ufft = be.fft.fft(u)
        # collect frequencies at each index
        ufreq = be.fft.fftfreq(ufft.size, d=dx)
        # compute derivative coefficient for each frequency
        dufreq = (2 * be.pi * 1j * ufreq) ** (self.order)
        # compute fourier domain derivatives
        dufft = dufreq * ufft
        # transform into 'normal' space
        rdufft = be.fft.ifft(dufft)
        return rdufft


class RFFT(object):
    """Uses `numpy.fft.rfft`, `numpy.fft.fftfreq` and `numpy.fft.irfft` to implement
    the differential operator."""

    name = "RFFT"

    def __init__(self, order, period):
        """Initialisation for this FFT differential operator.

        Parameters
        ---------
        order : int
            The order of the required derivative.
        period: float
        """
        self.order = order
        self.period = period

    def __call__(self, u, dx):
        """Returns the derivative.

        Parameters
        ---------
        u : tslice.TimeSlice
            The data to be differentiated.
        dx : float
            The spatial step size.
        """
        # get length of array
        n = u.shape[0]
        # transform into fourier space
        ufft = be.fft.rfft(u)
        # collect frequencies at each index
        ufreq = be.fft.fftfreq(ufft.size, d=dx)
        # compute derivative coefficient for each frequency
        dufreq = (2 * cmath.pi * cmath.sqrt(-1) * ufreq) ** (self.order)
        # compute fourier domain derivatives
        dufft = dufreq * ufft
        # transform into 'normal' space
        rdufft = be.fft.irfft(dufft, n)
        return rdufft


class FFT_diff_scipy(object):
    """Uses scipy.fftpack.diff to calculate the derivative."""

    name = "FFT_diff_scipy"

    def __init__(self, order, period):
        """Initialisation for this FFT differential operator.

        Parameters
        ---------
        order : int
            The order of the required derivative.
        period: float
        """
        self.order = order
        self.period = period

    def __call__(self, u, dx):
        """Returns the derivative.

        Parameters
        ---------
        u : tslice.TimeSlice
            The data to be differentiated.
        dx : float
            The spatial step size.
        """
        du = be.empty_like(u)
        du = fftpack.diff(u, self.order, self.period)
        return du


class FFTW(object):
    """Uses the fftw3 pacakge to calculate the derivative."""

    name = "FFTW"

    def __init__(self, order, period=None, fftw_flags=["estimate"]):
        self.order = order
        self.period = period
        self.fftw_flags = fftw_flags
        self.u = None
        self.dufreq = None
        # self.log = logging.getLogger("FFTW3")

    def __call__(self, u, dx):
        """returns the derivative.

        Parameters
        ---------
        u : tslice.timeslice
            the data to be differentiated.
        dx : float
            the spatial step size.
        """
        if self.u is None:
            self.u = be.empty_like(u, dtype=be.dtype(be.complex128))
            self.ufft = be.empty_like(self.u)
            self.dufft = be.empty_like(self.u)
            self.du = be.empty_like(self.u)
            self.fft = fftw.Plan(self.u, self.ufft, direction="forward")
            self.ifft = fftw.Plan(self.dufft, self.du, direction="backward")
        self.fft.guru_execute_dft(u.astype(be.dtype(be.complex128)), self.ufft)
        if self.dufreq is None:
            ufreq = be.array(
                [
                    self._compute_freq(i, self.ufft.shape[0], dx[0])
                    for i in range(self.ufft.shape[0])
                ]
            )
            self.dufreq = be.power(2 * be.pi * 1j * ufreq, self.order)
        self.dufft = self.dufreq * self.ufft
        self.ifft.guru_execute_dft(self.dufft, self.du)
        return self.du / u.shape[0]

    def _compute_freq(self, index, size, dx):
        if size % 2 == 0:
            mid = size / 2
        else:
            mid = (size - 1) / 2 + 1
        if index < mid:
            rfreq = index
        else:
            rfreq = -(size - index)
        return rfreq / (size * dx)


class FFTW_real(object):
    """Uses the fftw3 real fourier transform to calculate the derivative."""

    name = "FFTW_real"

    def __init__(self, order, period=None, fftw_flags=["estimate"]):
        """Initialisation for this FFT differential operator.

        Parameters
        ----------
        order : int
            The order of the required derivative.
        period: float, Optional
        fftw_flags: list of strings, Optional
        """
        self.order = order
        self.period = period
        self.fftw_flags = fftw_flags
        self.u = None
        self.dufreq = None

    def __call__(self, u, dx):
        """returns the derivative.

        Parameters
        ---------
        u : tslice.timeslice
            the data to be differentiated.
        dx : float
            the spatial step size.
        """
        if self.u is None:
            self.u = be.empty_like(u)
            self.ufft = be.empty(
                (math.floor((u.shape[0] / 2) + 1),), dtype=be.dtype(be.complex128)
            )
            self.dufft = be.empty_like(self.ufft)
            self.du = be.empty_like(u)
            self.fft = fftw.Plan(self.u, self.ufft, direction="forward")
            self.ifft = fftw.Plan(self.dufft, self.du, direction="backward")
        self.fft.guru_execute_dft(u, self.ufft)
        if self.dufreq is None:
            ufreq = be.array(
                [
                    self._compute_freq(i, self.ufft.shape[0], dx[0])
                    for i in range(self.ufft.shape[0])
                ]
            )
            self.dufreq = be.power(2 * be.pi * 1j * ufreq, self.order)
        self.dufft = self.dufreq * self.ufft
        self.ifft.guru_execute_dft(self.dufft, self.du)
        return self.du / u.shape[0]

    def _compute_freq(self, index, size, dx):
        if size % 2 == 0:
            mid = size / 2
        else:
            mid = (size - 1) / 2 + 1
        if index < mid:
            rfreq = index
        else:
            rfreq = -(size - index)
        return rfreq / (size * dx)


class FFT_scipy(object):
    """Uses `scipy.fftpack.fft` and `scipy.fftpack.ifft` to calculate the derivative."""

    name = "FFT_scipy"

    def __init__(self, order, period=None):
        """Initialisation for this FFT differential operator.

        Parameters
        ---------
        order : int
            The order of the required derivative.
        period: float, Optional
        """
        self.order = order
        self.period = period

    def __call__(self, u, dx):
        """Returns the derivative.

        Parameters
        ---------
        u : tslice.TimeSlice
            The data to be differentiated.
        dx : float
            The spatial step size.
        """
        ufft = fftpack.fft(u)
        ufreq = fftpack.fftfreq(ufft.size, d=dx)
        # self.log.debug("ufft = %s"%repr(ufft))
        # self.log.debug("ufreq = %s"%repr(ufreq))
        dufreq = be.power(2 * be.pi * 1j * ufreq, self.order)
        # self.log.debug("dufreq = %s"%repr(dufreq))
        dufft = dufreq * ufft
        # self.log.debug("dufft = %s"%repr(dufft))
        du = fftpack.ifft(dufft)
        # self.log.debug("du = %s"%repr(du))
        return du


class RFFT_scipy(object):
    """Uses the scipy fftpack real routines to calculate the derivative.

    The scipy implementation of RFFT doesn't seem to behave as well as
    the numpy RFFT implementation I suggest using FFT_diff_scipy instead."""

    name = "RFFT_scipy"

    def __init__(self, order, period):
        """Initialisation for this FFT differential operator.

        Parameters
        ---------
        order : int
            The order of the required derivative.
        period: float
        """
        self.order = order

    def __call__(self, u, dx):
        """Returns the derivative.

        Parameters
        ---------
        u : tslice.TimeSlice
            The data to be differentiated.
        dx : float
            The spatial step size.
        """
        n = u.shape[0]
        ufft = fftpack.rfft(u)
        ufreq = fftpack.rfftfreq(ufft.size)
        dufreq = (2 * be.pi * 1j * ufreq) ** (self.order)
        dufft = dufreq * ufft
        rdufft = fftpack.ifft(dufft, n)
        return rdufft / (dx**self.order)
