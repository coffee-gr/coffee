#!/usr/bin/env python
# encoding: utf-8
"""
The FreeData class describes the \`free data' of a system of equations.

This is a subclass of Abstract Base Class. It specifies an interface. All
classes that are used to calculate initial data or free data for a simulation
scheme are assumed to implement the methods below. In particular the object
actions.hdf_output.SimOutputType.Exact uses the exact method.

It is not, in general, needed to implement all of the following methods. What
you do implement should be determined by what additional objects interact with
your subclass of the system.System class.

This object is subclassed from abc inorder to make you think twice about what
you are doing. If your free data is sufficiently simple the code that
corresponds to free data could be included directly into system.System since
all calls to FreeData are mediated through system.System.
"""

import abc
from abc import ABCMeta

class FreeData(object, metaclass=ABCMeta):
    """This is a simple abstract base class for the implementation of exact
    solutions.

    All the 'virtual' functions must be defined, modulo the comments given in
    the module docstring.
    """

    @abc.abstractmethod
    def left_boundary(self, tslice):
        """Returns the boundary values given on the left boundary over a one dimensional grid.

        The currect FreeData object assumes that the simulation is being
        done over a one dimensional grid. Hence the call to 'left_boundary'.
        this method will eventually be altered to bring it into line with the
        grid method, \`boundary_slices' which specifies the slices of the domain
        corresponding to the boundary in a mannor which is independent of
        assumptions about the number of dimensions.

        Until then you'll just have to cope.

        Parameters
        ----------

        tslice : tslice.TimeSlice
            The timeslice of the data which left boundary values need to given.

        Returns
        -------

            The values of the functions being numerically simulated on
            left boundary. The type of these values is up to you,
            inconjunction with your implementation of the tslices.TimeSlice
            and system.System classes.
        """
        raise NotImplementedError("You must define the function left_boundary().")

    @abc.abstractmethod
    def right_boundary(self, tslice):
        """Returns the boundary values given on the right boundary over a one
        dimensional grid.

        Arguments and returned object are as for left_boundary.

        Parameters
        ----------

        tslice : tslice.TimeSlice
            The timeslice of the data which left boundary values need to given.

        Returns
        -------

            The values of the functions being numerically simulated on
              left boundary. The type of these values is up to you,
              inconjunction with your implementation of the tslices.TimeSlice
              and system.System classes.

        """
        raise NotImplementedError("You must define the function right_boundary().")

    @abc.abstractmethod
    def initial_data(self, t, grid):
        """Returns the initial data for a simulation based on the initial
        time and grid.

        Parameters
        ----------

        t : float
            The initial time at which the initial data is to be calculated.

        grid :
            The grid over which the initial data is to be calculated.

        Returns
        -------

        tslice.TimeSlice
            the initial data of the simulation.
        """
        raise NotImplementedError("You must define the function initial_data().")

    @abc.abstractmethod
    def exact(self, t, grid):
        """Returns the exact values of the functions being numerically found.

        The idea is that this data can be calculated and stored in the hdf file
        during the simulation as a result of a call to the
        actions.hdf_actions.SimOutputType.Exact object.

        Of course you don't have to use this method that way...

        Parameters
        ----------

        t : float
            The time at which to calculate the exact values.
        grid :
            The domain over which to calculate the exact values.

        Returns
        -------

            The exact values of the functions being solved for in the simulation.

        """
        raise NotImplementedError("You must define the function exact().")

    @abc.abstractmethod
    def dirichlet(self, u, intStep=None):
        """Return dirichlet boundary conditions.

        This is a poorly thought through addition to this class. The idea is
        that at certain steps during the evolution of data it might be
        necessary to call dirichlet boundary conditions. This method does that.
        Why is something different from left_boundary needed? I don't know.
        That's why it's poorly thought through. The bottom line is that this
        method might be removed from this class at some point.

        Originally this method was used to implement the intermediate boundary
        conditions for explicit RK routines during evolution. In order to do
        this knowledge of which intermediate step the boundary condition was
        being applied is necessary. This is the reason for the keyword
        argument.

        Parameters
        ----------

        u : tslice.TimeSlice
            Stores data for which dirichlet boundary conditions are needed.

        intStep : int
            Giving which intermediate step the RK routine is currently in.

        """
        raise NotImplementedError("You must define the function dirichlet().")
