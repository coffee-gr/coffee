#!/usr/bin/env python
# encoding: utf-8
"""
The abstract base class (abc) for System objects. 

This class serves to specify
the interface that System classes are assumed to have. These assumptions are
made by virtually all classes that call the system. For example, ibvp.IBVP,
subclasses of solver.Solver, almost all actions.hdf_actions.SimOutType
subclasses.

While you don't need to actually implement all the methods below, it really is
recommended. You may get errors about the methods missing as a consequence.

system.System has abc.ABCMeta as a meta class. This serves as a warning to
really think about what you are doing if you don't implement all methods below.
"""

from abc import ABCMeta


#############################################################################
class System(object, metaclass=ABCMeta):
    """The System class that specifies the interface for all other system
    classes.

    You don't need to implement everything below, but you can expect
    errors if you don't. The main purpose of this class is to document
    the API that classes that describe a system of equations are expected
    to provide.
    """

    def timestep(self, tslice):
        """Returns a number giving the time step for the next iteration.

        Parameters
        ----------

        tslice : tslices.TimeSlice
                A time slice with the current data.

        Returns
        -------

        float
                A number giving the dt for which the next data will be
                calculated in the next iteration.
        """
        return NotImplementedError("You need to implement this function")

    def evaluate(self, t, tslice):
        """Returns the data needed by the solver.Solver subclass needed to
        calculate the values of the functions at the next time step.

        For example, in a IBVP problem using an RK sover this method will
        return the values of the derivatives of the functions being evolved.

        Yes I know it all sounds a bit vague. But really, this method is the
        heart and soul of coffee. You need to know what you're doing, how the
        solver is implemented and what data the system is meant to pass to the
        solver.

        Parameters
        ----------

        t : float
            The time at which the 'derivatives' need to be calculated. NOTE
            THAT THIS TIME MAY BE DIFFERENT FROM THE TIME STORED IN u. The time
            stored in u is the time for the values of the functions stored in
            u, not the time at which the derivatives of the functions stored in
            u are needed. In the case of RK methods t and u.time can be
            different. YOU HAVE BEEN WARNED.

        tslice : tslice.TimeSlice
            Contains the data from which the
            'evaluate' method is meant to calculate data from.

        Returns
        -------

        tslice.TimeSlice
            data needed for the solver.Solver subclass to calculate the
            values of the functions being evolved at the next time slice.
        """
        return NotImplementedError("You need to implement this function")

    def initialValues(self, t, grid):
        """Returns initial data for the simulation.

        Parameters
        ----------

        t : float
            The time at which the initial data is to be calcualted.
        grid : ABCGrid
            The grid over which the data is to be calculated.

        Returns
        -------
        tslice.TimeSlice
            The initial values for the functions to be evolved.
        """
        return NotImplementedError("You need to implement this function")

    def left(self, t):
        """Returns the boundary data on the 'left' boundary.

        Yes this also makes the assumption that left is well defined. See the
        documentation for free_data.FreeData. To be honest I'm not sure that
        this method needs to be defined here. After all: does any class other
        than the system class need to access this method? If not there there is
        no need for it.

        As no grid object is passed as an argument the boundary value must
        be independent of current function data. This is an obvious defect.
        Also: this doesn't seem to be a problem: boundary data can be manually
        handled in the evaluate routine. Still this issue should be changed
        when it actually becomes a problem.

        Parameters
        ----------

        t : float
            The time for which the left data needs to be calculated.

        Returns
        -------

            The boundary data.
        """
        return NotImplementedError("You need to implement this function")

    def right(self, t):
        """Returns the boundary data on the 'right' boundary.

        See the documentation for the function system.System.left

        Parameters
        ----------

        t : float
            The time for which the left data needs to be calculated.

        Returns
        -------

            The boundary data.

        """
        return NotImplementedError("You need to implement this function")
