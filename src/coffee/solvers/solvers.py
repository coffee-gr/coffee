#!/usr/bin/env python
# encoding: utf-8
"""
This module contains the abstract base class for solvers as well as several
default implementations.

This default implementations all assume that tslice.TimeSlice objects can be
operated on by numerically (e.g. addition, multiplication, subtraction).
"""


from builtins import object
import abc

from coffee.tslices import tslices


################################################################################
# Solver Abstract Base Class
################################################################################
class ABCSolver(object):
    """The ABSolver abstract base class.

    This class provides the interface that is assumed of all other solvers.
    These methods are called in the ibvp.IBVP class.
    """

    def __init__(self, system, *args, **kwds):
        """Returns an instance of ABCSolver.

        This method should be called by any subclasses. It defines the system
        atribute and sets up a logging object.

        Parameters
        ----------

        system : system
            The system being used.
        """
        self.system = system
        super(ABCSolver, self).__init__(**kwds)

    @abc.abstractproperty
    def name(self):
        """This should return a name associated to the subclass.

        The name will need to be defined in the subclasses constructor. The
        name is used to identify the instantiated class in the logger.

        Returns
        -------

        string
            A string name for the class. This defaults to ABCSolver.

        """
        return "ABCSolver"

    @abc.abstractmethod
    def advance(self, t, u, dt):
        """Returns a tslice.TimeSlice containing data at time t + dt.

        Parameters
        ----------

        t : float
            The current time.

        u : tslice.TimeSlice
            The current tslice.TimeSlice. Note that the time in the timeslice
            is the time for the data in the timeslice, but not necessarily
            equal to t.

        dt : float
            The step in time required.

        Returns
        -------

        tslice.TimeSlice
            Contains the data at time t+dt.

        """
        raise NotImplementedError("This method needs to be implemented")

    def __repr__(self):
        return "<%s: system = %s>" % (self.name, self.system)


################################################################################
# First order methods
################################################################################
class Euler(ABCSolver):
    """An implementation of the first order Euler method."""

    name = "Euler"

    def advance(self, t, u, dt):
        """An Euler method implementation of ABCSolver.

        See the ABCSolver.advance method for documentation.
        """
        du = self.system.evaluate(t, u)
        r_time = t + dt
        r_slice = u + dt * du
        r_slice.time = r_time
        return (r_time, r_slice)


class ImplicitEuler(ABCSolver):
    """An implementation of the first order implicit Euler method."""

    name = "ImplicitEuler"

    def _NR(self, u_start, u_next, t, dt, domain):
        """A Newton-Raphson auxilliary routine for
        the implementation of the implicit Euler scheme.

        Parameters
        ----------
        u_start: numpy.ndarray
            Initial values of the function.

        u_next: numpy.ndarray
            The first guess at the next values of the function.

        t: float
            Current time.

        dt: float
            Current time step.

        domain: grid
            The domain that the function is being evaluated over.

        """

        # Starting guess for Newton-Raphson algorithm. Currently
        # chooses current tslice.
        # To be perhaps moved so user specifies this in system file or
        # perhaps in setup file when choosing implicit Euler.
        prev_u_next = u_next

        while True:
            # Evolution equation of u
            g = self.system.evaluate(
                t + dt, tslices.TimeSlice(prev_u_next, domain, time=t)
            ).data

            # Derivative of evolution equation of u w.r.t. evolved variable
            dg, TOL = self.system.implicit_method_jacobian(
                t + dt, tslices.TimeSlice(prev_u_next, domain, time=t)
            )

            dg = dg.data

            # NR is applied to f
            f = prev_u_next - u_start - dt * g
            df = 1.0 - dt * dg

            # Compute next value of u_next
            next_u_next = prev_u_next - (f / df)

            # If we get to within a certain accuracy of the root stop
            if all(abs(next_u_next - prev_u_next)[0] < TOL):
                break

            # Otherwise continue
            prev_u_next = next_u_next

        return next_u_next

    def advance(self, t, u, dt):
        """An implicit Euler method implementation of ABCSolver.

        See the ABCSolver.advance method for documentation.
        """

        # Current values of the evolved variables
        u_start = u.data

        # Initial guess for the Newton-Raphson method is the current value(s)
        # of the evolved variables
        u_next = u_start

        # prev_u_next = u_next
        prev_u_next = u_next

        # Returns u_next from applying Newton-Raphson
        u_ret = self._NR(u_start, u_next, t, dt, u.domain)

        r_time = t + dt
        r_slice = tslices.TimeSlice(u_ret, u.domain, r_time)
        return (r_time, r_slice)


################################################################################
# Fourth order methods
################################################################################
class RungeKutta4(ABCSolver):
    """A RungeKutta4 implementation of ABCSolver."""

    name = "RK4"

    def __init__(self, *args, **kwds):
        super(RungeKutta4, self).__init__(*args, **kwds)

    def advance(self, t0, u0, dt):
        """See the ABCSolver.advance method for documentation.

        Very simple minded implementation of the standard 4th order Runge-Kutta
        method to solve an ODE of the form fdot = rhs(t,f)
        """
        eqn = self.system
        u = u0
        t = t0
        k = eqn.evaluate(t, u)
        u1 = u0 + (dt / 6.0) * k

        u = u0 + dt / 2.0 * k
        t = t0 + dt / 2.0
        k = eqn.evaluate(t, u)
        u1 += dt / 3.0 * k

        u = u0 + dt / 2.0 * k
        t = t0 + dt / 2.0
        k = eqn.evaluate(t, u)
        u1 += dt / 3.0 * k

        u = u0 + dt * k
        t = t0 + dt
        k = eqn.evaluate(t, u)
        u1 += dt / 6.0 * k
        u1.time = t

        return (t, u1)


class RungeKutta4Dirichlet(ABCSolver):
    """An implementation of the Runge Kutta 4 routine that calls
    system.dirichlet_boundary."""

    def __init__(self, **kwds):
        super(RungeKutta4Dirichlet, self).__init__(**kwds)
        if not hasattr(self.system, "dirichlet_boundary"):
            raise Exception(
                "%s does not implement dirichlet_boundary method" % self.system
            )

    def advance(self, t0, u0, dt):
        """See the ABCSolver.advance method for documentation.

        Very simple minded implementation of the standard 4th order Runge-Kutta
        method to solve an ODE of the form
        fdot = rhs(t,f) that allows for the implementation of
        Dirichlet conditions during evaluation.

        Ensure that the corresponding system file has a method called,
        "dirichlet_boundary".
        """
        eqn = self.theEqn
        u = u0
        t = t0
        k = eqn.evaluate(t, u, intStep=1)
        u1 = u0 + (dt / 6.0) * k
        u1 = eqn.dirichlet_boundary(u1, intStep=1)

        u = u0 + dt / 2.0 * k
        t = t0 + dt / 2.0
        k = eqn.evaluate(t, u, intStep=2)
        u1 += dt / 3.0 * k
        u1 = eqn.dirichlet_boundary(u1, intStep=2)

        u = u0 + dt / 2.0 * k
        t = t0 + dt / 2.0
        k = eqn.evaluate(t, u, intStep=3)
        u1 += dt / 3.0 * k
        u1 = eqn.dirichlet_boundary(u1, intStep=3)

        u = u0 + dt * k
        t = t0 + dt
        k = eqn.evaluate(t, u, intStep=None)
        u1 += dt / 6.0 * k
        u1.time = t
        u1 = eqn.dirichlet_boundary(u1, intStep=None)

        return (t, u1)
