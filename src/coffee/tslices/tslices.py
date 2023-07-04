#!/usr/bin/env python
# encoding: utf-8
"""
A timeslice is the package which contains all data necessary to represent a
solution of the system at a given point in time.

This module contains the abstract base class (abc) for TimeSlice objects and a
default implementation that makes certain assumptions. It is highly recommended
that you subclass either of these classes.

As usual the abc is intended to make you think twice about doing things
differently. All other classes that interact with TimeSlice objects assume that
the interface is as given below.

Now it should be noted that the variables ABCTimeSlice.data,
ABCTimeSlice.domain and ABCTimeSlice.time should have been declared as
abc.abstractproperty. Doing this is possible, but because of a limitation in
the abc class for Python 2 a substantial amount of additional code is needed.
So this hasn't been done. In fact NONE of the methods or variables are declared
as abstract. 

What is recommended is that you subclass ABCTimeSlice and call it's functions
once appropriate checks have been carried out.

At some point this will need to be changed... I assume after the move to Python
3 is made.

TimeSlice objects should contain all information needed to interperet the
values of the functions being numerically evolved. In particular, aim to ensure
that you can restart a simulation if given a timeslice, system and solver.
"""

from builtins import object
import abc
from coffee.settings import be

###############################################################################
# TimeSlice Abstract Base Class
##############################################################################
class ABCTimeSlice(object, metaclass=abc.ABCMeta):
    """The abstract base class for TimeSlice objects.

    The interface below is assumed by all other classes that interact with
    TimeSlice objects (which is just about everything).

    """

    def __init__(self, data, domain, time, name=None, *args, **kwds):
        """Make a TimeSlice instance!

        It is recommended to call this method once your subclass has done
        appropriate vetting of data, domain and time. Not only because this
        ensures things are stored in the right places, but also because it sets
        up a logging object, logging.getLogger(self.name).

        Parameters
        ----------

        data :
            The values of the functions being solved for.
        domain :
            The grid over which the values are calcualted.
        time : float
            The time for which the values in data are valid.

        name : string
            The name of your subclass. This defaults to ABCTimeSlice and is
            used during logging.

        """
        self.data = data
        self.domain = domain
        self.time = time
        if name is None:
            self.name = "ABCTimeSlice"
        else:
            self.name = name
        super(ABCTimeSlice, self).__init__(*args, **kwds)

    def __repr__(self):
        """Returns a naive string representation of the slice.

        Returns
        -------
        string

        """
        s = "%s(data = %s, domain = %s, time = %s)" % (
            self.name,
            repr(self.data),
            repr(self.domain),
            repr(self.time),
        )
        return s

    def external_slices(self):
        """Returns a list of tuples that describe boundaries of grids for
        external boundaries.

        This is purely for convenience so that users don't need to
        write tslice.domain.external_slices().

        Returns
        -------
        list, slices
            The result of a call to self.domain.external_slices(self.data.shape)

        """
        return self.domain.external_slices(self.data.shape)

    def communicate(self, ghost_point_processor=None, data=None):
        """Returns a list of tuples that describe the boundaries of
        grids for internal boundaries.

        This is currently used when the full domain is divided into
        subdomains for use with mpi.
        This is purely for convenience so that users don't need to
        write tslice.domain.communicate().

        Parameters
        ----------
        ghost_point_processor : function
            A function with signature (data, list of slices). See the
            documentation for ABCGrid.communicate() for more information.

        data :
            The information to be communicated. It defaults to self.data.

        Returns
        -------
        list, slices
            The result of a call to self.domain.communicate(self.data)

        """
        return self.domain.communicate(
            data if data is not None else self.data,
            ghost_point_processor=ghost_point_processor,
        )

    def barrier(self):
        """A wrapper to self.domain.barrier()."""
        return self.domain.barrier()

    def collect_data(self):
        """A wrapper to self.domain.collect_data.

        Returns a timeslice containing the almalgum of all data and domains
        across all subdomains. It allocates new memory for the returned
        timeslice.
        This is currently used to pass non-distributed data to actions in the
        ibvp method.

        Returns
        -------
        tslice.TimeSlice
            A single timeslice which represents the complete data and
            domain for that point of time.

        """
        data_all = self.domain.collect_data(self.data)
        domain_all = self.domain.full_grid
        r_tslice = self.__class__(data_all, domain_all, self.time)
        return r_tslice

    def __add__(self, other):
        # It is very important that the rv is other + self.data not
        # self.data + other.
        # The reason (seems) is because as self.data can be an
        # nd.array the sum self.data + other
        # ends up being computed as the elements of self.other
        # plus other, itself an nd.array.
        # This screws up the array of elements.
        # Putting other first ensures that if other is a timeslice
        # then the sum isn't distributed over the elements of self.data.
        try:
            rv = other + self.data
        except:
            raise NotImplementedError(
                "Addition of %s and %s is not implemented" % (self, other)
            )
        if isinstance(rv, self.__class__):
            return rv
        else:
            return self.__class__(rv, self.domain, self.time, name=self.name)

    def __iadd__(self, other):
        if isinstance(other, self.__class__):
            self.data += other.data
        else:
            try:
                self.data += other
            except:
                raise NotImplementedError(
                    "Addition of %s and %s is not implemented" % (self, other)
                )
        return self

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        try:
            rv = self.data - other
        except:
            raise NotImplementedError(
                "Subtraction of %s and %s is not implemented" % (self, other)
            )
        return self.__class__(rv, self.domain, self.time)

    def __isub__(self, other):
        try:
            self.data -= other
        except:
            raise NotImplementedError(
                "In place subraction of %s and %s is not implemented" % (self, other)
            )
        return self

    def __rsub__(self, other):
        try:
            r_time_slice = other - self.data
        except:
            raise NotImplementedError(
                "Reflected subtraction of %s and %s is not implemented" % (other, self)
            )
        return r_time_slice

    def __mul__(self, other):
        try:
            rv = self.data * other
        except:
            raise NotImplementedError(
                "Multiplication of %s and %s is not implemented" % (self, other)
            )
        return self.__class__(rv, self.domain, self.time)

    def __imul__(self, other):
        try:
            self.data *= other
        except:
            raise NotImplementedError(
                "In place multiplicatio of %s and %s is not implemented" % (self, other)
            )
        return self

    def __rmul__(self, other):
        return self * other

    def __div__(self, other):
        try:
            rv = self.data / other
        except:
            raise NotImplementedError(
                "Division of %s and %s is not implemented" % (self, other)
            )
        return self.__class__(rv, self.domain, self.time)

    def __idiv__(self, other):
        try:
            self.data /= other
        except:
            raise NotImplementedError(
                "In place division of %s and %s is not implemented" % (self, other)
            )
        return self

    def __truediv__(self, other):
        return self.__div__(other)

    def __itruediv__(self, other):
        return self.idiv(other)

    def __pow__(self, other):
        try:
            rv = self.data**other
        except:
            raise NotImplementedError(
                "Exponentiation of %s by %s is not implemented" % (self, other)
            )
        return self.__class__(rv, self.domain, self.time)

    def __ipow__(self, other):
        try:
            self.data **= other
        except:
            raise NotImplementedError(
                "In place exponentiation of %s by %s is not implemented" % (self, other)
            )
        return self


###############################################################################
# Concrete implementations
###############################################################################
class TimeSlice(ABCTimeSlice):
    """A default subclass of ABCTimeSlice.

    This implementation assumes that all data are numpy arrays of the same
    shape. If your data is not like this then you should subclass ABCTimeSlice.
    This might be appropriate when working with spectral methods, for example.

    Due to the assumption that the data is represented as a numpy array, this
    class also implements a large number of methods which allow this object to
    be added, multiplied, etc...

    To be honest rather than implementing all the additional methods by hand
    it'd be easier just to make this default TimeSlice a subclass of be.ndarray
    itself and allow TimeSlice.data to access the underlying array. But this
    hasn't been done.
    """

    def __init__(self, data, *args, **kwds):
        """Instantiates a TimeSlice object.

        Calls the ABCTimeSlice __init__ method after handling the data.
        Thus you should refer to the documentation for ABCTimeSlice.__init__()
        for details of other parameters.

        Parameters
        ----------
        data :
            Converts data to a numpy array before passing it on to ABCTimeSlice.

        """
        data = be.array(data)
        if "name" not in kwds:
            kwds["name"] = "TimeSlice"
        super(TimeSlice, self).__init__(data, *args, **kwds)
