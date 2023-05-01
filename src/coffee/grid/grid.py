#!/usr/bin/env python
# encoding: utf-8

"""
Grid objects represent the computational, spatial, domain.

Each grid object understands it's own discretisation as well as how
"""


from builtins import range
from builtins import object
import math
import numpy as np
import abc
import logging

from coffee.mpi import mpiinterfaces
from future.utils import with_metaclass
from coffee.backend import be

################################################################################
# Base Boundary data class
################################################################################
class ABCBoundary(with_metaclass(abc.ABCMeta, object)):
    """The abstract base class for boundary classes.

    A boundary class manages data associated to ghost points and the slices
    specifying which data points should be communicated between MPI nodes.

    Suppose that our computational domain is a line, represented as the
    array::

        domain = be.linspace(0,1,num=11)

    This grid has two external boundaries at indices 0 and -1. 
    A subclass of ABCBoundary will generate a list of slices::

        [slice(0,None,None)

    that when applied to ``domain`` extract the values of the boundary points.

    Suppose that our computational domain is a line now represented as two
    sub-domains::

        domain_1 = be.linspace(0,1,num=11)[:6]
        domain_2 = be.linspace(0,1,num=11)[5:]

    Each array now has an external boundary, which represents a boundary of
    the global computational domain, and an internal boundary, which
    represents an artifical boundary generated by the representation of the
    global domain as two different arrays.

    A subclass of ABCBoundary should draw a distinction between internal
    boundaries and external boundaries. This is so that system objects
    can handle each type of boundary differently. The design goal
    is for the lists to be iterated over in the system object so that
    each boundary of each type can be appropriately treated.
    """

    # I use -1 and 1 to reduce friction with the mpi direction variable
    LEFT = -1
    RIGHT = 1

    DIRECTION_ERROR = ValueError("Direction must be +/- 1")

    @staticmethod
    def _direction_to_index(direction):
        """A utility that translates between the MPI direction variable and
        the index of a hypothetical
        numpy array that represents the data over a grid.

        This method assumes that the grid points increase in value with
        increasing index. Hence a direction of -1 (meaning left) returns 0
        the first element in the array and a direction of 1 (meaning right)
        returns -1 meaning the last element of the array.

        The method is used to build the list of slices that extract the
        relavant sub-array of data on the 'boundary' of the grid.

        Parameters
        ==========
        direction: int
            Can be either -1 or 1 in reference to the MPI direction.

        Returns
        =======
        int
            The index of a hypothetical array that contains the data in the
            'direction' end of the array. So, if the grid points inc

        Raises
        ======
        ValueError
            If direction is not -1 or 1.

        """
        if direction == -1:
            return 0
        elif direction == 1:
            return 1
        raise DIRECTION_ERROR

    @staticmethod
    def _empty_slice(dimensions):
        """A utilty to construct a list of empty slices with length dimensions.

        Parameters
        ==========
        dimensions: int
            The number of dimensions of the spatial domain of the function to
            be computed. Also: the number of dimensions in the grid.

        Returns
        =======
        list, slices
            Returns a list of empty slices of length dimension
        """
        return [slice(None, None, None) for i in range(dimensions)]

    def __init__(self, *args, **kwargs):
        pass

    @abc.abstractmethod
    def ghost_points(self, dimension, direction):
        """Returns the number of ghost points for the specified dimension and
        direction.

        Parameters
        ==========
        dimension: int
            Which dimension of the domain are we asking for the number of 
            ghost points?

        direction: int
            Of the one dimensional array representing the given dimension, which
            end of the array are we interested in knowing about the ghost points?

        Returns
        =======
        int
            The number of ghost points for the specified boundary.
        """

        return 0

    @abc.abstractmethod
    def internal_slice(self, shape, dimension, direction):
        """Returns a tuple of slices which, when applied
        to a data array, gives the data to be communicated to 
        neighbouring grids.

        Parameters
        ==========
        shape: tuple of ints
            The shape of the array to be communicated. For example the shape
            of a numpy array representing the computational domain.

        dimension: int
            Which dimension do we need the list of slices for?

        direction: int
            Which direction do we need the slice for?

        Returns
        =======
        tuple of slices
            A list containing the slices that extract the data on the boundary.
        """
        return self._empty_slice(shape)

    @abc.abstractmethod
    def external_slice(self, shape, dimension, direction):
        """
        Returns a tuple which, when applied
        to a data array, gives the data on the extenal boundaries of the grid.
        
        Parameters
        ==========
        shape: tuple of ints
            The shape of the array to be communicated. For example the shape
            of a numpy array representing the computational domain.

        dimension: int
            Which dimension do we need the list of slices for?

        direction: int
            Which direction do we need the slice for?

        Returns
        =======
        list of slices
            A list containing the slices that extract the data on the boundary.
        """
        slices = self._empty_slice(len(shape))

        # +1 because for some reason I decided that vector data should be in the
        # first dimension.
        dim_index = dimension + 1
        if direction == -1:
            slices[dim_index] = slice(None, 1, None)
            return slices

        if direction == 1:
            slices[dim_index] = slice(-1, None, None)
            return slices

        raise DIRECTION_ERROR

    def external_slices(self, shape):
        """Returns a list of tuples. Each tuple contains an integer,
        representing the dimension, a direction and the result of 
        calling external_slice(shape, dimension, direction).

        This is a convenience method to make iteration over external
        boundaries easy for the user.
        
        Parameters
        ==========
        shape: tuple of ints
            The shape of the computational domain. For example if the computational
            domain is a numpy array then the shape is the numpy shape of the array.

        Returns
        =======
        list of tuples of slices
            Each tuple of slices represents the sub array of the computational
            domain that made up of external boundary points.
        """

        # This implicitly assumes that the number of dimensions handled by
        # the grid object is the same as than managed by the mpi_comm
        dims = self.mpi_comm.Get_dim()
        if __debug__:
            self.log.debug("Shape = " + repr(shape))
        neg_slices = [
            (i, -1, self.external_slice(shape, i, -1)) 
            for i in range(dims) 
            if self.external_slice(shape, i , -1) != self._empty_slice(len(shape))
        ]
        pos_slices = [
            (i, 1, self.external_slice(shape, i, 1)) 
            for i in range(dims)
            if self.external_slice(shape, i , 1) != self._empty_slice(len(shape))
        ]
        if __debug__:
            self.log.debug("external slices are = " + repr(pos_slices + neg_slices))
        return pos_slices + neg_slices
 

################################################################################
# Base Boundary data class
################################################################################
class SingleCartesianGridBoundary(ABCBoundary):
    """Implements ABCBoundary for a simulation on a single grid with no 
    mpi dependence.

    That is a grid with no internal boundaries and every grid "edge" an 
    external boundary. There are no ghost points."""
    pass

class MPIBoundary(ABCBoundary):
    """MPIBoundary implements the ABCBoundary class, it assumes use of an underlying
    MPI.Cartcomm object wrapped by mpi4py.
    
    The numbers of ghost points
    and "internal points" can be specified directly using arrays of two tuples.
    If the grid is two dimensional the, for example, the array of ghost points
    and internal points will look like::
        
        [(a, b), (c, d)]

    The tuple ``(a, b)`` gives the ghost / internal points for the first dimension.
    The int ``a`` gives the number of ghost / internal points for the
    left direction and ``b`` gives the number for the right direction.
    Similarly for ``(c,d)``.
    
    
    """
    def __init__(
            self, 
            ghost_points, 
            internal_points, 
            mpi_comm=None, 
            *args, 
            **kwargs
        ):
        """ 
        The initialiser for MPIBoundary.

        The ghost points and internal points variables must both be an array 
        or tuple with one entry for each dimension. Each entry is a 2-tuple giving 
        the number of points in the negative and positive directions for 
        that axis. See class documentation for more information.

        Parameters
        ==========
        ghost_points: list of tuples of ints
            The number of ghost points in each direction for each dimension
        internal_points: list of tuples of ints
            The number of internal points in each direction for each dimension
        mpi_comm: MPI.Comm
            The mpi4py wrapper of MPI_COMM.

        """
        super(MPIBoundary, self).__init__(*args, **kwargs)
        self._ghost_points = ghost_points
        self._internal_points = internal_points
        self.mpi_comm = mpi_comm
        self.log = logging.getLogger("MPIBoundary")

    def ghost_points(self, dimension, direction):
        "See ABCBoundary.ghost_points() for documentation."""
        return self._ghost_points[dimension][self._direction_to_index(direction)]

    def source_and_dest(self, dimension, direction):
        """A wrapper for the MPI_Cart_Shift function via mpi4py.

        To be used in conjuction with a wrapper for MPI_Sendrecv.

        Parameters
        ==========
        dimension: int
            The dimension of the shift
        direction: int
            The direction of the shift

        Returns
        =======
        tuple of ints
            A tuple of ints, ``(source, destination)``, where
            ``source`` is the rank of the mpi node that is the
            source of the data to be sent
            and ``destination`` is the rank of the mpi node that is to
            receive the data. 
        """
        return self.mpi_comm.Shift(dimension, direction)

    def internal_slice(self, shape, dimension, direction):
        """See ABCBoundart.internal_slice() for documentation.
        
        The first tuple in the return of the function is the list of
        slices that identifies the data to be sent. The second tuple
        is the list of slices that identifies the data to be received.

        Note that the return variable for this method is inconsistent
        with the ABCBoundary method. The ABCBoundary method should be 
        updated.

        Returns
        =======
        A 2-tuple of lists of tuples of ints
        """
            
        source, dest = self.source_and_dest(dimension, direction)

        if dest < 0:
            rsend_slice = None
        else:
            send_slice = self._empty_slice(len(shape))
        if source < 0:
            rrecv_slice = None
        else:
            recv_slice = self._empty_slice(len(shape))

        i_point = self._internal_points[dimension][
            self._direction_to_index(direction)
        ]
        g_point = self._ghost_points[dimension][
            self._direction_to_index(direction)
        ]
        total_g_points = sum(self._ghost_points[dimension])

        # +1 because the first dimension of the data array is reserved for vector
        # valued data. This may not have been a good desicion
        #
        # The details of why the start and end points for the slices are as they
        # are is complicated. Double check things before you change anything.
        dim_index = dimension + 1

        # This functions are used to to avoid difficulties resulting from the
        # difference between slice(-1, 0) and slice(-1, None), we don't need one
        # for "pos_or_none" as slice(0, x) and slice(None, x) behave the same
        def neg_or_none(number):
            return number if number < 0 else None

        if direction == 1:
            if not dest < 0:
                send_slice[dim_index] = slice(
                    -total_g_points, 
                    neg_or_none(-total_g_points + i_point), 
                    None
                )
                rsend_slice = tuple(send_slice)
            if not source < 0:
                recv_slice[dim_index] = slice(
                    None,
                    i_point,
                    None
                )
                rrecv_slice = tuple(recv_slice)
        elif direction == -1:
            if not dest < 0:
                send_slice[dim_index] = slice(
                    max(0, total_g_points - i_point),
                    total_g_points, 
                    None
                )
                rsend_slice = tuple(send_slice)
            if not source < 0:
                recv_slice[dim_index] = slice(
                    -i_point,
                    None,
                    None
                )
                rrecv_slice = tuple(recv_slice)
        return rsend_slice, rrecv_slice

    def external_slice(self, shape, dimension, direction):
        """See ABCBoundary.external_slice() for documentation."""
        if self.mpi_comm.periods[dimension]:
            return self._empty_slice(len(shape))
        coords = self.mpi_comm.Get_coords(self.mpi_comm.rank)
        dim = self.mpi_comm.dims[dimension]
        if coords[dimension] == 0 and direction == -1:
            return super(MPIBoundary, self).external_slice(
                shape, 
                dimension, 
                -1
            )
        if coords[dimension] == dim - 1 and direction == 1:
            return super(MPIBoundary, self).external_slice(
                shape,
                dimension, 
                1
            )
        return self._empty_slice(len(shape))

class SimpleMPIBoundary(MPIBoundary):
    """SimpleMPIBoundary is for boundaries in mpi contexts with a fixed number of
    ghost points in every direction and dimension and the "internal points" equal
    to the ghost points."""

    def __init__(self, ghost_points, *args, **kwargs):
        """The initialiser for SimpleMPIBoundary.

        This class wraps MPIBoundary. It is essentially a short cut for the
        common situation where the number of ghost points is the same for 
        every direction and dimension and the number of internal points is
        the same as the number of ghost points.

        Parameters
        ==========
        ghost_points: int
            The number of ghost points.

        """
        number_of_dimensions = kwargs["number_of_dimensions"]
        points_tuple = [
            (ghost_points, ghost_points) for d in range(number_of_dimensions)
        ]
        super(SimpleMPIBoundary, self).__init__(
            points_tuple,
            points_tuple,
            *args,
            **kwargs
        )
    
class GeneralBoundary(ABCBoundary):
    """GeneralBoundary implements ABCBoundary for grids in an mpi context with arbitrary
    ghost points and arbitrary interal and external slices fixed at run time.
    
    The ghost points variable must be an array or tuple with one entry for each
    dimension. Each entry is a 2-tuple giving the number of ghost points in the
    negative and positive directions for that axis.

    The data slices variable is an array with one entry for each dimension. Each
    entry is a tuple of slices which when applied to the data stored in a tslice
    gives the data to be communicated / received.
    """

    def __init__(self, ghost_points, internal_slices, external_slices):
        """Inialiser for GeneralBoundary.

        Parameters
        ==========
        ghost_points : list of two tuples of ints
            Each entry of the list gives the number of ghost points for each direction
            for each dimension.
        internal_slices: list of tuples of slices
            Slices which identify the internal boundaries.
        external_slices: list of tuples of slices
            Slices which identify the external boundaries.
        """
        self._ghost_points = ghost_points
        self._interal_slices = internal_slices
        self._external_slices = external_slices

    def ghost_points(self, dimension, direction):
        """See ABCBoundary.ghost_points() for documentation."""
        return self._ghost_points[dimension][self._direction_to_index(direction)]

    def external_slices(self, shape, dimension, direction):
        """See ABCBoundary.external_slices() for documentation."""
        return self._external_slices[dimension][self._direction_to_index(direction)]


################################################################################
# Base Grid class
################################################################################
class ABCGrid(with_metaclass(abc.ABCMeta, object)):
    """The abstract base class for Grid objects.

    A grid class understands both the total computational domain as well as the
    domain associated to the particular MPI node that computation is performed 
    on. Each grid object wraps a boundary object and abstracts the boundary
    object api.
    
    """
    
    def __init__(self, 
            shape, bounds, name = "Grid", comparison = None,
            mpi = None, boundary_data = None, *args, **kwds
        ):
        """The initialiser for grid objects.

        Parameters
        ==========
        shape: tuple of ints, 
            The shape of a theoretical numpy array that represents the computational
            domain.
        bounds: list of tuples of floats
            For each dimension the tuple gives the maximum and minimum extent of the
            global computational domain
        name: string, Optional
            A name for this grid. Used when logging.
        comparison: Optional
            A argument used when performing comparison of different simulations
            on different grids. Stored in the hdf file during computation. It
            allows for easy identification of how different computation grids
            should be compared. 
        mpi: mpi.MPIInterface, Optional
            An instance of an mpi4py wrapped MPI_COMM object.
        boundary_data: ABCBoundary, Optional
            An implementation of ABCBoundary
        """
        self.mpi = mpi
        self.dim = len(shape)
        self.name = name
        self.log = logging.getLogger(name)
        self.comparison = comparison
        self.shape = shape
        self.bounds = bounds
        self.boundary_data = boundary_data
    
    def __strs__(self):
        return self.name

    def __repr__(self):
        return "<%s shape=%s, bounds=%s, comparison=%s, mpi=%s>"%(
            self.name, self.shape, self.bounds, self.comparison, self.mpi
            )

    def communicate(self, data, ghost_point_processor=None):
        """Wraps the mpi.mpiinterface.communicate() method.

        The ghost_point_processor is an arbitary function that
        takes the data and the returned information of the
        mpi.mpiinterfaces.communicate() method, called ``b_values`` here. 
        It must modify
        these data structure in place. It is called directly
        before the now modified data and ``b_values`` is returned.

        Parameters
        ==========
        data: 
            The data to be communicated
        ghost_point_processor: function(data, grid.ABCBoundary)

        Returns
        =======
        tuple:
            The first element of the tuple is the data that was passed
            to this method. The second element is the returned information
            to the mpi.mpiinterfaces.communicate() method.
        """
        if self.mpi is None:
            return data
        b_values = self.mpi.communicate(data, self.boundary_data)
        if ghost_point_processor:
            ghost_point_processor(data, b_values)
        return data, b_values

    def barrier(self):
        """Provides an easy way to call the MPI_Barrier() method via
        the mpi.mpiinterface object."""
        if self.mpi is None:
            return 
        return self.mpi.barrier()

    def external_slices(self, data_shape):
        """Wraps ABCBoundary.external_slices.

        See ABCBoundary.external_slices() for documentation.
        """
        if __debug__:
            self.log.debug("In grid.external_slices")
        return self.boundary_data.external_slices(data_shape)

    def collect_data(self, data):
        """Wraps mpi.mpiinterfaces.collect_data().

        See mpi.mpiinterfaces.collect_data() for documentation."""
        if self.mpi is None:
            return data
        return self.mpi.collect_data(data)
    
    @abc.abstractproperty
    def axes(self): pass
    
    @abc.abstractproperty
    def step_sizes(self): pass
    
    @property
    def meshes(self):
        """Generates numpy arrays of grid points.

        The method uses be.lib.stride_tricks.as_strided to cut down on computational
        time. It relies on information in the object itself to handle shape,
        number of points and values.

        The meshes are useful for vectorised functions when computing things like
        initial data for the system.

        Returns
        =======
        list of numpy.ndarray:
            A list of the coordinate values by dimension.
        """
        axes = self.axes
        grid_shape = tuple([axis.size for axis in axes])
        mesh = []
        for i, axis in enumerate(axes):
            strides = be.zeros((len(self.shape), ), dtype=int)
            strides[i] = axis.itemsize
            mesh += [be.lib_stride_tricks_as_strided(
                axis,
                grid_shape,
                strides,
                writeable=False
            )]
        return mesh

################################################################################
# Constructors for specific cases
################################################################################
class UniformCart(ABCGrid):
    """An implementation of ABCGrid that assumes that the grid has uniform
    steps and is Cartesian."""
    
    def __init__(self, 
            shape, 
            bounds, 
            mpi_comm=None, 
            comparison=None, 
            name=None, 
            *args, **kwds
        ):
        """Initialiser for UniformCart.

        Parameters
        ==========
        shape: tuple of ints, 
            The shape of a theoretical numpy array that represents the computational
            domain.
        bounds: list of tuples of floats
            For each dimension the tuple gives the maximum and minimum extent of the
            global computational domain
        mpi_comm: mpi.MPIInterface, Optional
            An instance of an mpi4py wrapped MPI_COMM object.
        comparison: Optional
            A argument used when performing comparison of different simulations
            on different grids. Stored in the hdf file during computation. It
            allows for easy identification of how different computation grids
            should be compared. 
        name: string
            The name of this object. Used in output to hdf files and logging.
        """

        _shape = tuple([s+1 for s in shape])
        if mpi_comm is None:
            mpi = None
        else:
            mpi = mpiinterfaces.EvenCart(
                _shape, 
                kwds.get("boundary_data", None),
                mpi_comm=mpi_comm, 
            )
        if name is None:
            name = "UniformCart%s%s%s"%(shape, bounds, comparison)
        super(UniformCart, self).__init__(
            shape, bounds, 
            name=name, comparison=comparison,
            mpi=mpi, *args, **kwds
            ) 
        _axes = [
            be.linspace(
                self.bounds[i][0], self.bounds[i][1], self.shape[i]+1
            )
            for i in range(len(self.bounds))
            ]
        self._step_sizes = [axis[1]-axis[0] for axis in _axes]
        if self.mpi is None:
            self._axes = _axes
        else:
            self._axes = [
                axis[self.mpi.subdomain[i]] 
                for i, axis in enumerate(_axes)
                ]

    @property
    def axes(self):
        """Returns the coordinate values over each dimension.

        Returns
        =======
        numpy.ndarray
        """
        return self._axes

    @property
    def full_grid(self):
        """Returns a UniformCart object that represents the total
        computational domain, rather than just the domain used by this MPI
        node.

        Returns
        =======
        UniformCart
        """
        if self.mpi is None:
            return self
        return UniformCart(
            self.shape, self.bounds, 
            comparison=self.comparison, name=self.name
            )
        
    @property
    def step_sizes(self):
        """Returns a numpy array containing the fixed step size in each dimension.

        Returns
        =======
        numpy.ndarray
        """
        return self._step_sizes

class GeneralGrid(ABCGrid):
    """An implementation of ABCGrid that allows for periodic dimensions.

    It currently does not allow for use in conjunction with MPI. Note that
    UniformCart and a period MPI network topology are compatible, however.
    """

    def __init__(self, 
            shape, 
            bounds, 
            periods,
            comparison=None, 
            name=None, 
            *args, 
            **kwds
        ):
        """The initialiser for GeneralGrid.

        Parameters
        ==========
        shape: tuple of ints, 
            The shape of a theoretical numpy array that represents the computational
            domain.
        bounds: list of tuples of floats
            For each dimension the tuple gives the maximum and minimum extent of the
            global computational domain
        periods: list of ints
            A non-zero value is assumed to indicated that the dimension at the
            given index is periodic. This is modelled on MPI periods.
        comparison: Optional
            A argument used when performing comparison of different simulations
            on different grids. Stored in the hdf file during computation. It
            allows for easy identification of how different computation grids
            should be compared. 
        name: string
            The name of this object. Used in output to hdf files and logging.
        """
        _shape = []
        for i,p in enumerate(periods):
            if p:
                _shape.append(shape[i])
            else:
                shape.append(shape[i]+1)
        _shape = tuple(_shape)
        mpi=None
        if name is None:
            name = "<GeneralGrid shape=%s, bounds=%s,periods=%s, comparison=%s>"%(
                shape, 
                bounds, 
                periods,
                comparison)
        super(GeneralGrid, self).__init__(
            shape, bounds, 
            name=name, 
            comparison=comparison,
            mpi=mpi, 
            *args, 
            **kwds
            ) 
        _axes = [
            be.linspace(
                self.bounds[i][0], self.bounds[i][1], self.shape[i]+1
            )
            for i in range(len(self.bounds))
            ]
        self._step_sizes = [axis[1]-axis[0] for axis in _axes]
        #self._axes = [axis[self.mpi.subdomain] for axis in _axes]

    @property
    def axes(self):
        """Returns the coordinate values over each dimension.

        Returns
        =======
        numpy.ndarray
        """
        return self._axes

    @property
    def full_grid(self):
        """Returns a GeneralGrid object that represents the total
        computational domain, rather than just the domain used by this MPI
        node.

        Of course, since GeneralGrid's don't use MPI, this method returns self.
        It is included for compatibility.

        Returns
        =======
        GeneralGrid
        """
        return self
        
    @property
    def step_sizes(self):
        """Returns a numpy array containing the fixed step size in each dimension.

        Returns
        =======
        numpy.ndarray
        """
        return self._step_sizes
