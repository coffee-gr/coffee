"""
Contains classes that abstract the MPI C API.

Classes in this module model the MPI node network and manage the relationship
of the complete domain to each 
"""
from builtins import range
from builtins import object
import abc
from coffee.settings import be


###############################################################################
# Abstract Base Class for MPI Interfaces
###############################################################################
class MPIInterface(object, metaclass=abc.ABCMeta):
    def __init__(self, mpi_comm, *args, **kwds):
        """Initialisaer for MPIInterface objects.

        Parameters
        ==========
        mpi_comm: mpi4py.Comm
            A wrapper to an MPI_COMM object via mpi4py.
        """
        self.comm = mpi_comm
        super(MPIInterface, self).__init__(*args, **kwds)

    @abc.abstractmethod
    def subdomain(self):
        """Returns a tuple of slice objects that represents the subdomain of the full
        grid that this MPI node should perform computations over.

        Returns
        =======
        tuple of slices
        """

    @abc.abstractmethod
    def communicate(self, data):
        """Returns a list of tuples. The first element of each tuple is a slice
        representing the portion of data which corresponds to the second
        element of the tuple. The second element of the tuple is an array of
        data which are the communicated values that correspond to the
        positions described by the slice."""

    def barrier(self):
        if self.comm is None:
            return
        return self.comm.barrier()


###############################################################################
# Concrete implementations
###############################################################################
class EvenCart(MPIInterface):
    def __init__(self, domain, boundary_data, *args, **kwds):
        """Initialiser for EvenCart objects

        An EvenCart object is an implementation of MPIInterface that
        assumes that the global topology uses the MPI Cartesian methods.

        Parameters
        ==========
        domain: numpy.ndarray
            A numpy array representing the values of points in the global
            computational domain
        boundary_data: grid.ABCBoundary
        """
        super(EvenCart, self).__init__(*args, **kwds)
        self.domain = domain
        self.domain_mapping = self._make_domain_mappings(domain, boundary_data)

    def _make_domain_mappings(self, domain, boundary_data):
        """A utility that constructs the tuples of slices that represent
        the sub-grids for each MPI node based on that nodes rank.

        The calculation assumes that each MPI node should handle roughly the
        same load as all others and that the MPI topology is Cartesian.

        Parameters
        ==========
        domain: numpy.ndarray
            A numpy array representing the values of points in the global
            computational domain
        boundary_data: grid.ABCBoundary

        Returns
        =======
        list of tuples of slices:
            Each tuple of slices in index i gives the sub grid for MPI node
            with rank i.

        """
        if self.comm is None:
            return tuple([slice(None, None, None) for dim in domain])
        r_map = []
        for rank in range(self.comm.size):
            coords = self.comm.Get_coords(rank)
            r_map += [
                tuple(
                    [
                        self._array_slice(
                            domain[i],
                            coord,
                            self.comm.dims[i],
                            boundary_data.ghost_points(i, -1),
                            boundary_data.ghost_points(i, 1),
                        )
                        for i, coord in enumerate(coords)
                    ]
                )
            ]
        return r_map

    def _neighbour_slices(self, shape, b_data):
        """A utility function that collects the source and destination
        ranks for neighbouring MPI nodes along with the slices needed
        to describe what data should be transferred.

        Parameters
        ==========
        shape: tuple of ints
            The shape of the data that will be communicated to neighbouring MPI
            nodes

        b_data: grid.ABCBoundary

        Returns
        =======
        list of tuples:
            Each tuple has three members, first the source MPI node rank,
            second the destination MPI node rank and third the list of tuples
            of slices that describes the data to be communicated.
        """
        if self.comm is None:
            return []
        dims = self.comm.Get_dim()
        pos_neighbours = [
            b_data.source_and_dest(d, 1) + b_data.internal_slice(shape, d, 1)
            for d in range(dims)
        ]
        neg_neighbours = [
            b_data.source_and_dest(d, -1) + b_data.internal_slice(shape, d, -1)
            for d in range(dims)
        ]
        neighbours = pos_neighbours + neg_neighbours
        return neighbours

    def _array_slice(
        self, array_length, rank, num_ranks, ghost_points_start, ghost_points_end
    ):
        """The utility method that calculates the slices that describe the
        sub-grid managed by the MPI node with rank ``rank``.

        The method assumes that it is doing the calculation for a single
        dimension. To build the full list of tuples of slices, all
        dimensions will need to be iterated over.

        Parameters
        ==========
        array_length: int
            The number of grid points in this particular dimension of the
            total computational domain.
        rank: int
            The MPI rank that the slices are computed for.
        num_ranks: int
            The number of MPI nodes that cover this particular dimension.
        ghost_points_start: int
            The number of ghost points at the start of the dimension.
        ghost_points_end: int
            The number of ghost points at the end of the dimension.

        Returns
        =======
        slice:
            A slice which gives the subdomain for the particular dimension and
            the needed MPI node rank.
        """
        # divide domain into appropriate parts
        q, r = divmod(array_length, num_ranks)
        # use the rank to details which part is relevant for this process
        s = rank * q + min(rank, r)
        e = s + q
        # Adjust e to account for min(rank, r) term, which spreads the remainder
        # over the appropriate number of processes
        if rank < r:
            e = e + 1
        # add in ghost_points if we can
        # this currently works for gp = 1 for an SAT boudnary method.
        # The code will need to be checked for consistency if this changes.
        if s - ghost_points_start > -1:
            s = s - ghost_points_start
        if e + ghost_points_end < array_length:
            e = e + ghost_points_end
        return slice(s, e, None)

    @property
    def subdomain(self):
        """Returns the slices which give the subdomain for the MPI node
        that this method is called on.

        Returns
        =======
        list of tuples of slices
        """
        if self.comm is None:
            return self.domain_mapping
        return tuple(self.domain_mapping[self.comm.rank])

    def communicate(self, data, b_data):
        """Communicates the given data to appropriate neighbouring nodes.

        Parameters
        ==========
        data: numpy.ndarray
            An array representing the values of the computed function over
            the sub-grid managed by the MPI node this code is being computed on.
        b_data: grid.ABCBoundary

        Returns
        =======
        list of tuples (list of slices, numpy.ndarray)
            The returned data is a list of tuples. The first element of the tuple
            describes which grid points the communicated data relates to.
            The second element is the data received from whichever MPI node
            manages the sub-grid which looks after the relevant grid points.
        """
        # if self.comm.size == 1:
        # return []
        nslices = self._neighbour_slices(data.shape, b_data)
        r_data = []
        for source, dest, send_slice, recv_slice in nslices:
            if dest < 0:
                send_data = None
            else:
                # I really did not want to copy this view. For large data
                # sets this feels like an unnecessary delay.
                # However mpi4py did not seem to work with views and using
                # self.comm.sendrecv to send and recv arbitrary objects.
                # Plus self.comm.Sendrecv requires contiguous data
                # so only in rare cases is it possible to send.
                # Hence the be.copy statement is required.
                send_data = be.array(data[send_slice], copy=True, order="C")
            if source < 0:
                recv_data = None
            else:
                # when testing the above assertion I suggest changing this
                # line to be.ones_like, rather than be.empty_like
                recv_data = be.empty_like(data[recv_slice])
            self.comm.Sendrecv(send_data, dest=dest, recvbuf=recv_data, source=source)
            if source >= 0:
                r_data += [(recv_slice, recv_data)]
        return r_data

    def collect_data(self, data):
        """A method that collects all computed function values over all
        sub-grids managed by all MPI nodes and returns the amalgamated data
        if the rank of this process is 0.

        This gives a way to get a single copy of all computed data so far. It does,
        however, assume that all subsequence processing occurs on a single MPI
        node. This is a bottle neck.

        Parameters
        ==========
        data: numpy.ndarray
            The data to be communicated.

        Returns
        =======
        None or numpy.ndarray
            The collated data if MPI_Comm_Rank is 0, none otherwise.
        """

        # Note that this method does not take account of ghost_points in
        # the domains. This does not cause a problem. It just means
        # that more data than necessary is passed.

        # If there are no partitions of the domain we don't need to collect
        if self.comm is None:
            return data
        # Gather the data
        fields = self.comm.gather(data, root=0)
        # If this process is root then collate the data and return
        if self.comm.rank == 0:
            rdata_edims_shape = (data.shape[0],) + data.shape[len(self.domain) + 1 :]
            rdata_edims_slice = tuple(
                [slice(None, None, None) for i in rdata_edims_shape]
            )
            rdata = be.zeros(
                (rdata_edims_shape[0],) + self.domain + rdata_edims_shape[1:],
                dtype=data.dtype,
            )
            for rank, field in enumerate(fields):
                dslice = self.domain_mapping[rank]
                rdata[(rdata_edims_slice[0],) + dslice + rdata_edims_slice[1:]] = field
            return rdata
        else:
            return None
