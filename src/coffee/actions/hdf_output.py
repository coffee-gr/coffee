"""This module writes the data of a timeslice to an hdf file.

The main class SimOutput uses objects of class SimOutput.SimOutputType to
actually write the data. By subclassing SimOutput.SimOutputType in the manner
specified below non - numpy array data can be stored conviently.

I think of this method as a plugin for the system. The user writes a new output
method for some kind of new data type in a sub class of SimOutputType and
'plugs' this method into the 'action' list of SimOutput.
"""
from builtins import str
from builtins import object
from coffee.settings import be

from .actions import Prototype
from ..io.simulation_data import dgTypes, systemD, DataGroup, sysDTypes
from functools import reduce


# A utility class to write to SimulationHDF via the
# actions.py framework.
class SimOutput(Prototype):
    """An action to handle output of data to hdf.

    Each piece of data to be output is passed, as an object, to this class in
    the array actionTypes.

    These objects are subclasses of SimOutputType a class which is accessible
    as an attribute of SimOutput.

    Note that there is a dependence on io.simulation data. In particular
    on the dgTypes dictionary which provides information on the appropriate
    names for created datagroups.
    For information about how to access this data please see the simulation_data
    module.

    Details on how to subclass SimOutputType are given in the documentation for
    the SimOutputType class.
    """

    def __init__(
        self,
        hdf_file,
        solver,
        theSystem,
        theInterval,
        actionTypes,
        cmp_=None,
        overwrite=True,
        name=None,
        **kwds,
    ):
        """The constructor for the SimOutput action.

        The goal for this object is to allow for complete replication of the
        simulation by storing the needed information in the hdf file. This
        requires that the system and solver being used have sensible string
        representations that allow for rebuilding exact duplicates of the
        corresponding objects used in the original simulation.

        Needless to say that this is not currently the case, but I think its a
        good idea.

        Parameters
        ==========
        hdf_file : h5py.File
            This is the file to which data will be written.
        solver : solvers.Solver
        theSystem : systems.System
        theInterval : grid.ABCGrid
        actiontypes : list of SimOutputType objects
                      Each time the SimOutput
                      action is run this list is iterated through and each
                      SimOutputType object is called. This object writes out
                      data according to it's own routine and stores it in the
                      given hdf file.
        cmp_ :
               Something that allows for comparison of simulations. In one
               dimensional simulations this is the number of grid points. This
               is used by the error script to work out which simulations should
               be compared.

        overwrite : bool
                    If you are outputting a simulation into an existing hdf
                    file do you want to allow the simulation to overwrite the
                    existing data?
        name : string
            What is the name of this simulation?
        """
        super(SimOutput, self).__init__(**kwds)
        if name == None:
            hour = str(time.localtime()[3])
            minute = str(time.localtime()[4])
            second = str(time.localtime()[5])
            name = hour + ":" + minute + ":" + second
        self.name = name
        self.hdf = hdf_file
        self.solver = solver
        self.system = theSystem
        self.domain = theInterval
        self.actions = actionTypes
        self.overwrite = overwrite
        self.cmp_ = cmp_
        for action in actionTypes:
            action.setup(self)

    def _doit(self, it, u):
        pass

    class SimOutputType(object):
        """If you want to customise the output to the hdf file subclass this
        class.

        This class handles the organisation of the hdf file structure so that
        you don't have to. It slices, it dices, it picks the group to store the
        datasets that will contain your data. That data group is stored in the
        atribute self.data_group.

        It also handles the inclusion of derived attribute in the .attrs
        variable of the relevant data_set. Sometimes this is needed if more
        than a numpy array needs to be retreived.

        """

        def __init__(self, derivedAttrs=None):
            """Initialiser for SimOutputType

            Parameters
            ==========
            derivedAttrs : A dictionary of keywords and functions, Optional
               The functions must have signature function(iteration,
               tslice, system) and must return an object that h5py
               knows how to store in an hdf file. That object is
               stored in self.data_group[it].attrs[key].
            """
            if derivedAttrs is None:
                self.derivedAttrs = {}
            else:
                self.derivedAttrs = derivedAttrs

        def setup(self, parent):
            """This method is called by SimOuput on construction.

            Overide this is additional setup as needed, in particular if data
            should be extracted from parent. See SimOutputType.System
            for an example.

            This method should always to called so if you do override it, make
            sure to include a call to super(YourClass, self).setup(partent).
            This is necessary because without the call the data_group may not
            be configured correctly.

            Parameters
            ==========
            SimOuput:
                parent the SimOutput class that this SimOutputType belongs to.
            """
            if parent.overwrite:
                self.data_group = DataGroup(
                    parent.hdf.require_group(self.groupname).require_group(parent.name)
                )
            else:
                self.data_group = DataGroup(
                    parent.hdf.create_group(self.groupname).create_group(parent.name)
                )
            self.parent = parent
            self.data_group.attrs["cmp"] = parent.cmp_

        def __call__(self, it, u):
            for key, value in list(self.derivedAttrs.items()):
                v = value(it, u, self.parent.system)
                self.data_group[it].attrs[key] = v

    class Data(SimOutputType):
        """Writes out tslices.TimeSlice.data.

        This class assumes that tslice.TimeSlice.data has a type that h5py
        understands.
        """

        groupname = dgTypes["raw"]

        def __call__(self, it, u):
            dg = self.data_group
            dg[it] = u.data
            super(SimOutput.Data, self).__call__(it, u)

    class Exact(SimOutputType):
        """Calls system.exact_value to write out the exact value of the system.

        This is useful if, during error calculation, you want the error rates
        against an exact solution.

        """

        groupname = dgTypes["exact"]

        def __call__(self, it, u):
            dg = self.data_group
            parent = self.parent
            dg[it] = parent.system.exact_value(u.time, u.domain).data
            super(SimOutput.Exact, self).__call__(it, u)

    class Times(SimOutputType):
        """Would you like to know what time the data at a particular data_set
        is meant to correspond to?

        Then you need this SimOutputType.

        """

        groupname = dgTypes["time"]

        def __call__(self, it, u):
            dg = self.data_group
            dg[it] = be.array([u.time])
            super(SimOutput.Times, self).__call__(it, u)

    class TimeStep(SimOutputType):
        """Just like time buts tells you what the calculated dt was.

        This can be calculated from the data in Times, only if you are
        recording the data in the next timestep...

        """

        groupname = dgTypes["dt"]

        def __call__(self, it, u):
            dg = self.data_group
            dg[it] = be.array([self.parent.system.timestep(u)])
            super(SimOutput.TimeStep, self).__call__(it, u)

    class Domains(SimOutputType):
        """Records domains of the simulation.

        This SimOutputType does not record the grid object, but rather the axes
        and comparison variables of the grid. Hence no mpi information is
        collected here.

        Currently this is because ibvp.IBVP before it runs actions collates all
        the data together so that, from the point of view of the hdf file there
        is only one process accessing it.

        This will need to be changed at some point. h5py has routines to allow
        more than one process to write to a file. But coffee does not
        currently capitalise on these routines.

        """

        groupname = dgTypes["domain"]

        def __call__(self, it, u):
            dg = self.data_group
            axes = u.domain.axes
            axes_shape = tuple([axis.size for axis in axes])
            axes_flat = be.empty(
                (reduce(lambda x, y: x + y, axes_shape),), dtype=u.domain.axes[0].dtype
            )
            start = 0
            for i, axis in enumerate(axes):
                axes_flat[start : start + axes_shape[i]] = axis
                start = start + axes_shape[i]
            dg[it] = axes_flat
            dg[it].attrs["axes_shape"] = axes_shape
            dg[it].attrs["shape"] = u.domain.shape
            dg[it].attrs["bounds"] = u.domain.bounds
            dg[it].attrs["comparison"] = u.domain.comparison
            super(SimOutput.Domains, self).__call__(it, u)

    class DerivedData(SimOutputType):
        """Runs the data through a user defined function before writing out to
        the data_set.

        The function must have the signature function(int,
        tslice.TimeSlice, system.System) and must return
        data of a type that h5py recognises.

        Note that there is a size limitation on the .attrs variable and
        therefore this SimOutputType maybe more appropriate.
        """

        def __init__(self, name, function, frequency=1, start=0, derivedAttrs=None):
            """The initialiser for DerivedData.

            Parameters
            ==========
            name : string
                The name for the hdf.datagroup.
            function : function(int, tslice.TimeSlice, system)
                The function that returns the data to be stored in
                the datagroup.
            frequency : int
                It dicates how many iterations pass for each call of this action.
            start : float
                The time from which to start performing the action.
            derivedAttrs :
                See the documentation for SimOutputType.
            """
            self.func = function
            self.groupname = name
            self.freq = frequency
            self.start = start
            super(SimOutput.DerivedData, self).__init__(derivedAttrs)

        def __call__(self, it, u):
            if it % self.freq == 0 and self.start < u.time:
                dg = self.data_group
                dg[it] = self.func(it, u, self.parent.system)
                super(SimOutput.DerivedData, self).__call__(it, u)

    class System(SimOutputType):
        """Attempts to write out enough information about the system to allow
        for exact reconstruction of the simulation that produced the data being
        stored.

        Of course most things arn't that simple... and this SimOutput type
        requires care during use.
        """

        groupname = systemD

        def setup(self, parent):
            super(SimOutput.System, self).setup(parent)
            g = self.data_group.group
            psystem = be.asarray(repr(parent.system)).astype("S")
            pgrid = be.asarray(repr(parent.domain)).astype("S")
            psolver = be.asarray(repr(parent.solver)).astype("S")
            pcmp = be.asarray(repr(parent.cmp_)).astype("S")
            pnumvar = be.asarray(repr(parent.system.numvar)).astype("S")

            g.require_dataset(
                sysDTypes["system"], psystem.shape, psystem.dtype, data=psystem
            )
            g.require_dataset(sysDTypes["grid"], pgrid.shape, pgrid.dtype, data=pgrid)
            g.require_dataset(
                sysDTypes["solver"], psolver.shape, psolver.dtype, data=psolver
            )
            g.require_dataset(sysDTypes["cmp"], pcmp.shape, pcmp.dtype, data=pcmp)
            g.require_dataset(
                sysDTypes["numvar"], pnumvar.shape, pnumvar.dtype, data=pnumvar
            )

        def __call__(self, it, u):
            pass
