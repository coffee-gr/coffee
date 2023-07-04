"""The io module provides a wrapper around h5py so as to translate between
system, solver, and timeslice objects in a COFFEE friendly way.

Datagroups in an hdf file are assumed to contain all data associated to a 
simulation. Each dataset in a data group 

The module hdf_output is closely related to this module. The dictionary
dgTypes maps data group names to SimOutputType class names. It is best
to read the hdf_output module documentation as well as this modules
documentation if the dgTypes dictionary needs to be altered.
"""

from builtins import zip
from builtins import map
from builtins import str
from builtins import range
from builtins import object
import functools
import h5py
from coffee.settings import be

NUMERICAL_TOLERANCE = 1e-14
# Important configuration is included after these two classes


# A wrapper class that helps ease iteration over datasets with
# str(int) indices going from 0,1,... upwards.
class DataGroup(object):
    """The DataGroup class wraps a h5py group so that the setter, getter
    and iter methods do sensible array-like things.

    Datagroups are assumed to represent an entire simulation. Datasets
    are assumed to represent each timeslice. Each dataset in a data group
    is the value of a keyword "i" where i is the iteration number of the
    timeslice that the dataset contains.

    This class abstracts these details so that a given a data group ``g`` we
    can write ``g[2]'' to get the timeslice of the second iteration
    and ``g[3]`` to get the timeslice of the following iteration.
    """

    global group

    @property
    def attrs(self):
        """Wraps the attrs dictionary of a hdf datagroup."""
        return self.group.attrs

    @property
    def name(self):
        """Wraps the hdf datagroup name."""
        return self.group.name

    def attrs_list(self, kwd):
        """Represents the attributes of each data set of the data group
        as a list.

        Returns
        =======
        list :
            Each entry in the list is a value of some keyword in the
            attrs dictionary of a sub data set.
        """
        list = []
        for data_set in self:
            list += data_set.attrs[kwd]
        return list

    def index_of_attr(
        self, attr, value, start_index=0, value_comparor=lambda x: x == value
    ):
        """Returns the index of the dataset in self whose attribute attr
        has the given value.

        The function value_comparor allows for fudging a
        little. It is a function that takes the attrs of a data set and
        returns true or false. Once true is found the index of that data
        set is returned.

        Parameters
        ==========
        attr :
            A keyword to be applied to the attrs attribute of each dataset.
        value :
            The value that is searched for.
        start_index : int, Optional
            The index at which to start the search.
        value_comparor : function(value)
            A function which returns true if the given value matched the
            desired value.

        Returns
        =======
        int :
            The index of the dataset whose attrs attributes matched.
        """
        index = -1
        for i in range(start_index, len(self)):
            if value_comparor(self[i].attrs[attr]):
                index = self[i].attrs["index"]
                break
        return index

    def __init__(self, grp, returnValue=False):
        """The group to behave like an array. It is assumed that
        the group has/will have a number of datasets with the labels
        '0','1', etc...

        Parameters
        ==========
        grp : h5py.Group
            The data group that this object will wrap.
        returnValue : bool
            If true __getitem__ will return the value of the dataset rather
            than the dataset itself.
        """
        self.group = grp
        self.rV = returnValue

    def __iter__(self):
        """Iterates in increasing numerical order 0,1,2,3,...
        across the datasets of the group. Returns the dataset
        at position 0,1,2,3...

        A generator that allows for iterative access to datasets.

        Yields
        ======
        h5py.dataset :
        """
        i = 0
        while True:
            try:
                yield self[i]
                i += 1
            except:
                return

    def __len__(self):
        """Returns the number of datasets.

        Returns
        =======
        int :
            The number of datasets in this group.
        """
        return len(self.group)

    def __setitem__(self, i, value):
        """Allows specification of the value of a dataset at a given index.

        Parameters
        ==========
        i : int
            The index at which value should be placed
        value :
            The value of the data set at index i. Must be an object that
            h5py can place into a dataset.
        """
        value = be.array(value)
        dataset = self.group.require_dataset(str(i), value.shape, value.dtype)
        dataset[:] = value
        dataset.attrs["index"] = i

    def __getitem__(self, i):
        """Returns either the dataset or the value at index i.

        Parameters
        ==========
        i : int
            The index of the dataset to be returned.

        Returns
        =======
        dataset or value :
            The value will be returned only if the returnValue variable is
            true
        """
        if self.rV:
            return self.group[str(i)][()]
        return self.group[str(i)]

    def __repr__(self):
        return r"<H5pyArray datagroup %s (%d)>" % (self.name, len(self))


class DomainDataGroup(DataGroup):
    """A DataGroup wrapper that handles the specific case of function data
    over the computational grid.
    """

    def __setitem__(self, i, value):
        value = be.array(value)
        dataset = self.group.require_dataset(str(i), value.shape, value.dtype)
        dataset[:] = value
        dataset.attrs["index"] = i

    def __getitem__(self, i):
        dataset = self.group[str(i)]
        axes_shape = dataset.attrs["axes_shape"]
        axes = []
        start = 0
        for i in range(len(axes_shape)):
            axes += [dataset[()][start : start + axes_shape[i]]]
            start = start + axes_shape[i]
        return axes

    def __repr__(self):
        return r"<H5pyArray domaindatagroup %s (%d)>" % (self.name, len(self))


# Do not change the keys!
#
# Each dgTypes (DataGroup Types) describes the group name for a data group
# structure in the hdf file. Each data group stores hdf datasets named by
# an index. The index across different dgtypes gives data for the corresponding
# iteration.
#
# The keys represent the names of the different dgtypes and are hard coded
# throughout the code. Therefore if you wish to rename a dgtype, just change
# the item, not the key.
#
# These represent the types of data that a simulation knows about.

dgTypes = {
    "raw": "Raw_Data",
    "constraints": "Constraints",
    "exact": "Exact_Data",
    "errorNum": "Error_Numeric",
    "errorExa": "Error_Exact",
    "IJ": "Weyl_Constants_IJ",
    "domain": "Domain",
    "time": "Time",
    "dt": "Time_Step",
    "scrif": "Scri+",
    "constraint": "Constraint",
    "mu": "mu",
    "mup": "mup",
}

"""The dgTypes dictionary maps data types that are produced during simulation
to the keys used in the hdf file to store that data. The dgTypesInv 
dictionary provide the reverse mapping.

If you need a new data type to be written out you can dynamically modify this
dictionary. Changes will also need to be reflected in the hdf output action."""

dgTypesInv = dict(list(zip(list(dgTypes.values()), list(dgTypes.keys()))))


dgTypes_DataGroups = {"domain": (None, DomainDataGroup)}

"""The dgTypes_DataGroups dictionary maps keys in the dgTypes dictionary to
a 2-tuple. 

The first entry in the two tuple is a module name that is imported
using  __import__, the second entry is the name of a class in that module that
provides DataGroup functionality suitable for the given dgType.

If the first entry is none then it is assumed that the appropriate data group
is in this module. The class name rather than the class string should be in
the second entry of the tuple.

If a dgType key is not in the dictionary the DataGroup class is used.
"""

# SystemDataTypes stores a list of all the subgroups in the system groups.

systemD = "System"

"""The systemD variable store the hdf key used to create datagroups that
store data for system objects."""

sysDTypes = {
    "system": systemD,
    "solver": "Solver",
    "grid": "Grid",
    "cmp": "cmp",
    "numvar": "NumVariables",
}

"""The sysDTypes dictionary lists the data types in a system object (keys)
against the key used in the hdf file to store that data (values).

Feel free to dynamically alter the dictionary. Logic is based on correctness
of the keys not the values. New data types for system objects should be added
in this dictionary as key value pairs.
"""


# An interface to ease interaction with the simulationHDF class when
# only a specific simulation is wanted. I expect this class to be used the
# most.
@functools.total_ordering
class Sim(object):
    """
    Represents the data associated to a simulation as stored in an HDF file.

    Designed to be accessed via a SimulationHDF object.
    """

    def __init__(self, simName, simHDF):
        """The initialiser for Sim objects.

        Parameters
        ==========
        simName : string
            The name of the simulation this object will represent.
        simHDF : SimulationHDF
            The SimulationHDF object that wraps the hdf file.
        """
        self.simHDF = simHDF
        self.name = simName
        existing_items = list(self.simHDF[systemD + "/" + self.name].keys())
        for key, item in list(sysDTypes.items()):
            if item in existing_items:
                setattr(self, key, self.simHDF[systemD + "/" + self.name][item][()])
        self.cmp = float(self.cmp)
        self.numvar = int(self.numvar)
        existing_items = list(self.simHDF.file.keys())
        for key, item in list(dgTypes.items()):
            if item in existing_items:
                if self.name in list(self.simHDF[item].keys()):
                    if key in dgTypes_DataGroups:
                        if dgTypes_DataGroups[key][0] is not None:
                            mod = __import__(
                                dgTypes_DataGroups[key][0],
                                fromlist=[dgTypes_DataGroups[key][1]],
                            )
                            dataGroup = getattr(mod, dgTypes_DataGroups[key][1])
                            setattr(
                                self,
                                key,
                                dataGroup(
                                    self.simHDF[item + "/" + self.name],
                                    returnValue=True,
                                ),
                            )
                        else:
                            setattr(
                                self,
                                key,
                                dgTypes_DataGroups[key][1](
                                    self.simHDF[item + "/" + self.name],
                                    returnValue=True,
                                ),
                            )
                    else:
                        setattr(
                            self,
                            key,
                            DataGroup(
                                self.simHDF[item + "/" + self.name], returnValue=True
                            ),
                        )
        setattr(
            self,
            "indices",
            sorted(
                map(
                    int,
                    list(
                        self.simHDF[
                            list(self.simHDF.file.keys())[0] + "/" + self.name
                        ].keys()
                    ),
                )
            ),
        )

    def tslice(self, i):
        """Wraps the SimulationHDF method of the same name.

        Assumes that the simulation used is the one represented by this object.

        Parameters
        ===========
        i : int
            The iteration number of the timeslice to return.
        sim: string
            The name of the simulation to access
        dgType : string, Optional
            The name of the type of data to retrieve. See the dgTypes dictionary
            for a list of possible values.

        Returns
        =======
        tslice.TimeSlice:
        """
        return self.simHDF.tslice(i, self.name, dgType=dgType["raw"])

    def indexOfTime(self, t):
        """Wraps the SimulationHDF method of the same name.

        Assumes that the simulation used is the one represented by this object.

        Parameters
        ==========
        t : float
            The time whose index is desired.
        sim : string
            The name of the simulation to search.

        Returns
        =======
        int:
            The iteration index that matches the given time up to the
            assumed NUMERICAL_TOLERANCE.
        """
        return self.simHDF.indexOfTime(t, self.name)

    def __eq__(self, other):
        """Rich comparison based on the cmp parameter.

        self.cmp is given by the comparison parameter specified during simulation.

        Returns
        =======
        bool:
        """
        return self.cmp == other.cmp

    def __lt__(self, other):
        """Rich comparison based on the cmp parameter.

        self.cmp is given by the comparison parameter specified during simulation.

        Returns
        =======
        bool:
        """
        return self.cmp < other.cmp

    def __str__(self):
        return self.name

    def write(self, dgType, it, data, name=None, derivedAttrs=None):
        """Wraps the SimulationHDF file of the same name.

        Assumes that the simulation name is given by this object.

        Parameters
        ==========
        dgType : string
            A key from the module level dgTypes dictionary.
        it : int
            The iteration number of the data to be written.
        data :
            The data to write. It must be able to be stored in an h5py.dataset.
        name : string, Optional
            A parameter used to create a sub-datagroup. See comments above.
        derivedAttrs : dictionary, Optional
            A dictionary of additional attributes to add to created datasets
        """
        self.simHDF.write(dgType, self.name, it, data, name, derivedAttrs)

    def getDgType(self, dgType):
        """Wraps the SimulationHDF method of the same name.

        Parameters
        ==========
        dgType: string
            A key from the dgTypes dictionary.

        Returns
        =======
        DataGroup:
        """
        return self.simHDF.getDgType(dgType, self.name)

    def getDgTypeAttr(self, dgType, attr, i):
        """Wraps the SimulationHDF class' method of the same name.

        Parameters
        ==========
        dgType : string
            A key of the dgTypes dictionary. The dgType whose attributes are desired.
        attr : string
            The attribute to return.
        i : int
            The index of the iteration whose attributes are desired.

        Returns
        =======
        DataGroup :
        """
        return self.simHDF.getDgTypeAttr(dgType, attr, i, self.name)

    class dsReturnValue(object):
        def __init__(self, dataset):
            self.ds = dataset

        def __getitem__(self, key):
            return self.ds[key][()]


# Allows for interaction with the hdf file without specific reference
# to a particular simulation. I expect that this class will only be used
# for easy access to the sim objects.
class SimulationHDF(object):
    """
    Represents the data associated to all simulations as stored in an HDF file.
    """

    def __init__(self, fileName, **kwds):
        """The initialiser for SimulationHDF.

        Parameters
        ==========
        filename: string
            The name of the hdf file to be wrapped.
        """
        self.file = h5py.File(fileName, "r+")

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.file.close()

    def sim(self, name):
        """Returns a Sim class wrapping the data for the appropriate simulation."""
        return Sim(name, self)

    def name(self):
        """The name of the file that this class wraps.

        Returns
        =======
        string :
            The filename.
        """
        return self.file.name

    def getSims(self):
        """Returns a list of Sim objects representing the simulations contained
        in the hdf file.

        Returns
        =======
        list of Sim classses:
            The Sim classes wrapping the simulation contained in this hdf file.
        """
        simArray = [self.sim(name) for name in list(self.file[systemD].keys())]
        return sorted(simArray)

    def __getitem__(self, key):
        """Provides direct access to the underlying data group with the given key.

        Returns
        =======
        h5py.Group:
        """

        return self.file[key]

    def getSimData(self, sim):
        """Returns a dictionary giving access to the hdf objects corresponding
        to the given simulation name.

        Parameters
        ==========
        sim : string
            The name of the simulation whose data to retrieve.

        Returns
        =======
        dictionary:
        """
        sg = self.file[dgTypes["system"] + sim]
        rl = {"name": sim}
        for key, item in systemData:
            rl[key] = sg[item][()]
        return rl

    def tslice(self, i, sim, dgType="raw"):
        """The timeslice of the given simulation at the given iteration.

        Parameters
        ===========
        i : int
            The iteration number of the timeslice to return.
        sim: string
            The name of the simulation to access
        dgType : string, Optional
            The name of the type of data to retrieve. See the dgTypes dictionary
            for a list of possible values.

        Returns
        =======
        tslice.TimeSlice:
        """
        rdata = self.file[dgTypes[dgtype] + sim][str(i)]
        time = self.file[dgTypes["times"] + sim][str(i)]
        domain = self.file[dgTypes["domains"] + sim][str(i)]
        return tslice(rdata, domain, time)

    def getDgType(self, dgType, sim):
        """Returns a DataGroup object wrapping the given simulations data
        for the given dgType.

        Parameters
        ==========
        dgType: string
            A key from the dgTypes dictionary.
        sim: string
            The name of the simulation whose data to access.

        Returns
        =======
        DataGroup:
        """
        return DataGroup(self.file[dgTypes[dgType] + "/" + sim])

    def getDgTypeAttr(dgType, attr, i, sim):
        """Returns a DataGroup object wrapping the given simulations data set
        attributes for the given dgType.

        Parameters
        ==========
        dgType: string
            A key from the dgTypes dictionary.
        sim: string
            The name of the simulation whose data to access.

        Returns
        =======
        DataGroup:
        """
        return DataGroup(self.file[dgTypes[dgType] + "/" + sim]).attr[attr]

    def indexOfTime(self, t, sim):
        """A utility method that supports the translation of a simulation time
        to the iteration index.

        The tolerance for the required float comparison is given in
        the module variable NUMERICAL_TOLERANCE.

        Parameters
        ==========
        t : float
            The time whose index is desired.
        sim : string
            The name of the simulation to search.

        Returns
        =======
        int:
            The iteration index that matches the given time up to the
            assumed NUMERICAL_TOLERANCE.
        """
        times_dg = DataGroup(self.file[dgTypes["time"] + "/" + sim])
        indices = sorted([int(index) for index in list(times_dg.group.keys())])
        # If only one index
        if len(indices) == 1:
            if abs(times_dg[indices[0]][(0)] - t) <= NUMERICAL_TOLERANCE:
                return indices[0]
            else:
                return -1

        # Check initial step
        time_dg = times_dg[indices[0]]
        if indices[1] - indices[0] == 1:
            dt = times_dg[indices[1]][()] - time_dg[()]
            if t <= time_dg[()] < t + (dt / 2.):
                return indices[0]
        else:
            if abs(time_dg[(0)] - t) <= NUMERICAL_TOLERANCE:
                return indices[0]

        # Check all other steps except final
        for i in range(1, len(indices) - 1):
            time_dg = times_dg[indices[i]]
            if indices[i + 1] - indices[i] == 1:
                dt = times_dg[indices[i + 1]][()] - time_dg[()]
                if t - (dt / 2.) <= time_dg[()] < t + (dt / 2.):
                    return indices[i]
            else:
                if abs(time_dg[(0)] - t) <= NUMERICAL_TOLERANCE:
                    return indices[i]

        i = len(indices) - 1
        time_dg = times_dg[indices[i]]
        if indices[i] - indices[i - 1] == 1:
            dt = time_dg[()] - times_dg[indices[i - 1]]
            if t - (dt / 2.) <= time_dg[()] <= t:
                return indices[i]
        else:
            if abs(time_dg[(0)] - t) <= NUMERICAL_TOLERANCE:
                return indices[i]
        return -1

    def write(
        self, dgType, sim, it, data, name=None, derivedAttrs=None, overwrite=True
    ):
        """This method allows for writing to SimulationHDF objects.

        Note that if
        dgType is an error type data group then name must be given. We recommend
        that its value be taken as the data group from which the error data was
        generated.
        To ensure that name is not / is needed refer to how the data is
        extracted.

        Parameters
        ==========
        dgType : string
            A key from the module level dgTypes dictionary.
        sim : string
            The name of the simulation to write to.
        it : int
            The iteration number of the data to be written.
        data :
            The data to write. It must be able to be stored in an h5py.dataset.
        name : string, Optional
            A parameter used to create a sub-datagroup. See comments above.
        derivedAttrs : dictionary, Optional
            A dictionary of additional attributes to add to created datasets
        overwrite: bool
            If true data will be overwritten, otherwise an error will be
            raised.
        """
        # Create empy derivedAttrs if no argument is passed
        if derivedAttrs is None:
            self.derivedAttrs = {}
        else:
            self.derivedAttrs = derivedAttrs

        # If dgType is an error type then name must be set.
        # if dgType == sd.dgTypes["errorNum"] or\
        #    dgType == sd.dgTypes["errorExa"]:
        #    if name is None:
        #        raise Exception("""If SimulationHDF.write() is based an error dgType the name keyword must be set. Suggested usage is that name = the dgType of the data from which the error data was calculated.""")

        # get name if not none
        if name is not None:
            dg_name = dgType + "/" + sim + "/" + name
        else:
            dg_name = dgType + "/" + sim

        # get dg
        if overwrite:
            if name is not None:
                dg = DataGroup(
                    self.file.require_group(dgType)
                    .require_group(sim)
                    .require_group(name)
                )
            else:
                dg = DataGroup(self.file.require_group(dgType).require_group(sim))
        else:
            if name is not None:
                dg = DataGroup(
                    self.file.create_group(dgType).create_group(sim).create_group(name)
                )
            else:
                dg = DataGroup(self.file.create_group(dgType).create_group(sim))

        # add data and derived attrs
        dg[it] = data
        for key, value in list(self.derivedAttrs.items()):
            dg[it].attrs[key] = value


def array_value_index_mapping(correct, comparison, compare_on_axes=1):
    """A utility function which is useful when comparing data.

    The function takes two arrays and returns a list of pairs
    of indices (index1, index2) so that correct[index1] = comparison[index2]
    (up to NUMERICAL_TOLERANCE)
    this is very useful when performing error calculations.

    The idea is that data over two arbitrary domains can be easily compared.

    Parameters
    ==========
    correct : numpy.ndaray
    comparison : numpy.ndarray
    compare_on_axes : int
        If 0 then the first axes is included in the indices

    Returns
    =======
    list of 2-tuples of tuples of ints:
        Each pair of tuples of ints represents the data values which are
        defined over the same point.
    """
    index_mapping = []
    cor_dims = len(correct.shape)
    com_dims = len(comparison.shape)
    if compare_on_axes == 0:
        cor_ind = [0 for i in range(cor_dims)]
        com_ind = [0 for i in range(com_dims)]
    else:
        cor_ind = [0 for i in range(cor_dims - 1)]
        com_ind = [0 for i in range(com_dims - 1)]
    return _array_value_index_mapping_recursive(
        correct, cor_ind, comparison, com_ind, index_mapping, compare_on_axes, 0
    )


def _array_value_index_mapping_recursive(
    correct, cor_ind, comparison, com_ind, index_mapping, compare_on_axes, depth
):
    while (
        cor_ind[depth] < correct.shape[depth]
        and com_ind[depth] < comparison.shape[depth]
    ):
        if compare_on_axes == 0:
            com = comparison[tuple(com_ind)]
            cor = correct[tuple(cor_ind)]
        else:
            com = comparison[tuple(com_ind)][depth]
            cor = correct[tuple(cor_ind)][depth]
        if com + NUMERICAL_TOLERANCE < cor:
            com_ind[depth] = com_ind[depth] + 1
        elif com > NUMERICAL_TOLERANCE + cor:
            cor_ind[depth] = cor_ind[depth] + 1
        elif abs(com - cor) < NUMERICAL_TOLERANCE:
            if depth == len(correct.shape) - 1 - compare_on_axes:
                index_mapping += [(tuple(cor_ind), tuple(com_ind))]
                com_ind[depth] = com_ind[depth] + 1
                cor_ind[depth] = cor_ind[depth] + 1
            else:
                index_mapping = _array_value_index_mapping_recursive(
                    correct,
                    cor_ind,
                    comparison,
                    com_ind,
                    index_mapping,
                    compare_on_axes,
                    depth + 1,
                )
                com_ind[depth] = com_ind[depth] + 1
                cor_ind[depth] = cor_ind[depth] + 1
                com_ind[depth + 1] = 0
                cor_ind[depth + 1] = 0
        else:
            raise Exception("Unable to compare %s and %s" % (com, cor))
    return index_mapping
