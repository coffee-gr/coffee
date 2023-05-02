"""This module makes use of the simulation_data class in coffee.io to support
output / input of salm objects in hdf files.
"""
# Standard Library imports
import numpy as np

# Coffee imports
from coffee.actions.hdf_output import SimOutput
from coffee.io import simulation_data
from coffee.swsh.spinsfastpy import salm

simulation_data.dgTypes["sralm"] = "sralm"


class Sralm_Out(SimOutput.SimOutputType):
    """An actions.hdf_output.SimOutput.SimOutputType object for the writing of sralm objects to
    hdf files.
    """

    groupname = simulation_data.dgTypes["sralm"]

    def __call__(self, it, u):
        dg = self.data_group
        _write_salm(dg, it, u.data)
        super(Sralm_Out, self).__call__(it, u)


class SralmDataGroup(simulation_data.DataGroup):
    """A data group object for sralm objects."""

    def __init__(self, *args, **kwds):
        """Initialisation of SralmDataGroup.

        See simulation_data.DataGroup.__init__() for parameters.
        """
        super(SralmDataGroup, self).__init__(*args, **kwds)

    def __setitem__(self, i, sralm):
        """Supports writing of data.

        See simulation_data.DataGroup.__setitem__() for parameters.
        """
        if len(sralm.shape) == 1:
            _write_salm(self.group, i, sralm)

    def __getitem__(self, i):
        """Supports reading of data.

        See simulation_data.DataGroup.__getitem__() for parameters.
        """
        if self.rV:
            cg_mod = __import__(
                self.group[str(i)].attrs["cg_module"],
                fromlist=[self.group[str(i)].attrs["cg_class"]],
            )
            cg_class = getattr(cg_mod, self.group[str(i)].attrs["cg_class"])
            return salm.sfpy_sralm(
                self.group[str(i)].value,
                self.group[str(i)].attrs["spins"],
                int(self.group[str(i)].attrs["lmax"]),
                cg=cg_class(),
                bandlimit_multiplication=self.group[str(i)].attrs["bl_mult"],
            )
        return self.group[str(i)]

    def __repr__(self):
        return r"<SralmDataGroup %s (%d)>" % (self.name, len(self))


def _write_salm(datagroup, index, salm):
    """A utility method that performs the writing of salm objects to hdf
    datagroups.

    Parameters
    ----------
    datagroup : hdf.datagroup
        The hdf data group to which the salm object will be written.
    index : int
        The index of the dataset in the data gropu that will contain the data.
    salm: swsh.salm.salm
        The spin weighted spherical harmonic components.
    """
    datagroup[str(index)] = be.asarray(salm)
    datagroup[str(index)].attrs["spins"] = be.array(salm.spins)

    lmax = salm[0].lmax
    for i in range(1, salm.shape[0]):
        if lmax is not salm[0].lmax:
            raise ValueError(
                "Unable to store salm objects with different\
                lmax in the same datagroup."
            )
    datagroup[str(index)].attrs["lmax"] = lmax

    cg = salm[0].cg
    for i in range(1, salm.shape[0]):
        if cg is not salm[0].cg:
            raise ValueError(
                "Unable to store salm objects with different\
                cg in the same datagroup."
            )
    datagroup[str(index)].attrs["cg_module"] = cg.__module__
    datagroup[str(index)].attrs["cg_class"] = cg.__class__.__name__

    bl_mult = salm[0].bl_mult
    for i in range(1, salm.shape[0]):
        if bl_mult is not salm[0].bl_mult:
            raise ValueError(
                "Unable to store salm objects with different\
                multiplication bandlimits in the same datagroup."
            )
    datagroup[str(index)].attrs["bl_mult"] = bl_mult


# Dynamic injection of the salm datagroup into the coffee.io module.
simulation_data.dgTypes_DataGroups["sralm"] = (
    "coffee.swsh.spinsfastpy.io",
    "SralmDataGroup",
)
