#! /usr/bin/env python


import argparse

from coffee.io import simulation_data as sd


def printTimes(simulationHDF):
    sims = simulationHDF.getSims()
    if args.s:
        sims = [sims[0]]
    for sim in sims:
        print("===========================")
        print(sim)
        print("Times: \n%s" % repr(sorted([sim.time[i] for i in sim.indices])))
        print("Time Step: \n%s" % repr(sorted([sim.dt[i] for i in sim.indices])))


def _validate(group):
    """If group is an instance of sd.SimulationHDF then group.file is return.
    Otherwise group is returned. This method validates the arguments of
    exploreKeys() and printHDF() to ensure that they are acting on
    HDf groups."""
    if isinstance(group, sd.SimulationHDF):
        return group.file
    return group


def exploreKeys(group):
    """Takes an HDF group and prints the group name followed by the group keys()"""
    group = _validate(group)
    print(group)
    print(list(group.keys()))


def printHDF(group):
    """Prints the first two levels of the passed HDF group. If the group is a simulation_data.py HDF file then the "System Data" and "Numerical Error" groups are further explored.

    group = hdf group to be printed"""
    group = _validate(group)
    print("===========================")
    exploreKeys(group)
    for key1 in list(group.keys()):
        print("===========================")
        exploreKeys(group[key1])
        for key2 in list(group[key1].keys()):
            if key1 == sd.dgTypes["errorNum"]:
                print("===========================")
                exploreKeys(group[key1][key2])
            elif key1 == sd.systemD:
                print("===========================")
                exploreKeys(group[key1][key2])
                for key3 in group[key1][key2]:
                    print(
                        "%s: %s"
                        % (group[key1][key2][key3], group[key1][key2][key3].value)
                    )
            elif key1 == sd.dgTypes["time"]:
                print(group[key1][key2])
                for key3 in group[key1][key2]:
                    print(
                        "%s: %s"
                        % (group[key1][key2][key3], group[key1][key2][key3].value)
                    )
            else:
                print(group[key1][key2])


def loadHDF(file):
    """Return a simulationHDF() object representing the HDF file. The original HDF file can be accessed via <returned object>.file. SimulationHDF() objects have a number of useful methods, including simplified access to individual runs via the .getSims() method.

    file = string giving location of the hdf file."""
    return sd.SimulationHDF(file)


def mathematica_output(array):
    s = "{"
    for i, v in enumerate(array):
        if len(v.shape) >= 1:
            s = s + mathematica_output(v)
            if i < array.shape[0] - 1:
                s = s + ","
        elif len(v.shape) == 0:
            if i < array.shape[0] - 1:
                s = s + "%.16f, " % v
            else:
                s = s + "%.16f" % v
    s = s + "}"
    return s


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This program outputs the keys and contents of the first two levels of an hdf file. If the file is produced from simulation_data.py then the "Numerical Error" and "System Data" groups 
    are further explored. The numerical values of data groups are not returned. It does not currently provide information on the attributes of the HDF groups.

This script is also designed to be loaded as a module. If done three functions, exploreKeys(), loadHDF() and printHDF() are made available. Remember to appropriately close the hdf file, otherwise corruption can result."""
    )
    parser.add_argument(
        "file", help="""The hdf file whose contents is to be printed to stdout."""
    )

    parser.add_argument(
        "-t",
        "-time",
        action="store_true",
        default=False,
        help="""Prints the start time, stop time and time step of each run.""",
    )

    parser.add_argument(
        "-s",
        "-small",
        action="store_true",
        default=False,
        help="""If set the script will only print the details for the smallest run as
    given by the cmp parameter.""",
    )

    parser.add_argument(
        "-g",
        "-groups",
        action="store_true",
        default=False,
        help="""Prints keys and group names.""",
    )

    args = parser.parse_args()

    with sd.SimulationHDF(args.file) as file:
        if args.t:
            printTimes(file)
        if args.g:
            printHDF(file.file)
