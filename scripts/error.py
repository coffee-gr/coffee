#! /usr/bin/env python


import argparse
import sys
from builtins import range
from builtins import str
from builtins import zip

# from coffee.swsh.spinsfastpy import io as swsh_io
from coffee.io import simulation_data as sd


def exact(args):
    if args.dg != ["raw"]:
        print("Exact errors can only be calculated for data type 'raw'.")
        sys.exit(1)

    with sd.SimulationHDF(args.file) as file:
        print("Initialising data.")
        sims = file.getSims()

        # needed for fixing domain size to smallest domain. This will cause
        # problems when plotting.
        # smallest_domain = sims[0].domain

        tSimNames = [sim.name for sim in sims]
        stepsizes = []
        # We make the assumption that the number of components for each run
        # is the same and that the values of the components are in the first
        # axis of raw.shape
        num_of_comps = sims[0].numvar
        tableE = be.zeros((len(args.Lp), len(tSimNames), num_of_comps), dtype=float)

        for time in args.t:
            print("==================================")
            print("Doing calculation for time = %f" % time)
            for i, sim in enumerate(sims):
                # print "Doing calculation for simulation %s"%sim.name
                it = sim.indexOfTime(time)
                domain = sim.domain[it]

                # needed for fixing domain size to smallest domain.
                # This will cause
                # problems when plotting.
                # if i == 0:
                # smallest_it = it
                # smallest_domain = domain

                print("Index for simulation %s is %i at time %f" % (sim.name, it, time))

                # needed for fixing domain size to smallest domain.
                # This will cause
                # problems when plotting.
                # if __debug__:
                # print "Index for smallest_domain is %i at time %f"%(
                # smallest_it, time)
                # print "smallest domain is: %s"%repr(smallest_domain)
                # axes_mappings = [
                # sd.array_value_index_mapping(
                # smallaxis,
                # axis,
                # compare_on_axes = 0
                # )
                # for smallaxis, axis in
                # zip(smallest_domain, domain)
                # ]
                # if __debug__: print "Mapping = "+str(axes_mappings)

                # Calculating difference in values of common domains

                # needed for fixing domain size to smallest domain.
                # This will cause
                # problems when plotting.
                # error = be.ones_like(sims[0].raw[smallest_it])
                # for from_tup, to_tup in map_generator(axes_mappings):
                # error[(slice(None),)+from_tup] = be.absolute(
                # sim.raw[it][(slice(None), ) + to_tup] -
                # sim.exact[it][(slice(None), ) + to_tup]
                # )

                error = be.absolute(sim.raw[it] - sim.exact[it])
                if __debug__:
                    print("Errors are: %s" % repr(error))
                sim.write(sd.dgTypes["errorExa"], it, error)

                # we assume that stepsizes is constant for each slice
                stepsizes += [be.asarray([axis[1] - axis[0] for axis in domain])]
                for j, p in enumerate(args.Lp):
                    tableE[j][i] = Lp(error, stepsizes[-1], p)
            for j, p in enumerate(args.Lp):
                rows = _printErrorConv(
                    tSimNames, list(range(num_of_comps)), tableE[j], stepsizes, time, p
                )
                with open(
                    "%s-e-exa_L%f_%f.csv" % (args.ofile_base, p, time), "wb"
                ) as file:
                    for row in rows:
                        file.write(row)


def numer(args):
    with sd.SimulationHDF(args.file) as file:
        # Get all the current simulations
        sims = file.getSims()
        tSimNames = [sim.name for sim in sims[:-1]]
        stepsizes = []

        # We make the assumption that the number of components for each run
        # is the same and that the values of the components are in the first
        # axis of raw.shape
        num_of_comps = sims[0].numvar
        tableE = be.zeros((len(args.Lp), len(tSimNames), num_of_comps), dtype=float)

        # Seperate them into the one with the best resolution
        # and the others to be compared to this one.
        # Note that the code assumes that sims is appropriately ordered
        subtrahend = sims[-1]
        minuends = sims[:-1]

        #        # We must now check how the domains have been constructed.
        #        # We need to tell the difference between [[i] for i in range]
        #        # and [i for i in range].
        #        # We make the assumption that the domain type remains the same
        #        # throughout the simulation.
        #        meshes = subtrahend.domain[0]
        #        dims = len(subtrahend.domain.group['0'].attrs['shape'])
        #        if dims is 1:
        #            compare_on_axes = 0
        #        else:
        #            index = tuple([0 for i in domain.shape])
        #            if len(be.array(domain[index]).shape) is 0 and be.array(domain[index[:-1]]).shape[0] == dims-1:
        #                compare_on_axes = 1
        #            elif len(be.array(domain[index])) is dims:
        #                compare_on_axes = 0
        #            else:
        #                raise Exception("Unable to determine domain type")

        for time in args.t:
            print("==================================")
            print("Doing calculation for time = %f" % time)

            # Get data for the subtrahead, the index of time, the domain.
            subtrahend_index = subtrahend.indexOfTime(time)
            print(
                "Subtrahend time = %f at index = %i"
                % (subtrahend.time[subtrahend_index], subtrahend_index)
            )
            subtrahend_domain = subtrahend.domain[subtrahend_index]
            if __debug__:
                print("Subtrahend domain = %s" % str(subtrahend_domain))

            # Need to collect the same information for each minuend and do the
            # comparison.
            for i, minuend in enumerate(minuends):
                print("Calculating %s - %s" % (minuend.name, subtrahend.name))

                # Finding index for correct time
                minuend_index = minuend.indexOfTime(time)
                print(
                    "Minuend time = %f at index = %i"
                    % (minuend.time[minuend_index], minuend_index)
                )

                # Comparing domains
                minuend_domain = minuend.domain[minuend_index]
                if __debug__:
                    print("Minuend domain = %s" % str(minuend_domain))

                #                # Get the mapping between the domains.
                #                compare_on_axes = True
                #                if len(minuend_domain.shape) == 1:
                #                    compare_on_axes = False
                #                mapping = sd.array_value_index_mapping(minuend_domain,\
                #                    subtrahend_domain,compare_on_axes = compare_on_axes)
                axes_mappings = [
                    sd.array_value_index_mapping(
                        minuend_axis, subtrahend_axis, compare_on_axes=0
                    )
                    for minuend_axis, subtrahend_axis in zip(
                        minuend_domain, subtrahend_domain
                    )
                ]
                if __debug__:
                    print("Mapping = " + str(axes_mappings))

                for dgType in args.dg:
                    print("Calculating error for dgType %s" % dgType)

                    # Collecting dgType data
                    subtrahend_dg = getattr(subtrahend, dgType)
                    minuend_dg = getattr(minuend, dgType)

                    # Calculating difference in values of common domains
                    # diff = be.ones_like(minuend_dg[minuend_index])
                    diff = be.ones(
                        (minuend_dg[minuend_index].shape[0], len(axes_mappings[0]))
                    )
                    # import pdb; pdb.set_trace()
                    for index, from_tup, to_tup in map_generator(axes_mappings):
                        diff[:, index] = (
                            minuend_dg[minuend_index][(slice(None),) + from_tup]
                            - subtrahend_dg[subtrahend_index][(slice(None),) + to_tup]
                        )

                    # Storing data
                    minuend.write(
                        sd.dgTypes["errorNum"],
                        minuend_index,
                        be.absolute(diff),  # derivedAttrs = {sd.dgTypes['time']:time})
                        # name=sd.dgTypes[dgType]
                    )

                    # we assume that dx is constant for each slice
                    stepsizes += [
                        be.asarray([axis[1] - axis[0] for axis in minuend_domain])
                    ]
                    for j, p in enumerate(args.Lp):
                        tableE[j][i] = Lp(diff, stepsizes, p)
                # end loop over args.dg
            for j, p in enumerate(args.Lp):
                rows = _printErrorConv(
                    tSimNames, list(range(num_of_comps)), tableE[j], stepsizes, time, p
                )
                with open(
                    "%s-e-num_L%f_%f.csv" % (args.ofile_base, p, time), "wb"
                ) as file:
                    for row in rows:
                        file.write(row)


def map_generator(array, i=0):
    data = array[i]
    index = 0
    for datum in data:
        from_res = (datum[0][0],)
        to_res = (datum[1][0],)
        if i == len(array) - 1:
            yield index, from_res, to_res
            index += 1
        else:
            for res in map_generator(array, i + 1):
                yield index, from_res + res[0], to_res + res[1]
                index += 1


################################################################################
# Routines for numerical error estimation and output
################################################################################

# returns descrete version of lp norms


def Lp(errors, stepsizes, p):
    errors = be.absolute(errors)
    if p == float("infinity"):
        rerrors = be.max(errors)
    else:
        rerrors = be.power(
            be.squeeze(
                be.apply_over_axes(
                    be.sum, be.power(errors, p), list(range(1, len(errors.shape)))
                )
            ),
            1 / p,
        )
    return rerrors


# helper method to print error conv data
def _printErrorConv(Sims, Indices, errorData, stepSizes, time, p):
    rString = []
    firC = 10
    secC = 10
    thiC = 10
    print("")
    head1 = "{0:24} {1}".format("Simulation at t = %.4f" % time, "Error data")
    head2 = "{0:24} ".format("")
    head3 = head2
    head4 = ""
    for i in Indices:
        head2 = head2 + "{0:^26} ".format("Component %i" % i)
    for i in Indices:
        head3 = head3 + "{0:15}  {1:10}".format("log(L_%.4f,2)" % p, "Conv R")
    print(head1)
    rString += [head1 + "\n"]
    print(head2)
    rString += [head2 + "\n"]
    print(head3)
    rString += [head3 + "\n"]
    print(head4)
    rString += [head4 + "\n"]
    for i, simName in enumerate(Sims):
        s = "{0:24} ".format(stepSizes[i])
        for j, comp in enumerate(Indices):
            if not i == 0:
                s = s + "{0:<15.6f}  {1:<10.4f}".format(
                    be.log2(errorData[i][j]),
                    _conv(
                        errorData[i - 1][j],
                        stepSizes[i - 1],
                        errorData[i][j],
                        stepSizes[i],
                    )[0],
                )
            else:
                s = s + "{0:<15.6f}  {1:10}".format(be.log2(errorData[i][j]), "")
        print(s)
        rString += [s + "\n"]
    return rString


def _conv(old, numSteps1, new, numSteps2):
    return be.log(old / new) / be.log(numSteps1 / numSteps2)


################################################################################
# Main parser
################################################################################
parser = argparse.ArgumentParser(
    description="""This program calculates the error for a particular data type based on comparison to the highest resolution simulation or to exact data. The highest resolution is selected as the simulation with the highest cmp parameter. See simulation_hdf.py for the use of cmp."""
)

subparsers = parser.add_subparsers(
    title="Subcommands",
    description="""The subcommands
allow for selection of numerical or exact calculation of error.""",
    help="""For more information on any of the subcommands please use 'error <subcommand> -h""",
)


parser.add_argument(
    "file",
    help="""The hdf file, produced by simulation_data.py that contains the data to be
analysed""",
)
parser.add_argument(
    "-t",
    "-time",
    type=float,
    metavar="TIME",
    action="append",
    required=True,
    help="""A time at which the error analysis should be run. Multiple times can be given""",
)
parser.add_argument(
    "-dg",
    "-dgtype",
    action="append",
    help="""A data group to be analysed. See simulation_data for information on
the different data groups. Defaults to "raw". Multiple data groups can be given.""",
)
parser.add_argument(
    "-Lp",
    metavar="P",
    action="append",
    help="""The p value of the norm to be used to calculate the error. The infinity
norm can be entered by using -Lp "inf".""",
)
parser.add_argument(
    "-O",
    default=False,
    action="store_true",
    help="""If set debugging information will be printed.""",
)

################################################################################
# Exact subcommand
################################################################################
parser_exact = subparsers.add_parser(
    "exact", help="""This command calculates error using exact data."""
)
parser_exact.set_defaults(func=exact)

################################################################################
# Exact subcommand
################################################################################
parser_numer = subparsers.add_parser(
    "numer",
    help="""This command calculates error using the highest resolution, as determined by cmp,
numerical data.""",
)
parser_numer.set_defaults(func=numer)


args = parser.parse_args()
if args.dg is None:
    args.dg = ["raw"]
print("Starting computation...")
if args.Lp is None:
    args.Lp = ["2"]
for i in range(len(args.Lp)):
    args.Lp[i] = float(args.Lp[i])
args.ofile_base, sep, exten = args.file.rpartition(".")
args.func(args)
print("...computation complete.-")
