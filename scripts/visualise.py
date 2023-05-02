#! /usr/bin/env python


from builtins import str
from builtins import range
from past.utils import old_div
import sys
import os
import Gnuplot
import argparse
import numpy as np
import math

from coffee.io import simulation_data as sd

# load default file
from hdfvis_gnuplot_defaults import *


################################################################################
# Routine to plot error plots
################################################################################
def error(args):
    # Initialise Gnuplot instence
    g = Gnuplot.Gnuplot()
    if args.g is None:
        args.g = hdfvis_gnuplot_error
    if not args.d:
        args.g += hdfvis_gnuplot_error_nodis
    for com in args.g:
        g(com)
    args.ofile_ext = hdfvis_gnuplot_error_nodis_ext

    args.dg = ["raw"]

    with sd.SimulationHDF(args.file) as file:
        # Collect groups and datagroups
        sims = file.getSims()

        # Setting argument defaults
        if args.c is None:
            args.c = list(range(sims[0].numvar))

        # Initialise array to calculate abs errors and convergence rates
        tSimNames = [sim.name for sim in sims[:-1]]
        tPhiIndices = args.c
        tableE = [[[] for i in args.c] for j in tSimNames]
        stepSizes = [sim.cmp for sim in sims[:-1]]

        print("Producing graph...")

        # Get data for additional data
        additional_data = []
        add_index = sims[-1].indexOfTime(args.t)
        for data in args.add:
            data_dg = sims[-1].getDgType(data)[add_index]
            add_domain = sims[-1].domain[add_index]
            add_range = None
            # scrif = sims[-1].scrif[baIndex]
            # scrifIndex = be.nonzero(be.absolute(scrif)==\
            #    be.min(be.absolute(scrif)))[0][0]
            # scrifDomain = be.array([baDomain[scrifIndex],baDomain[scrifIndex]])
            # scrifRange = be.array([-50,0])
            additional_data += [(data, index, data_dg, domain)]

        # Producing plots
        for j, phi_index in enumerate(args.c):
            print("Doing component %i" % phi_index)
            plot_data = []
            for i, sim in enumerate(sims):
                if i == len(sims) - 1 and not args.e:
                    break
                # Get data for plot
                index = sim.indexOfTime(args.t)
                if args.e:
                    error = sim.errorExa[index][phi_index, :]
                else:
                    error = sim.errorNum[index][phi_index, :]
                domain = sim.domain[index]
                name = sim.name
                if args.o:
                    if args.e:
                        plot_data += [
                            Gnuplot.Data(
                                domain[0],
                                be.log2(error),
                                title=name,
                                filename="%s-e-exa_%i_%f.gnup"
                                % (args.ofile_base, phi_index, args.t),
                            )
                        ]
                    else:
                        plot_data += [
                            Gnuplot.Data(
                                domain[0],
                                be.log2(error),
                                title=name,
                                filename="%s-e-num_%i_%f.gnup"
                                % (args.ofile_base, phi_index, args.t),
                            )
                        ]
                else:
                    plot_data += [Gnuplot.Data(domain[0], be.log2(error), title=name)]
            if "scri" in list(
                sims[-1].raw.group[str(sims[-1].indexOfTime(args.t))].attrs.keys()
            ):
                scri = sim.raw.group[str(index)].attrs["scri"]
                g(
                    r"set arrow from %f,graph(0,0) to %f,graph(1,1) nohead"
                    % (scri, scri)
                )
            for add_name, add_index, add_data, add_domain in additional_data:
                plot_data += [
                    Gnuplot.Data(
                        baDomain,
                        baRaw[phi_index, :],
                        axes="x1y2",
                        title="Phi %i" % phi_index,
                    )
                ]
                plot_data += [Gnuplot.Data(scrifDomain, scrifRange, title="Scri+")]
            g("set terminal postscript enhanced color solid")
            g(
                r'set title "Errors in {/Symbol \152}_{%i} at time %s" enhanced'
                % (phi_index, args.t)
            )
            # g(r'set title "Errors in phi %i at time %s\nfrom %s" enhanced' %\
            # (phi_index,args.t,args.file))
            if not args.d:
                if args.e:
                    g(
                        'set output "%s-e-exa_%i_%f.%s"'
                        % (args.ofile_base, phi_index, args.t, args.ofile_ext)
                    )
                else:
                    g(
                        'set output "%s-e-num_%i_%f.%s"'
                        % (args.ofile_base, phi_index, args.t, args.ofile_ext)
                    )
            g.plot(*plot_data)
        g.close()
        print("Done")


# Manipulate args for plot data
def plot(args):
    if args.g is None:
        args.g = hdfvis_gnuplot_plot
    if not args.d:
        args.g += hdfvis_gnuplot_plot_nodis
    args.ofile_ext = hdfvis_gnuplot_plot_nodis_ext

    args.t0 = args.t
    args.t1 = args.t
    args.ofile_base = args.ofile_base + "-p"

    _plot(args)


# Manipulate args for animation data
def anima(args):
    if args.g is None:
        args.g = hdfvis_gnuplot_anima
    if not args.d:
        args.g += hdfvis_gnuplot_anima_nodis
    args.ofile_ext = hdfvis_gnuplot_anima_nodis_ext

    args.t = repr([args.t0, args.t1])
    args.ofile_base = args.ofile_base + "-a"
    args.frame = 0

    _plot(args)


################################################################################
# Routine for plots and animations
################################################################################
def _plot(args):
    animationLength = 10
    framesPerSec = 60

    # Initialize gnuplot
    g = Gnuplot.Gnuplot()
    # for com in args.g:
    #    g(com)

    with sd.SimulationHDF(args.file) as file:
        print("Producing graphs...")
        # Collect groups and datagroups
        sims = file.getSims()
        if args.s is None:
            args.s = list(range(len(sims)))

        for sim in (sims[i] for i in args.s):
            print("Doing simulation %s" % sim.name)
            # Get x values
            domains = sim.getDgType("domain")

            # Get all times
            times = sim.getDgType("time")
            times.rV = True
            numOfFrames = len(times)
            frameSkip = int(old_div(numOfFrames, (animationLength * framesPerSec)))

            # Get data for additional data
            additional_data = [(data, sim.getDgType(data)) for data in args.add]

            # Iterate across group
            # Get starting and stoping index
            nextFrame_index = sim.indexOfTime(args.t0)
            stop_index = sim.indexOfTime(args.t1)

            for dg in args.dg:
                print("Doing dgType %s" % dg)
                g.reset()
                for com in args.g:
                    g(com)
                g("set yrange [-5:5]")
                if not args.d:
                    g(
                        'set output "%s_%s_%s.%s"'
                        % (args.ofile_base, sim.name, args.t, args.ofile_ext)
                    )

                group = sim.getDgType(dg)
                # While there is a next frame...
                while nextFrame_index < stop_index + 1:
                    # plot data
                    i = nextFrame_index
                    y = group[i]
                    g.title("Simulation %s at time %f" % (sim.name, times[i]))
                    plotItems = []
                    for j, row in enumerate(be.atleast_2d(y.value)):
                        if args.o and args.frame is None:
                            plotItems += [
                                Gnuplot.Data(
                                    domains[i],
                                    row,
                                    title="Component " + str(j),
                                    filename="%s_%s_%s_Comp%i.gnu"
                                    % (args.ofile_base, sim.name, args.t, j),
                                )
                            ]
                        elif args.o and args.frame is not None:
                            plotItems += [
                                Gnuplot.Data(
                                    domains[i],
                                    row,
                                    title="Component " + str(j),
                                    filename="%s_%s_%s_Comp%iframe%i.gnu"
                                    % (
                                        args.ofile_base,
                                        sim.name,
                                        args.t,
                                        j,
                                        args.frame,
                                    ),
                                )
                            ]
                        else:
                            plotItems += [
                                Gnuplot.Data(
                                    domains[i], row, title="Component " + str(j)
                                )
                            ]
                    if args.o and args.frame is None:
                        for add_data, add_dg in additional_data:
                            plotItems += [
                                Gnuplot.Data(
                                    domains[i],
                                    add_dg[i],
                                    title=add_data,
                                    filename="%s-%s_%s_%s.gnup"
                                    % (args.ofile_base, add_data, sim.name, args.t),
                                )
                            ]
                    elif args.o and args.frame is not None:
                        for add_data, add_dg in additional_data:
                            plotItems += [
                                Gnuplot.Data(
                                    domains[i],
                                    add_dg[i],
                                    title=add_data,
                                    filename="%s-%s_%s_%s_frame%i.gnup"
                                    % (
                                        args.ofile_base,
                                        add_data,
                                        sim.name,
                                        args.t,
                                        args.frame,
                                    ),
                                )
                            ]
                    else:
                        for add_data, add_dg in additional_data:
                            plotItems += [
                                Gnuplot.Data(domains[i], add_dg[i], title=add_data)
                            ]
                    g.plot(*plotItems)
                    # if there are not enough frames left set the frameSkip to 0
                    if nextFrame_index + frameSkip >= numOfFrames:
                        frameSkip = 0
                    nextFrame_index += 1 + frameSkip
                    if args.frame is not None:
                        args.frame = args.frame + 1
                    # if the next frame is larger than the stop_index then stop
                    # plotting
                    if nextFrame_index > stop_index:
                        break
        print("...Done.-")


################################################################################
# Main parser
################################################################################
parser = argparse.ArgumentParser(
    description="""This program allows different types of display 
of data contained in hdf files generated
using simulation_data.py."""
)
subparsers = parser.add_subparsers(
    title="Subcommands",
    description="""The hdfvis command allows for several different type of graphs and 
animations of hdf files produced by simulation_data.py hdf file. Each of the subcommands allows access to a particular type of graph or animation.""",
    help="""For more information on any of the subcommands please use 'hdfvis <subcommand> -h""",
)

parser.add_argument(
    "file",
    help="""The hdf file, produced by simulation_data.py that contains the data to be
plotted.""",
)

parser.add_argument(
    "-dg",
    "-dgtype",
    action="append",
    help="""The data group type to be plotted. See simulation_data for information on the different data groups. Defaults to "raw". 
Multiple data types can be given by specifying multiple '-dg' options.""",
)

parser.add_argument(
    "-d",
    "-display",
    action="store_true",
    default=False,
    help="""This flag indicates that the visualisations are to be displayed on screen. No files will be saved. 
If not specified the gnuplot terminal and the saved file extension will be set to the values given in 'hdfvis_gnuplot_defaults.py'. 
This file can be editted to change the defaults.""",
)

parser.add_argument(
    "-g",
    "-gnuplot",
    metavar="COMMAND",
    action="append",
    help="""A list of gnuplot commands. These commands are provided to a Gnuplot.py object. If not specified the defaults are read in from 'hdfvis_gnuplot_defaults.py'. 
This file can be editted to change the defaults. Multiple gnuplot commands can be given by specifying multiple '-g' commands.""",
)

parser.add_argument(
    "-o",
    "-output",
    action="store_true",
    default=False,
    help="""A flag to indicate that the data used for the plot is to be written to the HDD.""",
)

parser.add_argument(
    "-add",
    "-additional_data",
    default=[],
    action="append",
    help="""A list of data types for plotting in addition to the data specified by the
selection of a subcommand. Please inspect each hdf file to see which data types
are availale for plotting.""",
)

################################################################################
# Parser for animations
################################################################################
parser_ani = subparsers.add_parser(
    "ani",
    help="""This command produces an animation of the data over the given times against the relevant domain.""",
)
parser_ani.set_defaults(func=anima)
parser_ani.add_argument(
    "-s",
    "-sims",
    type=int,
    metavar="INDEX",
    nargs="+",
    help="""A list of the simulations to display information for given by the index of the simulation in the sorted list of simulation given by CMP. 
If not set an animation or graph for each simulation is given. """,
)
parser_ani.add_argument(
    "-t0",
    "-tstart",
    metavar="TIME",
    default=-float("Infinity"),
    type=float,
    help="""The start time of the animation or plot. Defaults to the earlist possible time.""",
)
parser_ani.add_argument(
    "-t1",
    "-tstop",
    type=float,
    metavar="TIME",
    default=float("Infinity"),
    help="""The stop time of an animation. Defaults to the latest possible time.""",
)

################################################################################
# Parser for plots
################################################################################
parser_plot = subparsers.add_parser(
    "plot",
    help="""The subcommand 'normal' produces a graph or animation (based on -t0 and -t1) of the values of functions at a particular time 
(or across a range of times for an animation) against the domain of those functions""",
)
parser_plot.set_defaults(func=plot)
parser_plot.add_argument(
    "-s",
    "-sims",
    type=int,
    metavar="INDEX",
    nargs="+",
    help="""A list of the simulations to display information for given by the index of the simulation in the sorted list of simulations given by increasing CMP. 
If not set an animation or graph for each simulation is given. """,
)
parser_plot.add_argument(
    "-t",
    "-time",
    metavar="TIME",
    type=float,
    help="""The time at which the data will be plotted.""",
)

################################################################################
# Parser for error graphs
################################################################################
parser_error = subparsers.add_parser(
    "err",
    help="""The subcommand 'error' plots the precalculated error of a particular component of a particular data type against the domain of the component at a particular time. 
The graph also include a vertical line indicating the position of scri and the graph (against the right hand 'y' axis of the most accurate approximation of the function. 
WARNING: The error command can currently only use '-dg "raw"'. Any other '-dg' setting will default to "raw".""",
)
parser_error.set_defaults(func=error)
parser_error.add_argument(
    "-c",
    "-components",
    nargs="+",
    type=int,
    help="""A list of the components for which you want error plots produced. If not set a graph for each component will be made.""",
)
parser_error.add_argument(
    "-t",
    "-time",
    metavar="TIME",
    type=float,
    help="""The time at which the errors for the given component should be displayed.""",
)
parser_error.add_argument(
    "-e",
    "-exact",
    action="store_true",
    default=False,
    help="""A flag that indicates that the errors calculated against the exact data
should be graphed.""",
)


# Main program
args = parser.parse_args()
if args.dg is None:
    args.dg = ["raw"]
args.ofile_base, sep, exten = args.file.rpartition(".")
args.frame = None
args.func(args)
