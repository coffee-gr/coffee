#! /usr/bin/env python

import argparse
import subprocess
import sys
import os

# Initialise parser
parser = argparse.ArgumentParser(
    description="""This program runs a complete simulation: the numerical calculations, the error caculations and the visulations. This requires the simulation to output data in to an hdf file using actions.hdf_output."""
)

# Parse files
parser.add_argument(
    "-f", "-file", required=True, help="""The name of the hdf file to be produced."""
)

parser.add_argument(
    "-s",
    "-simulation",
    help="""The simulation to be run. If not specified -f must be used.""",
)

# Parse times
parser.add_argument(
    "-terrN",
    "-times-error-numeric",
    nargs="+",
    default=[],
    help="""A list of all times for which error calculations should be performed against the highest resolution numeric data. This list will also be used for visualisation of error data.""",
)

parser.add_argument(
    "-terrE",
    "-times-error-exact",
    nargs="+",
    default=[],
    help="""A list of all times for which error calculations should be performed against the exact values of the function. This list will also be used for visualisation of error data.""",
)

parser.add_argument(
    "-tani",
    "-times-start-animation",
    nargs="+",
    default=[],
    help="""A list of time ranges in the format "[t0,t1]" to be used for the start and stop times of animations.""",
)

parser.add_argument(
    "-tplo",
    "-times-plot",
    nargs="+",
    default=[],
    help="""A list of all times for which plots of the specified data types should be made.""",
)

# Parse data types
parser.add_argument(
    "-dg",
    "-dgTypes",
    nargs="+",
    default=["raw"],
    help="""A list of the data group types for error calculation and visulisation. Currently only implemented for "raw".""",
)

# Collect args and set up defaults
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)

args = parser.parse_args()

# if args.s is None and args.f is None:
#    print "If the -s option is used -f must be set."
#    parser.print_usage()
#    sys.exit(1)

# if args.dg is None:
#    args.dg = ['raw']

if args.tani is not None:
    args.tani = [eval(tani) for tani in args.tani]
else:
    args.tani = []

print("Beginning setup of commands")

# Get directory
path = os.path.dirname(os.path.abspath(sys.argv[0]))


# Do conversions for commands that can take multiple times and dgs.
def bash_string(sarray, prefix):
    rs = ""
    for s in sarray:
        rs = rs + prefix + " '%s' " % s
    return rs.strip()


args.dg = bash_string(args.dg, "-dg")
args.terrN_bash = bash_string(args.terrN, "-t")
args.terrE_bash = bash_string(args.terrE, "-t")

errorNum_run = []
if args.terrN_bash is not "":
    command = os.path.join(path, "error.py")
    errorNum_run += [
        "python -O '%s' %s %s -Lp '2' numer '%s'"
        % (command, args.terrN_bash, args.dg, args.f)
    ]

hdfvis_error_num = []
command = os.path.join(path, "visualise.py")
for time in args.terrN:
    hdfvis_error_num += [
        'python "%s" %s err -t "%s" "%s"' % (command, args.dg, time, args.f)
    ]

errorExa_run = []
if args.terrE_bash is not "":
    command = os.path.join(path, "error.py")
    errorExa_run += [
        "python -O '%s' %s %s -Lp '2' exact '%s'"
        % (command, args.terrE_bash, args.dg, args.f)
    ]

hdfvis_error_exa = []
command = os.path.join(path, "visualise.py")
for time in args.terrE:
    hdfvis_error_exa += [
        'python "%s" %s err -e -t "%s" "%s"' % (command, args.dg, time, args.f)
    ]

hdfvis_plot = []
command = os.path.join(path, "visualise.py")
for time in args.tplo:
    hdfvis_plot += [
        'python "%s" %s plot -t "%s" "%s"' % (command, args.dg, time, args.f)
    ]

hdfvis_ani = []
command = os.path.join(path, "visualise.py")
for tani in args.tani:
    hdfvis_ani += [
        'python %s "%s" ani -t0 "%s" -t1 "%s" "%s"'
        % (command, args.dg, tani[0], tani[1], args.f)
    ]

print("Command set up complete")


if args.s is not None:
    # if args.f is None:
    #    main_run = "python -O %s"%args.s
    # else:
    #    main_run = "python -O %s -f %s"%(args.s,args.f)
    main_run = "python -O %s -f %s" % (args.s, args.f)
    print("MAIN CALCULATION: " + main_run)
    # args.mfile =subprocess.Popen(main_run,\
    #    shell=True,stdout=subprocess.PIPE).communicate()[0]
    # if args.f is None:
    #    args.f = args.mfile
    subprocess.call(main_run, shell=True)

for s in errorNum_run:
    print("ERROR CALCULATION: " + s)
    subprocess.call(s, shell=True)
for s in hdfvis_error_num:
    print("ERROR VISULISATION: " + s)
    subprocess.call(s, shell=True)
for s in errorExa_run:
    print("ERROR CALCULATION: " + s)
    subprocess.call(s, shell=True)
for s in hdfvis_error_exa:
    print("ERROR VISULISATION: " + s)
    subprocess.call(s, shell=True)
for s in hdfvis_plot:
    print("GRAPHING: " + s)
    subprocess.call(s, shell=True)
for s in hdfvis_ani:
    print("ANIMATION: " + s)
    subprocess.call(s, shell=True)
