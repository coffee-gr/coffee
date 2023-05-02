import argparse
import sys
from builtins import range

import backend as be

from coffee.io import simulation_data as sd

be.set_printoptions(threshold=be.inf)

parser = argparse.ArgumentParser(
    description="""This program compares the contents of two hdf files as produced by coffee."""
)

parser.add_argument("-f1", help="""The first file to be compared.""")
parser.add_argument("-f2", help="""The second file to be compared.""")
parser.add_argument(
    "-dg",
    "-dgtype",
    help="""A data group to be analysed. See simulation_data for information on
the different data groups. Defaults to "raw". Multiple data groups can be given.""",
)

args = parser.parse_args()
if args.dg is None:
    args.dg = "raw"

first = sd.SimulationHDF(args.f1)
second = sd.SimulationHDF(args.f2)

first_sims = first.getSims()
second_sims = second.getSims()

if first_sims != second_sims:
    print("Different simulations in files")
    sys.exit()


def get_sim(sims, name):
    for sim in sims:
        if sim.name == name:
            return sim
    raise ValueError("No sim named " + name)


for f_sim in first_sims:
    s_sim = get_sim(second_sims, f_sim.name)
    f_datagroup = f_sim.getDgType(args.dg)
    s_datagroup = s_sim.getDgType(args.dg)

    if len(f_datagroup) != len(s_datagroup):
        print(
            "The lengths of %s datagroup in simulations %s and %s do not match"
            % (args.dg, f_sim.name, s_sim.name)
        )
        sys.exit()

    index = 0
    while index < len(f_datagroup):
        f_data = f_datagroup[index][()]
        s_data = s_datagroup[index][()]
        if not be.array_equal(f_data, s_data):
            for i in range(len(f_data)):
                if not be.array_equal(f_data[i], s_data[i]):
                    print("Difference in component %s" % repr(i))
                    for j in range(len(f_data[i])):
                        if f_data[i][j] != s_data[i][j]:
                            print("Difference in element %s" % repr(j))
                            print("%s != %s" % (repr(f_data[i][j]), repr(s_data[i][j])))
                    print(
                        "Different data encountered at index %i:\n%s\n%s"
                        % (index, repr(f_data[i]), repr(s_data[i]))
                    )
            sys.exit()
        index += 1

print("DataType %s in files %s and %s is the same" % (args.dg, args.f1, args.f2))
