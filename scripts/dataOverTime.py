from builtins import input
import sys

sys.path.append("../../EvolutionSBP/")
import numpy as np
import h5py_array
import Gnuplot
import simulation_data as sd

# sys.argv[1] = 'exp_-20xx.hdf' - hdf data file
# sys.argv[2] = raw or weyl - dgType
# sys.argv[2] = 'gif' - output type (optional)
# sys.argv[3] = 'test.gif' - output file name (optional)
# sys.argv[4] = derived attr name

with sd.SimulationHDF(sys.argv[1]) as file:
    print("Producing graph...")
    # Get simulations
    sims = file.getSims()

    # Get domain for comparison of errors
    domain = sims[0].domain[0]

    # Set up Gnuplot
    g = Gnuplot.Gnuplot()
    g("set style data lines")
    g('set title "L2 of I and J from %s" enhanced' % (sys.argv[1],))
    g('set xlabel "time"')
    g('set ylabel "Log base 2 of sum of absolute value of I and J"')
    g("set key left")
    # g('set yrange [0:1]')
    # g('set xrange [0:1]')
    wait = True
    try:
        g("set terminal " + sys.argv[3])
        g('set out "' + sys.argv[4] + '"')
        wait = False
    except:
        pass

    # check if using derived attr

    # Collate data since we are graphing with respect to time
    # and want to compare numerical values of things we must make sure
    # to use the same domain
    plot_data = []
    for sim in sims:
        print("Calculating %s" % sim.name)
        data = []
        time = []
        mapping = sd.array_value_index_mapping(domain, sim.domain[0])
        for i, ds in enumerate(sim.getDgType(sys.argv[2])):
            dl = be.zeros((2,))
            for map in mapping:
                if sim.domain[0][map[1]] <= 0:
                    dl = dl + be.absolute(ds.value[:, map[1]])
            data += [be.sum(dl)]
            time += [sim.time[i][0]]
        plot_data += [Gnuplot.Data(time, be.log2(data), title=sim.name)]

    # Plot data
    g.plot(*plot_data)
    if wait:
        eval(input("Press key to finish programme"))
    print("...Done.-")
