import sys

sys.path.append("../../EvolutionSBP/")
import simulation_data

# code that helps to read the content of an hdf file
# sys.argv[1]: give the file.
# sys.argv[2]: give the grid points used in the simulation, i.e. 0 corresponds to 50 points,,1 corresponds to 100 points, etc...
# sys.argv[3]: give timeslice.
# sys.argv[4]: give component of Weyl.
# sys.argv[5]: give gridpoint.

with simulation_data.SimulationHDF(sys.argv[1]) as file:
    sims = file.getSims()
    run = int(sys.argv[2])  # number of grid point used in the simulation.
    timeslice = int(sys.argv[3])
    comp = int(sys.argv[4])
    gridp = int(sys.argv[5])

    print("Raw value = %f" % sims[run].raw[timeslice][comp][gridp])
    print(
        "Transformed value = %f"
        % (
            (sims[run].domain[timeslice][gridp] + 1) ** (3)
            * sims[run].raw[timeslice][comp][gridp]
        )
    )
    print(sims[run].domain[timeslice][gridp])
