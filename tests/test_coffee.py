import sys
import os
import numpy as np

from coffee.settings import init
from coffee.backend import backend as be
my_backend = be.set_backend('numpy')
init(my_backend)

from coffee.io import simulation_data as sd

# Test with OneDAdvectionMpi_sbp
srcd = os.path.dirname(os.path.realpath(__file__)) + "/"
cwd  = os.getcwd()

# Backends to test
backends = ["numpy", "numpy"]

# Python version
py_version = '3.11'

# Example to test with
exd        = 'TwoDAdvectionMpi_sbp'
setup_file = 'TwoDAdvection_setup.py'

# # Run the examples
for i in range(0, len(backends)):

    # Create the exec string
    exec_string = 'python' + py_version + ' ' + srcd + exd + '/' + \
        setup_file + ' ' + backends[i] + '_' + str(i) + '.hdf ' + backends[i]
    
    print(str(i) + ": Testing " + backends[i] + ' backend.')
    os.system(exec_string)

# Compare
exact_filename = backends[0] + '_' + str(0) + '.hdf'
exact = sd.SimulationHDF(exact_filename)
for i in range(1, len(backends)):
    test_filename = backends[i] + '_' + str(i) + '.hdf'
    test = sd.SimulationHDF(test_filename)

    exact_sim = exact.getSims()
    test_sim  = test.getSims()

    if exact_sim != test_sim:
        print("Different simulations in files")
        sys.exit()

    def get_sim(sims, name):
        for sim in sims:
            if sim.name == name:
                return sim
        raise ValueError("No sim named " + name)

    for e_sim in exact_sim:
        t_sim = get_sim(exact_sim, e_sim.name)
        e_datagroup = e_sim.getDgType('raw')
        t_datagroup = t_sim.getDgType('raw')
        
        if len(e_datagroup) != len(t_datagroup):
            print("The lengths of %s datagroup in simulations %s and %s do not match"%\
            ('raw', e_sim.name, t_sim.name))
            sys.exit()

        rtol = 0.
        atol = 1e-14

        index = 0
        while index < len(e_datagroup):
            e_data = e_datagroup[index][()]
            t_data = t_datagroup[index][()]

            if not np.allclose(e_data, t_data.real, rtol=rtol, atol=atol):
                for i in range(len(e_data)):
                    if not np.allclose(e_data[i], t_data[i].real, rtol=rtol, atol=atol):
                        print("Difference in x-component %s"%repr(i))
                        for k in range(len(e_data)):
                            if not np.allclose(e_data[i][k], t_data[i][k].real, rtol=rtol, atol=atol):
                                print("Difference in y-component %s"%repr(i))
                                for j in range(len(f_data[i][k])):
                                    if e_data[i][k][j] != t_data[i][k][j].real:
                                        print("Difference in element %s"%repr(j))
                                        print("%s != %s"%(repr(e_data[i][k][j]), repr(t_data[i][k][j].real)))
                                print("Different data encountered at index %i:\n%s\n%s"%\
                                (index, repr(e_data[i]), repr(t_data[i].real)))
                sys.exit()
            index += 1

    print("DataType %s in files %s and %s is the same up to a tolerance of %.1e"%\
    ('raw', exact_filename, test_filename, atol))