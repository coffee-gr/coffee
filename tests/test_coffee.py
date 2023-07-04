import h5py
import os

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

for i in range(0, len(backends)):

    # Create the exec string
    exec_string = 'python' + py_version + ' ' + srcd + exd + '/' + \
        setup_file + ' ' + backends[i] + '_' + str(i) + '.hdf ' + backends[i]
    
    print(str(i) + ": Testing " + backends[i] + ' backend.')
    os.system(exec_string)