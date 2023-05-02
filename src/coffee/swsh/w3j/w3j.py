"""This module contains classes and routines to compute 
Wigner 3j functions.
"""


import ctypes
import logging
import math
import numpy as np
import os
import inspect
import h5py
import subprocess
from numpy import ctypeslib, typeDict
import sys

NUMERICAL_ERROR_TOLERANCE = 1e-15
"""The numerical tolerance used in float comparisons."""


def valid(j1, j2, j3, m1, m2, m3):
    """A check if the ji's and mi's specify a valid 3j symbol.

    For the ji's and mi's to specify a valid 3j symbol it is necessary for
    the ji's to be integer or half-integers and the mi's must be such that
    mi is one of -ji, -ji+1, ..., 0, ..., ji-1, ji.

    We perform the calculation up to a numerical tolerance. This
    tolerance can be specified by setting the module variable
    NUMERICAL_ERROR_TOLERANCE. The default value is 1e-15.

    Parameters
    ----------
    j1 : float
        A float representation of an integer or half-integer.
    j2 : float
        A float representation of an integer or half-integer.
    j3 : float
        A float representation of an integer or half-integer.
    m1 : float
        A float representation of an integer or half-integer.
    m2 : float
        A float representation of an integer or half-integer.
    m3 : float
        A float representation of an integer or half-integer.

    Returns
    -------
    bool:
        If if the parameters correspond to a valid w3j symbol.
    """
    # are the ji's integers or half-integers?
    if not (
        2 * j1 % 1 < NUMERICAL_ERROR_TOLERANCE
        and 2 * j2 % 1 < NUMERICAL_ERROR_TOLERANCE
        and 2 * j3 % 1 < NUMERICAL_ERROR_TOLERANCE
    ):
        return False
    # Are the ms in the correct arange?
    js = be.array(
        [
            be.arange(-j1, j1 + 1, 1),
            be.arange(-j2, j2 + 1, 1),
            be.arange(-j3, j3 + 1, 1),
        ]
    )
    ms = be.array([m1, m2, m3])
    for i in range(0, 3):
        js[i] = be.absolute(js[i] - ms[i])
    if not (
        any(js[0] < NUMERICAL_ERROR_TOLERANCE)
        and any(js[1] < NUMERICAL_ERROR_TOLERANCE)
        and any(js[2] < NUMERICAL_ERROR_TOLERANCE)
    ):
        return False
    return True


def trivial_zero(j1, j2, j3, m1, m2, m3):
    """A check if the ji's and mi's specify a trival zero.

    Parameters
    ----------
    j1 : float
        A float representation of an integer or half-integer.
    j2 : float
        A float representation of an integer or half-integer.
    j3 : float
        A float representation of an integer or half-integer.
    m1 : float
        A float representation of an integer or half-integer.
    m2 : float
        A float representation of an integer or half-integer.
    m3 : float
        A float representation of an integer or half-integer.

    Returns
    -------
    bool:
        If if the parameters describe a trivial zero w3j symbol.
    """
    WI = j1 + j2 + j3
    if not WI >= 0 and not WI % 1 < NUMERICAL_ERROR_TOLERANCE:
        return True
    WII = abs(m1 + m2 + m3)
    if not WII < NUMERICAL_ERROR_TOLERANCE:
        return True
    if not math.fabs(j1 - j2) <= j3 <= j1 + j2:
        if not math.fabs(j1 - j3) <= j2 <= j1 + j3:
            if not math.fabs(j3 - j2) <= j1 <= j3 + j2:
                return True
    return False


class W3jStone(object):
    """Provides python access to Anthony Stone's RRF code (v 3.2), via the program
    rrfcalc, allowing computation of w3j symbols.

    The program rrfcalc must be installed for this class to work.

    For some reason wrapping his fortran code directly with f2py
    ends with errors that are difficult to diagnose.
    """

    def __call__(self, j1, j2, j3, m1, m2, m3):
        """Returns the double precision value of the exact calculation of the
        3j symbol:

            (j1 j2 j3)
            (m1 m2 m3)

        Parameters
        ----------
        j1 : float
            A float representation of an integer or half-integer.
        j2 : float
            A float representation of an integer or half-integer.
        j3 : float
            A float representation of an integer or half-integer.
        m1 : float
            A float representation of an integer or half-integer,
            m1 is one of -j1, -j1+1, ..., 0, ..., j1-1, j1.
        m2 : float
            A float representation of an integer or half-integer,
            m2 is one of -j2, -j2+1, ..., 0, ..., j2-1, j2.
        m3 : float
            A float representation of an integer or half-integer,
            m3 is one of -j3, -j3+1, ..., 0, ..., j3-1, j3.

        Returns
        -------
        float:
            The value of the 3j symbol,

                (j1, j2, j3)
                (m1, m2, m3)
        """
        if not valid(j1, j2, j3, m1, m2, m3):
            raise ValueError(
                "j1, j2, j3 must be integers or half integers the \
            mi's must be so that mi is one of -ji, -ji+1, ..., 0, ..., ji-1, ji \
            where i is 1, 2 or 3."
            )
        if trivial_zero(j1, j2, j3, m1, m2, m3):
            return 0.0
        rrfcalc = os.path.abspath(
            os.path.join(
                os.path.dirname(inspect.getfile(inspect.currentframe())),
                "..",
                "lib",
                "rrfcalc",
            )
        )
        exact = subprocess.Popen(
            [rrfcalc], stdin=subprocess.PIPE, stdout=subprocess.PIPE
        )
        argument = "3j %d/2 %d/2 %d/2 %d/2 %d/2 %d/2" % (
            2 * j1,
            2 * j2,
            2 * j3,
            2 * m1,
            2 * m2,
            2 * m3,
        )
        stdout, stderr = exact.communicate(argument)
        rv = eval(stdout[147:-4].replace("sqrt", "math.sqrt"))
        return rv

    def exact(self, j1, j2, j3, m1, m2, m3):
        """Returns the exact value of the
        3j symbol

            (j1 j2 j3)
            (m1 m2 m3)

        as a string.

        The format of the string is
        <+/-><integer>/<integer>*math.sqrt(<integer>/<integer>). The double
        precision value represented by this string can be retrieved by using the
        eval() builtin.

        Parameters
        ----------
        j1 : float
            A float representation of an integer or half-integer.
        j2 : float
            A float representation of an integer or half-integer.
        j3 : float
            A float representation of an integer or half-integer.
        m1 : float
            A float representation of an integer or half-integer,
            m1 is one of -j1, -j1+1, ..., 0, ..., j1-1, j1.
        m2 : float
            A float representation of an integer or half-integer,
            m2 is one of -j2, -j2+1, ..., 0, ..., j2-1, j2.
        m3 : float
            A float representation of an integer or half-integer,
            m3 is one of -j3, -j3+1, ..., 0, ..., j3-1, j3.

        Returns
        -------
        string:
            Returns a string representation of the value of the 3j symbol:

                (j1, j2, j3)
                (m1, m2, m3)

            as a string. The format of the string is
            <+/-><integer>/<integer>*math.sqrt(<integer>/<integer>). The double
            precison value represented by this string can be retrieved by using the
            eval() builtin.
        """
        if not valid(j1, j2, j3, m1, m2, m3):
            raise ValueError(
                "j1, j2, j3 must be integers or half integers the \
            mi's must be so that mi is one of -ji, -ji+1, ..., 0, ..., ji-1, ji \
            where i is 1, 2 or 3."
            )
        if trivial_zero(j1, j2, j3, m1, m2, m3):
            return 0.0
        rrfcalc = os.path.abspath(
            os.path.join(
                os.path.dirname(inspect.getfile(inspect.currentframe())),
                "..",
                "lib",
                "rrfcalc",
            )
        )
        exact = subprocess.Popen(
            [rrfcalc], stdin=subprocess.PIPE, stdout=subprocess.PIPE
        )
        argument = "3j %d/2 %d/2 %d/2 %d/2 %d/2 %d/2" % (
            2 * j1,
            2 * j2,
            2 * j3,
            2 * m1,
            2 * m2,
            2 * m3,
        )
        stdout = exact.communicate(argument)[0][147:-4].replace("sqrt", "math.sqrt")
        return stdout


class W3jBen(object):
    """Provides a wrapper to 'lower level' functions in Boris' code.

    As a consequence this class offers some capabilities that W3jBoris does not.
    Please note, however, that there maybe some degradation in performance as a
    result. This class is recommended when it is necessary to conserve memory.

    Specifically this class only calls hashRegge of mathWigner3j.c and
    calcCoef of mathWigner3jRecursion.c. The class uses an hdf file to store
    the values of 3j symbols. When a 3j symbol is requested the file is checked
    to see if the value is stored. If the value is not stored then the recursive
    algorithm of coefCalc (in mathWigner3jRecursion.c) is used to calculate it.
    This has the side effect that a number of other 3j symbols are also
    calculated. All calculated symbols are stored in the hdf file.

    There is currently no caching facility implemented.

    Boris' code can be found in the folder c_code/boris. It must be compiled
    as a shared library 'libboris.so.1.0.1' in the folder lib/.
    """

    def __init__(self, hdf_file):
        """Construct a W3jBen object.

        Parameters
        ----------
        hdf_file : string
            The name of the hdf_file which stores, or will be used to
            store the values of needed 3j symbols.
        """
        # Set up api for boris' code
        self._cg = ctypes.CDLL(
            os.path.abspath(
                os.path.join(
                    os.path.dirname(inspect.getfile(inspect.currentframe())),
                    "..",
                    "lib",
                    "libboris.so.1.0.1",
                )
            )
        )

        def _hashRegge_check(result, func, arguments):
            a = be.fromiter(result, dtype=typeDict["ulonglong"], count=3)
            return a

        self._cg.hashRegge.restype = ctypes.POINTER(ctypes.c_ulonglong)
        self._cg.hashRegge.errcheck = _hashRegge_check
        self._cg.hashRegge.argtypes = [
            ctypeslib.ndpointer(
                dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
            ),
            ctypeslib.ndpointer(
                dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
            ),
        ]

        self._cg.calcCoeff.restype = ctypes.c_void_p
        self._cg.calcCoeff.argtypes = [
            ctypeslib.ndpointer(
                dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
            ),
            ctypeslib.ndpointer(
                dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
            ),
        ]

        # Initialise internal variables
        self.hdf_file = h5py.File(hdf_file)
        self.log = logging.getLogger("W3jBen")

    def __del__(self):
        """Ensures that the hdf file is closed when this object is garbage \
        collected."""
        self.hdf_file.close()

    def __call__(self, j1, j2, j3, m1, m2, m3):
        """Returns the value of a 3j symbol.

        Checks if the value of the 3j symbol is stored in the file. If so the
        value is returned. If not then the value is computed, along with a
        number of other 3j symbol values. All values are stored (if they are
        not already present in the file) and the needed value is returned.

        Parameters
        ----------
        j1 : float
            A float representation of an integer or half-integer.
        j2 : float
            A float representation of an integer or half-integer.
        j3 : float
            A float representation of an integer or half-integer.
        m1 : float
            A float representation of an integer or half-integer,
            m1 is one of -j1, -j1+1, ..., 0, ..., j1-1, j1.
        m2 : float
            A float representation of an integer or half-integer,
            m2 is one of -j2, -j2+1, ..., 0, ..., j2-1, j2.
        m3 : float
            A float representation of an integer or half-integer,
            m3 is one of -j3, -j3+1, ..., 0, ..., j3-1, j3.

        Returns
        -------
        float:
            The value of the 3j symbol,

                (j1, j2, j3)
                (m1, m2, m3)
        """
        if not valid(j1, j2, j3, m1, m2, m3):
            raise ValueError(
                "j1, j2, j3 must be integers or half integers the \
            mi's must be so that mi is one of -ji, -ji+1, ..., 0, ..., ji-1, ji \
            where i is 1, 2 or 3."
            )
        if trivial_zero(j1, j2, j3, m1, m2, m3):
            return 0.0
        js = be.array([j1, j2, j3], dtype=typeDict["double"])
        ms = be.array([m1, m2, m3], dtype=typeDict["double"])
        i = self._cg.hashRegge(js, ms)
        regge_hash = i[0]
        odd_flip = i[1]
        J = i[2]
        data_set = self.hdf_file.get(repr(regge_hash))
        if data_set is not None:
            rv = data_set[0]
            if odd_flip and J % 2 == 1:
                rv = -rv
        else:
            w3j_arr, j1_min, j1_num, j1_max = self._compute_values(j2, j3, m1, m2, m3)
            self._write_w3j_symbols(w3j_arr, j1_min, j1_num, j1_max, j2, j3, m1, m2, m3)
            index = int(j1 - j1_min)
            rv = w3j_arr[index]
        return rv

    def _write_w3j_symbols(
        self, w3j_symbols, j1_min, j1_num, j1_max, j2, j3, m1, m2, m3
    ):
        """Writes an array of 3j symbol values to the hdf file, skipping those values
        that already exist.

        Parameters
        ----------
        w3j_symbols : list of floats
            Each element of the list is the value of a w3j symbol indexed with
            respect to the `j1` entry.
        j1_min : float
            A float representation of an integer or half-integer, the minimum
            value for the `j1` entry.
        j1_num : float
            A float representation of an integer or half-integer, the number of
            entries in the the w3j_symbols list.
        j1_max : float
            A float representation of an integer or half-integer, the maximum
            value for the `j1`.
        j2 : float
            A float representation of an integer or half-integer.
        j3 : float
            A float representation of an integer or half-integer.
        m1 : float
            A float representation of an integer or half-integer,
            m1 is one of -j1, -j1+1, ..., 0, ..., j1-1, j1.
        m2 : float
            A float representation of an integer or half-integer,
            m2 is one of -j2, -j2+1, ..., 0, ..., j2-1, j2.
        m3 : float
            A float representation of an integer or half-integer,
            m3 is one of -j3, -j3+1, ..., 0, ..., j3-1, j3.

        """
        ms = be.array([m1, m2, m3], dtype=typeDict["double"])
        js = be.array([j1_min, j2, j3], dtype=typeDict["double"])
        for i in range(0, j1_num):
            array = self._cg.hashRegge(js, ms)
            regge_hash = array[0]
            data_set = self.hdf_file.get(repr(regge_hash))
            if data_set is None:
                odd_flip = array[1]
                J = array[2]
                data = w3j_symbols[i]
                if odd_flip and J % 2 == 1:
                    data = -data
                data_set = self.hdf_file.create_dataset(repr(regge_hash), data=[data])
            js[0] = js[0] + 1

    def _compute_values(self, j2, j3, m1, m2, m3):
        """Computes a range of 3j symbols.

        Parameters
        ----------
        j2 : float
            A float representation of an integer or half-integer.
        j3 : float
            A float representation of an integer or half-integer.
        m1 : float
            A float representation of an integer or half-integer,
            m1 is one of -j1, -j1+1, ..., 0, ..., j1-1, j1.
        m2 : float
            A float representation of an integer or half-integer,
            m2 is one of -j2, -j2+1, ..., 0, ..., j2-1, j2.
        m3 : float
            A float representation of an integer or half-integer,
            m3 is one of -j3, -j3+1, ..., 0, ..., j3-1, j3.

        Returns
        -------
        a tuple : (list of floats, float, int, float)
            The list of floats contains the values of w3j symbols, the index
            of each entry corresponds to the `j1` entry. The first float
            is the minimum `j1` entry, the second float is the maximum
            `j1` entry. The int is the number of `j1` entries, that is the
            length of the list.
        """
        jms = be.array([j2, j3, m1, m2, m3], dtype=typeDict["double"])
        j1_max = j2 + j3
        j1_min = max(math.fabs(j2 - j3), math.fabs(m1))
        j1_num = int((j1_max - j1_min + 1) // 1)
        w3j_arr = be.empty((j1_num,), dtype=typeDict["double"])
        self._cg.calcCoeff(jms, w3j_arr)
        return w3j_arr, j1_min, j1_num, j1_max


class W3jBoris(object):
    """This class provides a wrapper to the portions of Boris' code that
    compute w3j symbols.

    Each object is passed a bandwidth limitor, L, and a step size (either 1 or
    0.5). The data giving the values of all 3j symbols with respect to L and dl
    can then be either loaded from the file data.h5 or calculated. The values
    are then accessed by calling the object. For example to get the value
    of the 3j symbol:

        ( j1, j2, j3 )
        ( m1, m2, m3 )

    from object, w, one uses the code:

        w(j1, j2, j3, m1, m2, m3).

    The name of the file is currently hard coded into Boris' code and cannot be
    changed.

    Calculated values can be written to the file data.h5 by using the method
    save().

    If at some later data point in time it is necessary to increase the
    bandwidth limit, this object requires recalculation of all w3j values and
    storing this information separately in the data.h5 file.

    As a consequence some thought is required before hand to ensure that the
    file does not become bloated with multiple copies of the same data.

    Boris' code can be found in the folder c_code/boris. It must be compiled
    as a shard library 'libboris.so.1.0.1' in the folder lib/.
    """

    def __init__(self, L, dl):
        """Construct a w3j_boris object.

        Parameters
        ----------
        L : float
            The bandwidth limitor. j2, j3, m1 and m2 are all restricted to be
            less than or equal to L. The value m3 is always taken as -(m1+m2)
            and, therefore, is such that -2L <= m3 <= 0. The value j1 is such
            that `max( |j2-j3|, |m1| ) <= j1 <= j2 + j3`.
        dl : float
            If only integer values of j1, j2, j3, m1, m2 and m3 are required
            then set equal to 1. If half-integer values are required then set
            equal to 0.5.
        """
        # Set up api for boris' code
        self._cg = ctypes.CDLL(
            os.path.abspath(
                os.path.join(
                    os.path.dirname(inspect.getfile(inspect.currentframe())),
                    "..",
                    "lib",
                    "libboris.so.1.0.1",
                )
            )
        )

        self._cg.calcw3jSymbols.restype = ctypes.c_void_p
        self._cg.calcw3jSymbols.argtypes = [ctypes.c_double, ctypes.c_double]

        self._cg.getW3jPrecalc.restype = ctypes.c_double
        self._cg.getW3jPrecalc.argtypes = [
            ctypeslib.ndpointer(
                dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
            ),
            ctypeslib.ndpointer(
                dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
            ),
            ctypes.c_void_p,
        ]

        self._cg.writeCw3jStruct.restype = ctypes.c_void_p
        self._cg.writeCw3jStruct.argtypes = [ctypes.c_void_p]

        self._cg.checkCw3jStructExists.restype = ctypes.c_int
        self._cg.checkCw3jStructExists.argtypes = [ctypes.c_double, ctypes.c_double]

        self._cg.readcplCw3jStruct.restype = ctypes.c_void_p
        self._cg.readcplCw3jStruct.argtypes = [ctypes.c_double, ctypes.c_double]

        # Initialise internal variables
        self.L = L
        self.dl = dl
        self.log = logging.getLogger("W3jBoris")
        self._cplCw3jStruct_ptr = None

    def load(self):
        """Loads precalculated data from the file data.h5.

        The name of this file is hard coded into Boris' code.
        """
        self._cg.hdf5InterfaceInitialise()
        if not self._cg.checkCw3jStructExists(self.L, self.dl):
            raise Exception("Data.h5 doesnot contain the requested data")
        self._cplCw3jStruct_ptr = self._cg.readcplCw3jStruct(self.L, self.dl)

    def save(self):
        """Save the calculated w3j symbols to the file data.h5 using Boris'
        data format.

        The name of this file is hard coded.
        """
        self._cg.hdf5InterfaceInitialise()
        if self._cplCw3jStruct_ptr is None:
            raise Exception(
                "The values cannot be saved as they have not been calculated"
            )
        if self._cg.checkCw3jStructExists(self.L, self.dl):
            return
        self._cg.writeCw3jStruct(self._cplCw3jStruct_ptr)

    def calculate(self):
        """Calculates all w3j symbols below the bandwidth limit L."""
        self.log.info("Performing calculation of w3j symbols")
        self._cplCw3jStruct_ptr = self._cg.calcw3jSymbols(self.L, self.dl)
        self.log.info("Calculation of w3j symbols completed")

    def __call__(self, j1, j2, j3, m1, m2, m3):
        """Returns the value of a 3j symbol.

        Parameters
        ----------
        j1 : float
            A float representation of an integer or half-integer.
        j2 : float
            A float representation of an integer or half-integer.
        j3 : float
            A float representation of an integer or half-integer.
        m1 : float
            A float representation of an integer or half-integer,
            m1 is one of -j1, -j1+1, ..., 0, ..., j1-1, j1.
        m2 : float
            A float representation of an integer or half-integer,
            m2 is one of -j2, -j2+1, ..., 0, ..., j2-1, j2.
        m3 : float
            A float representation of an integer or half-integer,
            m3 is one of -j3, -j3+1, ..., 0, ..., j3-1, j3.

        Returns
        -------
        float:
            The value of the 3j symbol,

                (j1, j2, j3)
                (m1, m2, m3)
        """
        if not valid(j1, j2, j3, m1, m2, m3):
            raise ValueError(
                "j1, j2, j3 must be integers or half integers the \
            mi's must be so that mi is one of-ji, -ji+1, ..., 0, ..., ji-1, ji \
            where i is 1, 2 or 3."
            )
        if trivial_zero(j1, j2, j3, m1, m2, m3):
            return 0.0
        if self._cplCw3jStruct_ptr is None:
            raise Exception("Precalculation of the symbols has not been done")
        js = be.array([j1, j2, j3], dtype=typeDict["double"])
        ms = be.array([m1, m2, m3], dtype=typeDict["double"])
        rv = self._cg.getW3jPrecalc(js, ms, self._cplCw3jStruct_ptr)
        return rv


# if __name__ == "__main__":
# wBen = W3jBen("test.hdf")
# print wBen(0, 1, 1, 0.0, 0, 0)
# sys.exit(0)
# x = 2
# wStone = W3jStone()
# wBen = W3jBen("test.hdf")
# wBoris = W3jBoris(x, 0.5)
# wBoris.calculate()
# ms = be.array([(u,v,w)
# for u in be.arange(-x,x,0.5)
# for v in be.arange(-x,x,0.5)
# for w in be.arange(-x,x,0.5)])
# js = be.array([(u,v)
# for u in be.arange(0,x,0.5)
# for v in be.arange(0,x,0.5)])
# print "Objects created"
# print "Test commensing"
# count_range = range(0, 16 * x + 1, 2 * x)
# for j1 in be.arange(0,x,0.5):
# print "\nj1 is %f"%j1
# for j2, j3 in js:
# for m1, m2, m3 in ms:
##sys.stdout.write('.')
# if valid(j1, j2, j3, m1, m2, m3):
##print "collecting 3j values"
# a = wStone(j1, j2, j3, m1, m2, m3)
# b = wBen(j1, j2, j3, m1, m2, m3)
# c = wBoris(j1, j2, j3, m1, m2, m3)
# if (abs(a-b)< NUMERICAL_ERROR_TOLERANCE
# and abs(b-c) < NUMERICAL_ERROR_TOLERANCE
# and abs(a-c) < NUMERICAL_ERROR_TOLERANCE):
# continue
# else:
# print "\nFail: %f, %f, %f, %f, %f, %f"%(j1, j2, j3, m1, m2, m3)
# print "abs(wStone - wBen) = %.15f"%(abs(a-b))
# print "abs(wStone - wBoris) = %.15f"%(abs(a-c))
# print "abs(wBen - wBoris) = %.15f"%(abs(b-c))
# sys.stdout.flush()
# print "\nTest completed"
