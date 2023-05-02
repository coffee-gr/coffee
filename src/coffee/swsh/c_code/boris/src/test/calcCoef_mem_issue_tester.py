import ctypes
from ctypes import Structure
import logging
import math
import numpy as np
import os
import inspect
import h5py
from numpy import ctypeslib, typeDict


_cg = ctypes.CDLL(
    os.path.abspath(
        os.path.join(
            os.path.dirname(inspect.getfile(inspect.currentframe())),
            "..",
            "lib",
            "libw3jboris.so.1",
        )
    )
)

_cg.hashRegge.restype = ctypes.c_ulonglong * 3
_cg.hashRegge.argtypes = [
    ctypeslib.ndpointer(
        dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
    ),
    ctypeslib.ndpointer(
        dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
    ),
]

_cg.calcCoeff.restype = ctypes.c_void_p
_cg.calcCoeff.argtypes = [
    ctypeslib.ndpointer(
        dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
    ),
    ctypeslib.ndpointer(
        dtype=typeDict["double"], ndim=1, flags="CONTIGUOUS, WRITEABLE, ALIGNED"
    ),
]


if __name__ == "__main__":
    x = 20
    jms = be.array([x, x, 0, 0, 0], dtype=typeDict["double"])
    j1_max = x + x
    j1_min = max(math.fabs(x - x), math.fabs(0))
    j1_num = j1_max - j1_min + 1
    w3j_arr = be.empty((j1_num,), dtype=typeDict["double"])
    print(w3j_arr)
    print("Entering calcCoeff")
    _cg.calcCoeff(jms, w3j_arr)
    print("Exit calcCoeff")
    print(w3j_arr, j1_min, j1_num, j1_max)
