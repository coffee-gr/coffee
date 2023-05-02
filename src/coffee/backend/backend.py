#!/usr/bin/env python
# encoding: utf-8

"""
Methods that wrap common computations
"""
import numpy as np

################################################################################
# NumPy implementation
################################################################################

class BE:
    """
    A class that provides a generic method for calling NumPy functions.
    """

    def __getattr__(self, name):
        """
        Returns a method that calls the corresponding NumPy method with the same name.
        """
        if hasattr(np, name):
            return lambda *args, **kwargs: getattr(np, name)(*args, **kwargs)
        else:
            raise AttributeError(f"{self.__class__.__name__} object has no attribute '{name}'")

    def lib_stride_tricks_as_strided(self, x, shape=None, strides=None, writeable=True):
        return np.lib.stride_tricks.as_strided(x, shape=shape, strides=strides, writeable=writeable)

be = BE()