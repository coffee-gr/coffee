#!/usr/bin/env python
# encoding: utf-8

"""
Methods that wrap common computations.
"""

class Backend:
    def __init__(self):
        self.backend_name = None
        self.backend = None

    def set_backend(self, backend_name="numpy"):
        self.backend_name = backend_name
        self.backend = self.get_backend(backend_name)

        # Set constants
        self.pi = self.backend.pi

        # Return backend instance
        return self.backend

    def get_backend(self, backend_name):
        """
        Returns an instance of the specified backend class.

        This function creates and returns an instance of the backend class specified by the
        `backend_name` argument. Currently supported backends are "numpy" and others can be
        added as needed.

        Args:
            backend_name (str, optional): The name of the backend for which an instance is to be created.
                Defaults to "numpy".

        Raises:
            ValueError: If the specified `backend_name` is not supported.

        Returns:
            object: An instance of the specified backend class.
        """
        backends = {
            "numpy": NumpyBackend,
            # Add other backends here as needed
        }

        if backend_name not in backends:
            raise ValueError(f"Unsupported backend: {backend_name}")

        return backends[backend_name]()

    def abs(self, *args, **kwargs):
        return self.backend.abs(*args, **kwargs)

    def absolute(self, *args, **kwargs):
        return self.backend.absolute(*args, **kwargs)

    def any(self, *args, **kwargs):
        return self.backend.any(*args, **kwargs)

    def apply_along_axis(self, *args, **kwargs):
        return self.backend.apply_along_axis(*args, **kwargs)

    def apply_over_axes(self, *args, **kwargs):
        return self.backend.apply_over_axes(*args, **kwargs)

    def arange(self, *args, **kwargs):
        return self.backend.arange(*args, **kwargs)

    def array(self, *args, **kwargs):
        return self.backend.array(*args, **kwargs)

    def array_equal(self, *args, **kwargs):
        return self.backend.array_equal(*args, **kwargs)

    def asarray(self, *args, **kwargs):
        return self.backend.asarray(*args, **kwargs)

    def atleast_1d(self, *args, **kwargs):
        return self.backend.atleast_1d(*args, **kwargs)

    def atleast_2d(self, *args, **kwargs):
        return self.backend.atleast_2d(*args, **kwargs)

    def convolve(self, *args, **kwargs):
        return self.backend.convolve(*args, **kwargs)
    
    def cos(self, *args, **kwargs):
        return self.backend.cos(*args, **kwargs)

    def diag(self, *args, **kwargs):
        return self.backend.diag(*args, **kwargs)

    def dot(self, *args, **kwargs):
        return self.backend.dot(*args, **kwargs)

    def dtype(self, *args, **kwargs):
        return self.backend.dtype(*args, **kwargs)

    def empty(self, *args, **kwargs):
        return self.backend.empty(*args, **kwargs)

    def empty_like(self, *args, **kwargs):
        return self.backend.empty_like(*args, **kwargs)

    def exp(self, *args, **kwargs):
        return self.backend.exp(*args, **kwargs)

    def fromiter(self, *args, **kwargs):
        return self.backend.fromiter(*args, **kwargs)

    def lib_stride_tricks_as_strided(self, *args, **kwargs):
        return self.backend.lib_stride_tricks_as_strided(*args, **kwargs)

    def linspace(self, *args, **kwargs):
        return self.backend.linspace(*args, **kwargs)

    def log(self, *args, **kwargs):
        return self.backend.log(*args, **kwargs)

    def log2(self, *args, **kwargs):
        return self.backend.log2(*args, **kwargs)

    def mat(self, *args, **kwargs):
        return self.backend.mat(*args, **kwargs)

    def max(self, *args, **kwargs):
        return self.backend.max(*args, **kwargs)

    def min(self, *args, **kwargs):
        return self.backend.min(*args, **kwargs)

    def ndarray(self, *args, **kwargs):
        return self.backend.ndarray(*args, **kwargs)

    def nonzero(self, *args, **kwargs):
        return self.backend.nonzero(*args, **kwargs)

    def ones(self, *args, **kwargs):
        return self.backend.ones(*args, **kwargs)

    def ones_like(self, *args, **kwargs):
        return self.backend.ones_like(*args, **kwargs)

    def power(self, *args, **kwargs):
        return self.backend.power(*args, **kwargs)

    def savetxt(self, *args, **kwargs):
        return self.backend.savetxt(*args, **kwargs)

    def set_printoptions(self, *args, **kwargs):
        return self.backend.set_printoptions(*args, **kwargs)
    
    def sin(self, *args, **kwargs):
        return self.backend.sin(*args, **kwargs)

    def squeeze(self, *args, **kwargs):
        return self.backend.squeeze(*args, **kwargs)

    def sum(self, *args, **kwargs):
        return self.backend.sum(*args, **kwargs)

    def tan(self, *args, **kwargs):
        return self.backend.tan(*args, **kwargs)

    def union1d(self, *args, **kwargs):
        return self.backend.union1d(*args, **kwargs)

    def vectorize(self, *args, **kwargs):
        return self.backend.vectorize(*args, **kwargs)

    def where(self, *args, **kwargs):
        return self.backend.where(*args, **kwargs)

    def zeros(self, *args, **kwargs):
        return self.backend.zeros(*args, **kwargs)

    def zeros_like(self, *args, **kwargs):
        return self.backend.zeros_like(*args, **kwargs)


class BackendBase:
    def __init__(self):
        pass


# ################################################################################
# # NumPy implementation
# ################################################################################


class NumpyBackend(BackendBase):
    def __init__(self):
        # print("Initialising Numpy Backend")
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "NumPy library is not installed. Please install NumPy to use this function."
            )

        self.np = np

        # Constants
        self.pi = np.pi

    def abs(self, *args, **kwargs):
        return self.np.abs(*args, **kwargs)

    def absolute(self, *args, **kwargs):
        return self.np.absolute(*args, **kwargs)

    def any(self, *args, **kwargs):
        return self.np.any(*args, **kwargs)

    def apply_along_axis(self, *args, **kwargs):
        return self.np.apply_along_axis(*args, **kwargs)

    def apply_over_axes(self, *args, **kwargs):
        return self.np.apply_over_axes(*args, **kwargs)

    def arange(self, *args, **kwargs):
        return self.np.arange(*args, **kwargs)

    def array(self, *args, **kwargs):
        return self.np.array(*args, **kwargs)

    def array_equal(self, *args, **kwargs):
        return self.np.array_equal(*args, **kwargs)

    def asarray(self, *args, **kwargs):
        return self.np.asarray(*args, **kwargs)

    def atleast_1d(self, *args, **kwargs):
        return self.np.atleast_1d(*args, **kwargs)

    def atleast_2d(self, *args, **kwargs):
        return self.np.atleast_2d(*args, **kwargs)

    def convolve(self, *args, **kwargs):
        return self.np.convolve(*args, **kwargs)
    
    def cos(self, *args, **kwargs):
        return self.np.cos(*args, **kwargs)

    def diag(self, *args, **kwargs):
        return self.np.diag(*args, **kwargs)

    def dot(self, *args, **kwargs):
        return self.np.dot(*args, **kwargs)

    def dtype(self, *args, **kwargs):
        return self.np.dtype(*args, **kwargs)

    def empty(self, *args, **kwargs):
        return self.np.empty(*args, **kwargs)

    def empty_like(self, *args, **kwargs):
        return self.np.empty_like(*args, **kwargs)

    def exp(self, *args, **kwargs):
        return self.np.exp(*args, **kwargs)

    def fromiter(self, *args, **kwargs):
        return self.np.fromiter(*args, **kwargs)

    def lib_stride_tricks_as_strided(self, *args, **kwargs):
        return self.np.lib_stride_tricks_as_strided(*args, **kwargs)

    def linspace(self, *args, **kwargs):
        return self.np.linspace(*args, **kwargs)

    def log(self, *args, **kwargs):
        return self.np.log(*args, **kwargs)

    def log2(self, *args, **kwargs):
        return self.np.log2(*args, **kwargs)

    def mat(self, *args, **kwargs):
        return self.np.mat(*args, **kwargs)

    def max(self, *args, **kwargs):
        return self.np.max(*args, **kwargs)

    def min(self, *args, **kwargs):
        return self.np.min(*args, **kwargs)

    def ndarray(self, *args, **kwargs):
        return self.np.ndarray(*args, **kwargs)

    def nonzero(self, *args, **kwargs):
        return self.np.nonzero(*args, **kwargs)

    def ones(self, *args, **kwargs):
        return self.np.ones(*args, **kwargs)

    def ones_like(self, *args, **kwargs):
        return self.np.ones_like(*args, **kwargs)

    def power(self, *args, **kwargs):
        return self.np.power(*args, **kwargs)

    def savetxt(self, *args, **kwargs):
        return self.np.savetxt(*args, **kwargs)

    def set_printoptions(self, *args, **kwargs):
        return self.np.set_printoptions(*args, **kwargs)
    
    def sin(self, *args, **kwargs):
        return self.np.sin(*args, **kwargs)

    def squeeze(self, *args, **kwargs):
        return self.np.squeeze(*args, **kwargs)

    def sum(self, *args, **kwargs):
        return self.np.sum(*args, **kwargs)

    def tan(self, *args, **kwargs):
        return self.np.tan(*args, **kwargs)

    def union1d(self, *args, **kwargs):
        return self.np.union1d(*args, **kwargs)

    def vectorize(self, *args, **kwargs):
        return self.np.vectorize(*args, **kwargs)

    def where(self, *args, **kwargs):
        return self.np.where(*args, **kwargs)

    def zeros(self, *args, **kwargs):
        return self.np.zeros(*args, **kwargs)

    def zeros_like(self, *args, **kwargs):
        return self.np.zeros_like(*args, **kwargs)

    def lib_stride_tricks_as_strided(self, x, shape=None, strides=None, writeable=True):
        return self.np.lib.stride_tricks.as_strided(
            x, shape=shape, strides=strides, writeable=writeable
        )
