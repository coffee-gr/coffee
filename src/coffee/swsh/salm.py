""" 
This module attempts to hide implementation details regarding spectral decompositions
of spin weighted spherical harmonic functions.
"""
# Standard package imports
from abc import ABCMeta, abstractmethod


class Salm(object, metaclass=ABCMeta):
    """The Abstract Base Class for representations of the spin weighted
    spherical decomposition of a function with bandwidth limit."""

    @abstractmethod
    def __repr__(self):
        """A string representation of the concrete class."""

    @abstractmethod
    def multiplication_bandlimit(self, bool):
        """Set the bandlimit behaviour for multiplication of salm objects.

        Parameters
        ----------
        bool - True = the minimum bandlimit of both salm objects is used.
               False = the bandlimit is the sum of the bandlimits of both
                       salm objects"""
        self.bl_mult = bool

    @abstractmethod
    def __getitem__(self, key):
        """Retrieve the coefficients of the decomposition of the function.

        All salm objects should be accessible as salm[s,l,m]. We leave it to
        the developers discretion if slicing, and what type of slicing, is
        implemented.

        Parameters
        ----------
        key :
            An object that specifies which component or components should
            be considered. See implementations for key details.

        Returns
        -------
        float :
            The coefficient of the sYlm component of the function.
        """
        pass

    def check_lm(l, m):
        """Determines if the l and m values are valid.

        Parameters
        ----------
        l : float
            A float representation of an integer or half-integer.
        m : float
            A float representation of an integer or half-integer.

        Returns
        -------
        bool :
            If l and m are valid then return True, else return False.
        """
        if l > self.lmax:
            raise ValueError("l must be less than or equal to lmax")
        if math.abs(m) > l:
            raise ValueError("|m| must be less than or equal to l")

    @abstractmethod
    def __setitem__(self, key, value):
        """Gets the coefficients of the decomposition of the function.

        All salm objects should be accessible as salm[s,l,m]. We leave it to
        the developers discretion if slicing, and what type of slicing, is
        implemented.

        Parameters
        ----------
        key :
            An object that specifies which component or components should
            be considered. See implementations for key details.

        Returns
        -------
        float :
            The coefficient of the sYlm component of the function.
        """
        pass

    @abstractmethod
    def __add__(self, other):
        """Adds two salm objects.

        Parameters
        ----------
        other:
            The thing to add to this object. See implementations for valid
            values.
        """
        pass

    @abstractmethod
    def __sub__(self, other):
        """Subtracts two salm objects.

        Parameters
        ----------
        other:
            The thing to add to this object. See implementations for valid
            values.
        """
        pass

    @abstractmethod
    def __mul__(self, other):
        """Multiplies two salm objects.

        Parameters
        ----------
        other:
            The thing to add to this object. See implementations for valid
            values.
        """
        pass
