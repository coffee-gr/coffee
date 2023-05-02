""" This module contains an object which represents 
the spherical harmonic coefficients of a 
function defined on S^2. 

Internally it provides python bindings for revision 104
of the methods listed in alm.h from Huffenberger's & Wandelt's spinsfast code. 
See "Fast and exact spin-s Spherical Harmonic Transformations" in Astrophysical 
Journal Supplement Series (2010) 189:255-260. All references to sections, 
appendices or equations in the documentation below refer to equations given in 
this paper. The original spinsfast code can be found at
http://astrophysics.physics.fsu.edu/~huffenbe/research/spinsfast/index.html.

The object is a subclass of numpy.ndarray which contains the spin spherical 
harmonic coefficients along with some meta data.

The module also contains bindings to the methods given in the file alm.h 
contained in Huffenberger's & Wandelt's spinsfast code.
"""
# future imports


# Standard libraries
import numpy as np
import math
import logging

# Package internal imports
from coffee.swsh import w3j
from coffee.swsh import clebschgordan as cg_mod
from coffee.swsh import salm
from coffee.swsh.spinsfastpy.index_utilities import lm_ind, ind_lm, lmax_Nlm


class sfpy_salm(be.ndarray):
    """Represents the spin spherical harmonic coefficients of a
    function defined on S^2.

    The coefficient of sYlm is accessed as sfpy_salm[s,l,m]. Basic slicing
    is implemented.
    """

    cg_default = cg_mod.CGStone()

    def __new__(cls, array, spins, lmax, cg=None, bandlimit_multiplication=False):
        """Returns an instance of alm.salm.

        It must be the case that:

            len(spins) == len(array[:,0]).

        This is not checked.

        While it is possible to construct instances of this class directly, see
        the file example_mulispin.c of spinsfast, the more usual construction
        will be as the object returned by spinsfastpy.forward.forward.

        Parameters
        ----------
        array : numpy.ndarray
            The alm array with the structure described in
            the doc string for the class.

        spins : numpy.ndarray
            An array of the spin values.

        lmax : float
            The maximum bandwidth.

        cg : Optional
            An object that computes Clebsch Gordan coefficients.

        bandlimit_multiplication : bool, Optional
            If true the bandwidth maximum during multiplication is the sum
            of the limits for the two series. Otherwise the bandwidth limit
            during multiplication is the maximum of the limits of the
            two multiplicands.

        Returns
        -------
        salm : an instance of class alm.salm
        """
        # no error checking is currently done to test if s>=l. Be warned.
        spins = be.asarray(spins)
        if len(array.shape) == 1:
            if spins.shape is () and array.shape[0] == lmax_Nlm(lmax):
                obj = be.asarray(array).view(cls)
            else:
                raise ValueError("Mal-formed array and shape for lmax")
        elif len(array.shape) == 2:
            if spins.shape[0] == array.shape[0] and array.shape[1] == lmax_Nlm(lmax):
                obj = be.asarray(array).view(cls)
            else:
                raise ValueError("Mal-formed array and shape for lmax")
        else:
            raise ValueError("Array shape has too many dimensions")
        obj.spins = spins
        obj.lmax = lmax
        obj.bl_mult = bandlimit_multiplication
        obj.log = logging.getLogger("sfpy_salm")
        # cg should never be None
        if cg is None:
            obj.cg = sfpy_salm.cg_default
        else:
            obj.cg = cg
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.spins = getattr(obj, "spins", None)
        self.lmax = getattr(obj, "lmax", None)
        self.bl_mult = getattr(obj, "bl_mult", None)
        self.cg = getattr(obj, "cg", sfpy_salm.cg_default)
        self.log = getattr(obj, "log", None)

    def __str__(self):
        s = "spins = %s,\n" % repr(self.spins)
        if self.spins.shape is ():
            for l in range(self.lmax + 1):
                s += "(%f, %f): %s\n" % (self.spins, l, repr(self[l].view(be.ndarray)))
        else:
            for spin in self.spins:
                for l in range(self.lmax + 1):
                    s += "(%f, %f): %s\n" % (
                        spin,
                        l,
                        repr(self[spin, l].view(be.ndarray)),
                    )
        return s

    def __repr__(self):
        s = (
            "sfpy_salm object:\nlmax = %d,\nspins = %s,\nbandlimit_multiplication = %s,\ncg = %s,\ncoefficients = \n%s"
            % (
                self.lmax,
                repr(self.spins),
                repr(self.bl_mult),
                repr(self.cg),
                repr(self.view(be.ndarray)),
            )
        )
        return s

    def multiplication_bandlimit(self, bool):
        """Sets the multiplication bandwidth limit.

        Parameters
        ----------
        bool : bool
        """
        self.bl_mult = bool

    def __getslice__(self, start, stop):
        """This solves a subtle bug, where __getitem__ is not called, and all
        the dimensional checking not done, when a slice of only the first
        dimension is taken, e.g. a[1:3]. From the Python docs:
        Deprecated since version 2.0: Support slice objects as parameters
        to the __getitem__() method. (However, built-in types in CPython
        currently still implement __getslice__(). Therefore, you have to
        override it in derived classes when implementing slicing.)

        See: https://stackoverflow.com/questions/14553485/numpy-getitem-delayed-evaluation-and-a-1-not-the-same-as-aslice-1-none
        """
        return self.__getitem__(slice(start, stop))

    def __getitem__(self, key):
        """Returns the salm coefficient given by the key.

        See convert_key for information about key.

        Parameters
        ----------
        key:
            Specifies s, l, m.
        """
        alt_key = self.convert_key(key)
        if len(alt_key) == 1:
            if self.spins.shape is ():
                return self.view(be.ndarray)[alt_key]
            else:
                return sfpy_salm(
                    self.view(be.ndarray)[alt_key],
                    self.spins[alt_key],
                    self.lmax,
                    self.cg,
                    self.bl_mult,
                )
        else:
            return self.view(be.ndarray)[alt_key]

    def __setitem__(self, key, value):
        """Sets the value of the coefficient corresponding to the given values
        of s, l, m in the key.

        See convert_key for information about key.

        Parameters
        ----------
        key:
            Specifies s, l, m.
        """
        alt_key = self.convert_key(key)
        self.view(be.ndarray)[alt_key] = value

    def convert_key(self, key):
        """Converts a 'key' representing (s, l, m) to a tuple of
        indices whose array value is the coefficient for the slm component.

        Valid forms of key are: int or slice, or a two / three tuple of
        ints or slices. If len(key) == 1 it is assumed that key gives l.
        If len(key) == 2 it is assumed that key gives l, m.
        If len(key) == 3 it is assumed that key gives s, l, m.
        Sensible defaults for s and m are assumed if not given.

        Parameters
        ----------
        key:
            Represents the s,l,m coefficient.

        Returns
        -------
        two tuple of ints:
            The first int gives and spinor index and the second gives the
            index corresponding to l and m.
        """
        # convert key into a tuple of ints and slices
        key = be.index_exp[key]

        # perform the transformation from salm to indices
        if len(key) == 1:
            if self.spins.shape is ():
                return (_convert_lm_key(key, self.lmax),)
            else:
                return (_convert_spin_key(key[0], self.spins),)
        elif len(key) >= 2:
            if self.spins.shape is ():
                return (_convert_lm_key(key, self.lmax),)

            else:
                spin_key = _convert_spin_key(key[0], self.spins)
                order_key = _convert_lm_key(key[1:], self.lmax)
                return spin_key, order_key
        else:
            raise IndexError("Too many indices")

    def _addsub(self, other, add):
        """A utility that supports addition or subtraction of salm.salm
        arrays.

        Parameters
        ----------
        other : salm.salm
        add : bool
            True for addition, False for subtraction.

        Returns
        -------
        salm.salm:
            The result of the addition or subtraction.
        """
        # check object types and delegate addition if not sfpy_salm
        if not isinstance(other, sfpy_salm):
            if __debug__:
                self.log.debug("self is = %s" % repr(self))
                self.log.debug("other is = %s" % repr(other))
                self.log.debug("self.spins is = %s" % repr(self.spins))
                self.log.debug("self.lmax is = %s" % repr(self.lmax))
            if add:
                if __debug__:
                    self.log.debug(
                        "self.view(be.adarray) + other = %s"
                        % repr(self.view(be.ndarray) + other)
                    )
                return sfpy_salm(
                    self.view(be.ndarray) + other,
                    self.spins,
                    self.lmax,
                    self.cg,
                    self.bl_mult,
                )
            else:
                if __debug__:
                    self.log.debug(
                        "self.view(be.adarray) - other = %s"
                        % repr(self.view(be.ndarray) - other)
                    )
                return sfpy_salm(
                    self.view(be.ndarray) - other,
                    self.spins,
                    self.lmax,
                    self.cg,
                    self.bl_mult,
                )

        # Perform set up
        # Remember that we have three cases to deal with:
        # 1) other.spins == () and self.spins == ()
        # 2) other.spins == () and self.spins is not ()
        # 3) other.spins is not () and self.spins is not ()
        # I deal with these cases by converting everything to
        # arrays and then formatting the output appropriately.
        lmax = max(self.lmax, other.lmax)
        bl_mult = self.bl_mult or other.bl_mult
        cg = self.cg
        if self.cg is None and other.cg is not None:
            cg = other.cg
        s_spins = self.spins
        s_array = be.asarray(self)
        if self.spins.shape == ():
            s_spins = be.atleast_1d(s_spins)
            s_array = be.atleast_2d(s_array)
        o_spins = other.spins
        o_array = be.asarray(other)
        if other.spins.shape is ():
            o_spins = be.atleast_1d(o_spins)
            o_array = be.atleast_2d(o_array)
        spins = be.union1d(s_spins, o_spins)
        if self.dtype <= other.dtype:
            dtype = other.dtype
        elif self.dtype > other.dtype:
            dtype = self.dtype
        else:
            raise ValueError(
                "Unable to infer the dtype of the result of\
            addition"
            )
        array = be.zeros((spins.shape[0], lmax_Nlm(lmax)), dtype=dtype)
        s_len = lmax_Nlm(self.lmax)
        o_len = lmax_Nlm(other.lmax)

        # Do computation
        for i, spin in enumerate(spins):
            s_spins_select = s_spins == spin
            self_has_spin = any(s_spins_select)
            o_spins_select = o_spins == spin
            other_has_spin = any(o_spins_select)
            if self_has_spin:
                self_index = be.where(s_spins_select)[0][0]
                array[i][:s_len] = s_array[self_index]
                if other_has_spin:
                    other_index = be.where(o_spins_select)[0][0]
                    if add:
                        array[i][:o_len] += o_array[other_index]
                    else:
                        array[i][:o_len] -= o_array[other_index]
            elif other_has_spin:
                other_index = be.where(o_spins_select)[0][0]
                if add:
                    array[i][:o_len] += o_array[other_index]
                else:
                    array[i][:o_len] -= o_array[other_index]
            else:
                raise Exception(
                    "A spin has been encountered that is not in \
                either summand."
                )

        # Now we need to appropriately format the output
        if (
            self.spins.shape is ()
            and other.spins.shape is ()
            and self.spins == other.spins
        ):
            spins = spins[0]
            array = array[0]
        return sfpy_salm(array, spins, lmax, cg, bl_mult)

    def __add__(self, other):
        """Addition.

        Parameters
        ----------
        other : salm.salm

        Returns
        -------
        salm.salm :
        """
        return self._addsub(other, True)

    def __sub__(self, other):
        """Subtraction.

        Parameters
        ----------
        other : salm.salm

        Returns
        -------
        salm.salm :
        """
        return self._addsub(other, False)

    def __rmul__(self, other):
        """rmul.

        Parameters
        ----------
        other : salm.salm

        Returns
        -------
        salm.salm :
        """
        return sfpy_salm(
            other * self.view(be.ndarray), self.spins, self.lmax, self.cg, self.bl_mult
        )

    def __mul__(self, other):
        """Multiplication.

        Parameters
        ----------
        other : salm.salm

        Returns
        -------
        salm.salm :
        """
        # Do we need to treat this multiplication as between salm objects?
        if not isinstance(other, sfpy_salm):
            a = sfpy_salm(
                self.view(be.ndarray) * other,
                self.spins,
                self.lmax,
                self.cg,
                self.bl_mult,
            )
            return a
        # Set up work
        # Transform everything to arrays
        # Then format back when returning the result
        bl_mult = self.bl_mult or other.bl_mult
        cg = self.cg
        if self.cg is None and other.cg is not None:
            cg = other.cg
        if self.cg is None:
            raise Exception(
                "alm.cg must be set to a valid clebschgordan object before multiplication can be done."
            )
        if bl_mult:
            lmax = min(self.lmax, other.lmax)
        else:
            lmax = self.lmax + other.lmax
        s_spins = self.spins
        s_array = be.asarray(self)
        if self.spins.shape == ():
            s_spins = be.atleast_1d(s_spins)
            s_array = be.atleast_2d(s_array)
        self_sorted_spins = sorted(s_spins)
        o_spins = other.spins
        o_array = be.asarray(other)
        if other.spins.shape is ():
            o_spins = be.atleast_1d(o_spins)
            o_array = be.atleast_2d(o_array)
        other_sorted_spins = sorted(o_spins)
        min_spin = self_sorted_spins[0] + other_sorted_spins[0]
        max_spin = self_sorted_spins[-1] + other_sorted_spins[-1]
        spins = be.arange(min_spin, max_spin + 1, 1)
        if self.dtype <= other.dtype:
            dtype = other.dtype
        elif self.dtype > other.dtype:
            dtype = self.dtype
        else:
            raise ValueError(
                "Unable to infer the dtype of the result of\
            multiplication"
            )
        array = be.zeros((spins.shape[0], lmax_Nlm(lmax)), dtype=dtype)
        # Do the calculation already!
        for self_s_index, self_s in enumerate(s_spins):
            for other_s_index, other_s in enumerate(o_spins):
                s = self_s + other_s
                s_index = be.where(spins == s)[0][0]
                for k, self_salm in enumerate(s_array[self_s_index]):
                    self_j, self_m = ind_lm(k)
                    if self_j < abs(self_s):
                        continue
                    for l, other_salm in enumerate(o_array[other_s_index]):
                        other_j, other_m = ind_lm(l)
                        if other_j < abs(other_s):
                            continue
                        m = self_m + other_m
                        jmin = max(
                            abs(self_j - other_j),
                            abs(self_s + other_s),
                            abs(self_m + other_m),
                        )
                        jmax = min(self_j + other_j, lmax)
                        js = be.arange(jmin, jmax + 1, 1)
                        first_f = 0.5 * math.sqrt(
                            (2 * self_j + 1) * (2 * other_j + 1) / math.pi
                        )
                        for j in js:
                            second_f = 1 / math.sqrt(2 * j + 1)
                            cg_fact = self.cg(
                                self_j, self_m, other_j, other_m, j, m
                            ) * self.cg(self_j, -self_s, other_j, -other_s, j, -s)
                            jm_index = lm_ind(j, m)
                            #                            print "(j,m) = (%d, %d)"%(j,m)
                            #                            print "(s_index, jm_index) = (%d, %d)"%(s_index,jm_index)
                            #                            print "value is %f"%(first_f * second_f \
                            #                                * cg_fact * self_salm * other_salm)
                            # Change this calculation to avoid this sum...
                            array[s_index, jm_index] = (
                                array[s_index, jm_index]
                                + first_f * second_f * cg_fact * self_salm * other_salm
                            )
        #                            print "array is now %s"%repr(array)
        if self.spins.shape is () and other.spins.shape is ():
            array = array[0]
            spins = spins[0]
        return sfpy_salm(
            array, spins, lmax, cg=self.cg, bandlimit_multiplication=bl_mult
        )


# salm.Salm.register(sfpy_salm)


def _convert_lm_key(key, lmax):
    """A utility function that maps the variable `key` to
    a slice or index representing the array location of
    the lm coefficient given by the key.

    The method attempts to make a reasonable attempt of interpreting the key.
    Valid values of key are: int, slice, tuple of ints or slices or both.

    Parameters
    ----------
    key : a two tuple
        Represents the l and m variables for a harmonic representation of
        a function.

    Returns
    -------
    tuple of ints or slices:

    Raises
    ------
    IndexError:
        This is raised when an invalid key is given.
    """
    if isinstance(key, int):
        l_ind = lm_ind(key[0], -key[0])
        return slice(l_ind, l_ind + 2 * key[0] + 1)
    elif isinstance(key, slice):
        if key == slice(None, None, None):
            return key
        else:
            raise IndexError("The order index cannot be a slice")
    else:
        if isinstance(key[0], slice):
            if key[0] == slice(None, None, None):
                if len(key) == 1:
                    return key[0]
                elif len(key) == 2 and key[1] == slice(None, None, None):
                    return key
                else:
                    raise IndexError("The order/degree index cannot be a slice")
            raise IndexError("The order index cannot be a slice")
        elif key[0] <= lmax and key[0] >= 0:
            if len(key) == 1:
                l_ind = lm_ind(key[0], -key[0])
                return slice(l_ind, l_ind + 2 * key[0] + 1)
            elif len(key) == 2:
                if isinstance(key[1], slice):
                    if key[1] == slice(None, None, None):
                        l_ind = lm_ind(key[0], -key[0])
                        return slice(l_ind, l_ind + 2 * key[0] + 1)
                    else:
                        raise IndexError("The degree index cannot be a slice")
                elif abs(key[1]) <= key[0]:
                    return lm_ind(key[0], key[1])
                else:
                    raise IndexError("degree out of bounds")
        else:
            raise IndexError("order out of bounds")


def _convert_spin_key(key, spins):
    """A utility method that converts the key variable representing
    s,l,m values to the array indices that contain the slm coefficient.

    Valid values of `key` are int, slice.

    Parameters
    ----------
    key : int or slice

    Returns
    -------
    int or slice:

    Raises
    ------
    IndexError:
        Raised when the key is invalid.
    """
    if isinstance(key, int):
        try:
            spin_key = be.where(spins == key)[0][0]
        except IndexError:
            raise IndexError("Spins %d not found in spins" % key)
    elif isinstance(key, slice):
        if key.start is None:
            start = 0
        else:
            start = be.where(spins == key.start)[0][0]
        if key.stop is None:
            stop = spins.size
        else:
            stop = be.where(spins == key.stop)[0][0]
        spin_key = slice(start, stop, key.step)
    else:
        raise IndexError("spin index may only be an integer or slice")
    return spin_key


class sfpy_sralm(be.ndarray):
    """
    Represents the spin spherical harmonic coefficients of a function defined on
    [a,b] x S^2. It is assume that this array will be stored in an array
    corresponding to spins. Thus, while this object knows it's own spin
    it is not possible for this object to store multiple spins.

    The coefficient of sYlm is accessed as sfpy_salm[:,l,m]. Basic slicing
    is implemented.
    """

    cg_default = None

    def __new__(cls, array, spins, lmax, cg=None, bandlimit_multiplication=False):
        """Returns an instance of salm.sralm.

        The array must have the following structure: len(array.shape)=2 or 3,
        If .... balch
        array.shape[0] = number of spins
        array.shape[1] = number of grid points in the interval
        array.shape[2] = salm.lmax_Nlm(lmax)

        Parameters
        ----------
        array : numpy.ndarray
        spins : numpy.ndarray
            an array of ints giving the spin values.

        Returns
        -------
        sfpy_sralm :
        """
        spins = be.asarray(spins)
        if len(array.shape) == 2:
            if spins.shape is () and array.shape[1] == lmax_Nlm(lmax):
                obj = be.asarray(array).view(cls)
            else:
                raise ValueError("Mal-formed array and shape for lmax")
        elif len(array.shape) == 3:
            if spins.shape[0] == array.shape[0] and array.shape[2] == lmax_Nlm(lmax):
                obj = be.asarray(array).view(cls)
            else:
                raise ValueError("Mal-formed array and shape for lmax")
        else:
            raise ValueError("Array shape has too many or too few dimensions")
        obj.spins = spins
        obj.lmax = lmax
        obj.bl_mult = bandlimit_multiplication
        # cg should never be none
        if cg is None:
            obj.cg = sfpy_salm.cg_default
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.spins = getattr(obj, "spins", None)
        self.lmax = getattr(obj, "lmax", None)
        self.bl_mult = getattr(obj, "bl_mult", None)
        self.cg = getattr(obj, "cg", sfpy_salm.cg_default)

    def __str__(self):
        s = "spins = %s,\n" % repr(self.spins)
        s += "lmax = %d,\n" % self.lmax
        if self.spins.shape is ():
            # for i in range(self.shape[0]):
            for i in [0, 1, 2, self.shape[0] - 3, self.shape[0] - 2, self.shape[0] - 1]:
                for l in range(self.lmax + 1):
                    s += "(%.1f, %.1f) at %d r value: %s\n" % (
                        int(self.spins),
                        l,
                        i,
                        repr(self[i, l]),
                    )
        else:
            for j in self.spins:
                for i in range(self.shape[1]):
                    for l in range(self.lmax + 1):
                        s += "(%.1f, %.1f) at %d r value: %s\n" % (
                            j,
                            l,
                            i,
                            repr(self[j, i, l]),
                        )
        return s

    def __repr__(self):
        s = (
            "<sfpy_sralm lmax = %d, spins = %s, bandlimit_multiplication =\
        %s, cg = %s, coefficients = %s>"
            % (
                self.lmax,
                repr(self.spins),
                repr(self.bl_mult),
                repr(self.cg),
                repr(self.view(be.ndarray)),
            )
        )
        return s

    def __getslice__(self, start, stop):
        """This solves a subtle bug, where __getitem__ is not called, and all
        the dimensional checking not done, when a slice of only the first
        dimension is taken, e.g. a[1:3]. From the Python docs:
        Deprecated since version 2.0: Support slice objects as parameters
        to the __getitem__() method. (However, built-in types in CPython
        currently still implement __getslice__(). Therefore, you have to
        override it in derived classes when implementing slicing.)
        """
        return self.__getitem__(slice(start, stop))

    def __setitem__(self, key, value):
        """Sets the value of the coefficient corresponding to the given values
        of s, l, m in the key.

        See convert_key for information about key.

        Parameters
        ----------
        key:
            Specifies s, l, m.
        """
        alt_key = self.convert_key(key)
        self.view(be.ndarray)[alt_key] = value

    def __getitem__(self, key):
        """Returns the salm coefficient given by the key.

        See convert_key for information about key.

        Parameters
        ----------
        key:
            Specifies s, l, m.
        """
        key = self.convert_key(key)
        if self.spins.shape is ():
            return self._getitem_no_spin(key)
        else:
            return self._getitem_spin(key)

    def _getitem_no_spin(self, key):
        rv = be.asarray(self)[key]
        if len(key) == 1:
            if len(rv.shape) == 2:
                if rv.shape[0] == 1:
                    return sfpy_salm(
                        rv[0], self.spins, self.lmax, self.cg, self.bl_mult
                    )
                else:
                    return sfpy_sralm(rv, self.spins, self.lmax, self.cg, self.bl_mult)
            else:
                return sfpy_salm(rv, self.spins, self.lmax, self.cg, self.bl_mult)
        elif len(key) == 2:
            if (
                not isinstance(rv, be.ndarray)
                or len(rv.shape) < 2
                or rv.shape[1] != lmax_Nlm(self.lmax)
            ):
                return rv
            else:
                if rv.shape[0] == 1:
                    return sfpy_salm(
                        rv[0], self.spins, self.lmax, self.cg, self.bl_mult
                    )
                else:
                    return sfpy_sralm(rv, self.spins, self.lmax, self.cg, self.bl_mult)
        raise IndexError("Unable to process selection")

    def _getitem_spin(self, key):
        rv = be.asarray(self)[key]
        rs = self.spins[key[0]]
        if len(key) < 2:
            return sfpy_sralm(rv, rs, self.lmax, self.cg, self.bl_mult)
        elif len(key) == 3:
            return rv
        elif len(key) == 2:
            # This is the awkward case. Did we get rid of s or r dependence?
            # note that r dependence can be removed in two ways, either the
            # rv array has reduced shape length or length of the dimension
            # representing r is 1.
            if rs.shape is ():
                # spins has disapeared we know that len(rv.shape) <=2
                if len(rv.shape) < 2:
                    # since the key is selecting for only s and r we know that
                    # rv must contain all l values for the selected s and r
                    # this is an salm object
                    return sfpy_salm(rv, self.spins, self.lmax, self.cg, self.bl_mult)
                else:
                    # we need to test to see if r dependence has been removed
                    if rv.shape[0] == 1:
                        return sfpy_salm(
                            rv[0], self.spins, self.lmax, self.cg, self.bl_mult
                        )
                    else:
                        return sfpy_sralm(rv, rs, self.lmax, self.cg, self.bl_mult)
            else:
                # spin dependence still exists
                if len(rv.shape) <= 2:
                    # the r dependence has been removed
                    return sfpy_salm(rv, self.spins, self.lmax, self.cg, self.bl_mult)
                elif rv.shape[1] == 1:
                    return sfpy_salm(
                        rv[:, 0, :], self.spins, self.lmax, self.cg, self.bl_mult
                    )
                else:
                    return sfpy_sralm(rv, rs, self.lmax, self.cg, self.bl_mult)
        # if nothing was selected return an error
        raise IndexError("Unable to process selection")

    def convert_key(self, key):
        """Converts a 'key' representing (s, l, m) to a tuple of
        indices whose array value is the coefficient for the slm component.

        Valid forms of key are: int or slice, or a two / three tuple of
        ints or slices. If len(key) == 1 it is assumed that key gives l.
        If len(key) == 2 it is assumed that key gives l, m.
        If len(key) == 3 it is assumed that key gives s, l, m.
        Sensible defaults for s and m are assumed if not given.

        Parameters
        ----------
        key:
            Represents the s,l,m coefficient.

        Returns
        -------
        three tuple of ints:
            The first int gives and spinor index and the second gives the
            index corresponding to l and m.
        """
        key = be.index_exp[key]
        if len(key) == 1:
            if self.spins.shape is ():
                return key
            return (_convert_spin_key(key[0], self.spins),)
        elif len(key) == 2:
            if self.spins.shape is ():
                return key[0], _convert_lm_key(key[1:], self.lmax)
            return _convert_spin_key(key[0], self.spins), key[1]
        elif len(key) >= 3:
            if self.spins.shape is ():
                return key[0], _convert_lm_key(key[1:], self.lmax)
            else:
                return (
                    _convert_spin_key(key[0], self.spins),
                    key[1],
                    _convert_lm_key(key[2:], self.lmax),
                )
        else:
            return IndexError("too many indices")

    def __add__(self, other):
        """Addition.

        Parameters
        ----------
        other:

        Returns
        -------
        salm.salm:
        """
        if isinstance(other, sfpy_sralm):
            if self.spins.shape == ():
                rv_temp = self[0] + other[0]
                rv = be.empty((self.shape[0],) + (rv_temp.shape[0],), dtype=self.dtype)
                rv[0] = rv_temp
                for i in range(1, rv.shape[0]):
                    rv[i] = self[i] + other[i]
            else:
                rv_temp = self[:, 0] + other[:, 0]
                rv = be.empty(
                    (rv_temp.shape[0],) + (self.shape[1],) + (rv_temp.shape[1],),
                    dtype=self.dtype,
                )
                rv[:, 0] = rv_temp
                for i in range(1, rv.shape[1]):
                    rv[:, i] = self[:, i] + other[:, i]
            return sfpy_sralm(
                rv,
                rv_temp.spins,
                rv_temp.lmax,
                cg=rv_temp.cg,
                bandlimit_multiplication=rv_temp.bl_mult,
            )
        else:
            return sfpy_sralm(
                self.view(be.ndarray) + other,
                self.spins,
                self.lmax,
                cg=self.cg,
                bandlimit_multiplication=self.bl_mult,
            )

    def __radd__(self, other):
        """radd.

        Parameters
        ----------
        other:

        Returns
        -------
        salm.salm:
        """
        return sfpy_sralm(
            other + be.asarray(self),
            self.spins,
            self.lmax,
            cg=self.cg,
            bandlimit_multiplication=self.bl_mult,
        )

    def __sub__(self, other):
        """Subtraction.

        Parameters
        ----------
        other:

        Returns
        -------
        salm.salm:
        """
        if isinstance(other, sfpy_sralm):
            if self.spins.shape == ():
                rv_temp = self[0] - other[0]
                rv = be.empty((self.shape[0],) + (rv_temp.shape[0],), dtype=self.dtype)
                rv[0] = rv_temp
                for i in range(1, rv.shape[0]):
                    rv[i] = self[i] - other[i]
            else:
                rv_temp = self[:, 0] - other[:, 0]
                rv = be.empty(
                    (rv_temp.shape[0],) + (self.shape[1],) + (rv_temp.shape[1],),
                    dtype=self.dtype,
                )
                rv[:, 0] = rv_temp
                for i in range(1, rv.shape[0]):
                    rv[:, i] = self[:, i] - other[:, i]
            return sfpy_sralm(
                rv,
                rv_temp.spins,
                rv_temp.lmax,
                cg=rv_temp.cg,
                bandlimit_multiplication=rv_temp.bl_mult,
            )
        else:
            return sfpy_sralm(
                self.view(be.ndarray) - other,
                self.spins,
                self.lmax,
                cg=self.cg,
                bandlimit_multiplication=self.bl_mult,
            )

    def __rsub__(self, other):
        """rsub.

        Parameters
        ----------
        other:

        Returns
        -------
        salm.salm:
        """
        return sfpy_sralm(
            other - be.asarray(self),
            self.spins,
            self.lmax,
            cg=self.cg,
            bandlimit_multiplication=self.bl_mult,
        )

    def __mul__(self, other):
        """Multiplication.

        Parameters
        ----------
        other:

        Returns
        -------
        salm.salm:
        """
        if isinstance(other, sfpy_sralm):
            rv_temp = self[:, 0] * other[:, 0]
            rv = be.empty(
                (rv_temp.shape[0],) + (self.shape[1],) + (rv_temp.shape[1],),
                dtype=self.dtype,
            )
            rv[:, 0] = rv_temp
            for i in range(1, rv.shape[0]):
                rv[:, i] = self[:, i] * other[:, i]
            return sfpy_sralm(
                be.ndarray(rv),
                rv_temp.spins,
                rv_temp.lmax,
                cg=rv_temp.cg,
                bandlimit_multiplication=rv_temp.bl_mult,
            )
        else:
            return sfpy_sralm(
                self.view(be.ndarray) * other,
                self.spins,
                self.lmax,
                cg=self.cg,
                bandlimit_multiplication=self.bl_mult,
            )

    def __rmul__(self, other):
        """rmul.

        Parameters
        ----------
        other:

        Returns
        -------
        salm.salm:
        """
        return sfpy_sralm(
            other * self.view(be.ndarray),
            self.spins,
            self.lmax,
            cg=self.cg,
            bandlimit_multiplication=self.bl_mult,
        )


# salm.Salm.register(sfpy_sralm)

# if __name__ == "__main__":
# lmax = 2
# a = be.ones((2,lmax_Nlm(lmax)))
# spins = 2
# p = sfpy_sralm(a, spins, lmax)
# p[-1]
# print "\n\n"
# p[-1:]
# print "\n\n"
# p[slice(-1,None,None)]
# from coffee.swsh import clebschgordan as cg
# cg_object = cg.CGStone()

# lmax = 3
# spins = -1
# Nlmax = lmax_Nlm(lmax)
# a = be.arange(Nlmax)
# a[:] =1
##a = be.array([a])
# a_salm = sfpy_salm(a, spins, lmax, cg=cg_object)
# print "a info"
# print a_salm
# print a_salm.spins
# print a_salm.shape

# lmax = 3
# spins = -1
# Nlmax = lmax_Nlm(lmax)
# b = be.arange(Nlmax)
# b = be.array(b)
# b_salm = sfpy_salm(b, spins, lmax, cg=cg_object)
# print "\n\nb info"
# print b_salm
# print b_salm.spins
# print b_salm.shape

# c_salm = a_salm * b_salm
# print "\n\nc info"
# print c_salm
# print c_salm.spins
# print c_salm.shape

# print "\n\n"
# print c_salm[5,-2]
