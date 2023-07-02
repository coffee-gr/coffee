#!/usr/bin/env python
# encoding: utf-8
"""
A module that contains dissipation operators for SBP operators.
"""

import math
import abc

################################################################################
# Base classes for SBP artificial dissipation operators
################################################################################
class Diss(object):
    """The parent class for SBP dissipation operators."""

    name = "Dissipation"

    def __init__(self, *args, **kwds):
        pass

    @abc.abstractmethod
    def __call__(self, u, dx, boundary_ID=None):
        pass

    def _consistancy_check(self, u, shape):
        r, c = shape
        if u.shape[0] + 1 <= 2 * c:
            raise ValueError("Domain too small for application of operator")

    def __str__(self):
        return "SBP Dissipation operator" % self.name


class DissDiag(Diss):
    """A class that represents dissipation operators that are diagonal."""

    name = "Dissipation Diagonal Norm"

    def __init__(self):
        self.bdyRegion = self.Ql.shape

    def __call__(self, u, dx, boundary_ID=None):
        super(DissDiag, self)._consistancy_check(u, self.bdyRegion)
        r, c = self.bdyRegion
        diss_u = be.zeros_like(u)
        diss_u = be.convolve(u, self.A, mode="same")
        if boundary_ID is None:
            diss_u[0:r] = be.dot(self.Ql, u[0:c])
            diss_u[-r:] = be.dot(self.Qr, u[-c:])
        elif boundary_ID == grid.LEFT:
            diss_u[0:r] = be.dot(self.Ql, u[0:c])
        elif boundary_ID == grid.RIGHT:
            diss_u[-r:] = be.dot(self.Qr, u[-c:])
        factor = math.pow(0.5, 2 * self.p)
        return factor * diss_u


class DissRestFull(Diss):
    """A class that represents dissipation operators for restricted full norms."""

    name = "Dissipation Restricted Full Norm"

    def __init__(self, p):
        """The Initialiser for this class.

        Parameters
        ----------
        p : int
            The power of the operator.
        """
        self.bdyRegion = self.Ql.shape
        self.p = p

    def __call__(self, u, dx, boundary_ID=None):
        """Applies the operator.

        We follow the notation of Diener Dorband Schnetter and Tiglio. We
        make the assumption that the differencing operator D_p is always of the
        form of some central finite difference stencil in the interior and some
        boundary matrix on the boundary.

        An important assumption is that D decomposes into the sum of two
        matrices A, Q. Where A represents the application of a centred
        finite difference operator on the interior and Q the application of the
        centred difference on the boundary.

        If the stencil for A is [-1,3, -3, 1] then it is assumed that
        D has the form:

            [[-1, 3, -3, 1, 0,...], [-1, 3, -3, 1, 0,...], [-1, 3, -3, 1, 0,...],
             [0, -1, 3, -3, 1, 0,...], [0, 0, -1, 3, -3, 1, 0, ...], ... ]

        and that the number of times the stencil is repeated on the boundary
        is the stencil length (4 in this case) minus 1 (to give 3 repeats).
        Hence Q is assumed to have the form:

            [[-1, 3, -3, 1, 0,...], [-1, 3, -3, 1, 0,...], [0,0,0,0,0,...],
             [0,0,0,0,0,0,...], ... ]

        and A is assumed to be D-Q. This ensures that if
        r,c = self.bdyRegion then r = self.A.shape[0]//2. This is an important
        assumption that is applied in the formula below to calculate
        Transpose[A].B.A.u where B is the non-constant diagonal matrix required
        to ensure numerical stability (see the paper).

        The explicit form of A and Q here ensures that Transpose[A].Q = 0
        and Transpose[Q].A = 0. This implies that
        Transpose[D].B.D.u = Transpose[Q].B.Q.u + Transpose[A].B.A.u
        This equation reduces the need for further calculations involving
        Transpose[A].Q and Transpose[Q].A

        In principle it is possible to remove these restrictions, i.e. allow
        r neq self.A.shape[0]//2 and to allow Transpose[A].Q neq 0
        and Transpose[Q].A neq 0, but the code below WILL need to be
        rewritten.

        Parameters
        ----------
        u: tslice.TimeSlice
            The data
        ds : float
            The spatial step size
        boundary_ID: int
            See the SBP operator module.

        Returns
        -------
        numpy.ndarray:
        """
        # This performs a consistance check
        super(DissRestFull, self)._consistancy_check(u, self.Ql.shape)

        # Get useful variables
        size = u.shape[0]
        r, c = self.Ql.shape

        # First do the internal convolution.
        # The array A is the central difference stencil that makes up the
        # "interior" portion of D_p. By assumption A is a matrix with the first
        # r rows full of zeros and the first non-zero row starting with the
        # stencil. That is if the stencil is [1, -2, 1] then the top left corner
        # of A is [[0,0,0,...], [0,0,0,...], [1, -2, 1, 0,...], ...
        # This implies that the convolvution operator used to implement matrix
        # multiplication must take account of the zero rows.
        # It turns out that this can be done by using the 'valid' mode in
        # the be.convolve method. This will reduce the size of the output by r
        # elements on the end and beginning of the output array. The values that
        # should be in these positions are all zero due to the r zero'd rows in
        # the matrix A.
        diss_u_int = be.zeros_like(u)
        diss_u_int[r:-r] = be.convolve(u, self.A, mode="valid")
        for i in range(size):
            diss_u_int[i] = self.B(i, dx, size) * diss_u_int[i]

        # To multiply my Transpose[A], requires some care. The zero'd rows now
        # become columns and the stencil is reversed. I make here the assumption
        # that the stencil is always symmetric, up to sign! For stencils with
        # odd length (p is even) the sign does not change, but for stencils with
        # even length (p is odd) the sign does change. For the stencil
        # [-1,3,-3,1] the matrix A is
        # [[0,0,0,0,0,...], [0,0,0,0,0,...], [-1, 3, -3, 1, 0,...],
        # [0, -1, 3, -3, 1, 0,...], ...]
        # Hence Transpose[A] is
        # [[0, 0, -1, 0,...], [0, 0, 3, -1, 0,...], [0, 0, -3, 3, -1, 0,...],
        # [0, 0, 1, -3, 3, -1, 0,...],... ]
        # This is incontrast to the A matrix generated by the stencil [1,-2,1]
        # where Transpose[A] is
        # [[0, 1, 0,...], [0, -2, 1, 0,...], [0, 1, -2, 1, 0,...], ...]
        # This is why the term math.pow(-1, p) is incorporated.
        #
        # On top of this be.convolve performs the convolution iterating from the
        # end of the stencil. This makes a difference when boundary affects of
        # convolution are desirable (as in this case). This is why there is
        # a ::-1 in the iteration for u in convolve. This did not need to be
        # accounted for in the calculation above because of the use of the
        #'valid' mode.
        #
        # Lastly the differences in vector length diss_u[r:-r] and diss_u is
        # inculded to account for the vector lengthening of the 'full' mode
        # of the convolve method which adds r number of terms to the vector
        # at both ends.
        diss_u_int = be.convolve(
            diss_u_int[r:-r], math.pow(-1, self.p) * self.A[::-1], mode="full"
        )

        # Second do the boundary convolution
        # The check at the beginning of the method ensures that the
        # two array's diss_u_b[0:r] and diss_u_b[-r:] do not share any points
        # of diss_u_b
        diss_u_b = be.zeros_like(u)
        if boundary_ID is None:
            diss_u_b[0:r] = be.dot(self.Ql, u[0:c])
            diss_u_b[-r:] = be.dot(self.Qr, u[-c:])
            for i in range(size):
                diss_u_b[i] = self.B(i, dx, size) * diss_u_b[i]
            diss_u_b[0:c] = be.dot(diss_u_b[0:r], self.Ql)
            diss_u_b[-c:] = be.dot(diss_u_b[-r:], self.Qr)
        elif boundary_ID == grid.LEFT:
            diss_u[0:r] = be.dot(self.Ql, u[0:c])
            for i in range(size):
                diss_u_b[i] = self.B(i, dx, size) * diss_u_b[i]
            diss_u_b[0:c] = be.dot(diss_u_b[0:r], self.Ql)
        elif boundary_ID == grid.RIGHT:
            diss_u[-r:] = be.dot(self.Qr, u[-c:])
            for i in range(size):
                diss_u_b[i] = self.B(i, dx, size) * diss_u_b[i]
            diss_u_b[-c:] = be.dot(diss_u_b[-r:], self.Qr)

        # Add the two parts together and multiply by the appropriate
        # numerical factor.
        diss_u = diss_u_int + diss_u_b
        factor = math.pow(0.5, 2 * self.p) * math.pow(dx, 2 * self.p - 2)
        diss_u = -factor * diss_u

        # Multiply by the inverse of the norm
        super(DissRestFull, self)._consistancy_check(u, self.norm_inv.shape)
        r, c = self.norm_inv.shape
        diss_u[0:r] = be.dot(self.norm_inv, diss_u[0:c])
        diss_u[-r:] = be.dot(self.norm_inv[::-1, ::-1], diss_u[-c:])

        # return the result
        return diss_u


################################################################################
# Dissipation for diagonal norm SBP operators.
################################################################################
class Diss21_DDST(DissDiag):
    """A dissipation operator.

    Taken from Diener, Dorband, Schnetter and Tiglio. Due to consistence
    conditions on SBP dissipation operators and the derivative operators this
    dissipation operator can only be used with the D21 operator from the same
    paper. Note that D21 is unique so the operator sbp.D21_CNG is the
    correct operator.
    """

    def __init__(self, *args, **kwds):
        self.A = be.array([1.0, -2.0, 1.0])
        self.name = "Diss21_DDST"
        self.p = 2

        self.Ql = be.zeros((2, 3))

        self.Ql[0, 0] = -2
        self.Ql[0, 1] = 2.0
        self.Ql[0, 2] = 0.0
        self.Ql[1, 0] = 1.0
        self.Ql[1, 1] = -2.0
        self.Ql[1, 2] = 1.0

        self.Qr = -self.Ql[::-1, ::-1]
        self.bdyRegion = (2, 3)


class Diss42_DDST(DissDiag):
    """A dissipation operator.

    Taken from Diener, Dorband, Schnetter and Tiglio. Due to consistency
    conditions on SBP dissipation operators and the derivative operators this
    disspation operator can only be used with the sbp.D42 operator which
    is from the same
    paper.
    """

    def __init__(self, *args, **kwds):
        self.A = be.array([-1.0, 4.0, -6.0, 4.0, -1.0])
        self.name = "Diss42_DDST"
        self.p = 4

        Q = be.zeros((4, 6))

        Q[0, 0] = -2.8235294117647058823529411764705882352941176470588
        Q[0, 1] = 5.6470588235294117647058823529411764705882352941176
        Q[0, 2] = -2.8235294117647058823529411764705882352941176470588
        Q[0, 3] = 0.0
        Q[0, 4] = 0.0
        Q[0, 5] = 0.0

        Q[1, 0] = 1.6271186440677966101694915254237288135593220338983
        Q[1, 1] = -4.0677966101694915254237288135593220338983050847458
        Q[1, 2] = 3.2542372881355932203389830508474576271186440677966
        Q[1, 3] = -0.81355932203389830508474576271186440677966101694915
        Q[1, 4] = 0.0
        Q[1, 5] = 0.0

        Q[2, 0] = -1.1162790697674418604651162790697674418604651162791
        Q[2, 1] = 4.4651162790697674418604651162790697674418604651163
        Q[2, 2] = -6.6976744186046511627906976744186046511627906976744
        Q[2, 3] = 4.4651162790697674418604651162790697674418604651163
        Q[2, 4] = -1.1162790697674418604651162790697674418604651162791
        Q[2, 5] = 0.0

        Q[3, 0] = 0.0
        Q[3, 1] = -0.97959183673469387755102040816326530612244897959184
        Q[3, 2] = 3.9183673469387755102040816326530612244897959183673
        Q[3, 3] = -5.8775510204081632653061224489795918367346938775510
        Q[3, 4] = 3.9183673469387755102040816326530612244897959183673
        Q[3, 5] = -0.97959183673469387755102040816326530612244897959184

        self.Ql = Q
        self.Qr = -self.Ql[::-1, ::-1]
        self.bdyRegion = (4, 6)
        super(Diss42_DDST, self).__init__()


################################################################################
# Dissipation for restricted full norm SBP operators.
################################################################################
class Diss43_DDST(DissRestFull):
    """A dissipation operator.

    Compatibility considerations require that this disspation operator only
    used with sbp.D43_Tiglioetal.
    """

    def __init__(self, bdy_percent):
        """Iniitialisation for Diss43_DDST.

        Parameters
        ----------
        bdy_percent : float
            A tuning parameter for the cut_off parameter used in the matrix
            B for numerical stability.
        """
        self.bdy_percent = bdy_percent
        self.p = 2

        self.A = be.array([1, -2, 1])
        self.Ql = be.array([[1, -2, 1]])
        self.Qr = self.Ql

        self.norm_inv = be.zeros((5, 5))
        self.norm_inv[0, 0] = 4.186595370392226897362216859769846226369
        self.norm_inv[0, 1] = 0
        self.norm_inv[0, 2] = 0
        self.norm_inv[0, 3] = 0
        self.norm_inv[0, 4] = 0
        self.norm_inv[1, 0] = 0
        self.norm_inv[1, 1] = 0.6725191921225620731888714836983116420871
        self.norm_inv[1, 2] = 0.3613418181134949259370502966736306984367
        self.norm_inv[1, 3] = -0.2021316117293899791481674539631879662707
        self.norm_inv[1, 4] = 0.03455320708729270824077678274955265350304
        self.norm_inv[2, 0] = 0
        self.norm_inv[2, 1] = 0.3613418181134949259370502966736306984367
        self.norm_inv[2, 2] = 0.7206133711630147057720442098623847362950
        self.norm_inv[2, 3] = 0.1376472340546569368321616389764958792591
        self.norm_inv[2, 4] = -0.04136405531324488624637892257286207044784
        self.norm_inv[3, 0] = 0
        self.norm_inv[3, 1] = -0.2021316117293899791481674539631879662707
        self.norm_inv[3, 2] = 0.1376472340546569368321616389764958792591
        self.norm_inv[3, 3] = 0.9578653607931026822074133441449909689509
        self.norm_inv[3, 4] = 0.02069353627247161734563597102894256809696
        self.norm_inv[4, 0] = 0
        self.norm_inv[4, 1] = 0.03455320708729270824077678274955265350304
        self.norm_inv[4, 2] = -0.04136405531324488624637892257286207044784
        self.norm_inv[4, 3] = 0.02069353627247161734563597102894256809696
        self.norm_inv[4, 4] = 0.9908272703370861473007798925906968380654

        super(Diss43_DDST, self).__init__(self.p)

    def B(self, i, dx, size):
        """Returns an element of the matrix `B` which ensures numerical stability.

        Parameters
        ----------
        i: int
            The row and column, of the diagonal matrix, to return
        dx: float
            The step size
        size: int
            Length of the data
        """
        cut_off = math.floor(size * self.bdy_percent)
        if i < cut_off:
            return dx + (i / cut_off) * (1 - dx)
        elif i > size - 1 - cut_off:
            return 1 + (cut_off + 1 - size + i) * (dx - 1) / cut_off
        else:
            return 1