#!/usr/bin/env python
# encoding: utf-8
"""
A module to manage summation by parts finite difference operators.
"""

import math
import logging
from coffee.settings import be

################################################################################
# Base class for SBP operators
################################################################################
BOUNDARY_TYPE_SAT = 0
"""Specifies that the simulataneous approximation term method will be used
for boundaries."""
BOUNDARY_TYPE_GHOST_POINTS = 1
"""Specifies that ghost points will be used for boundaries."""


def validate_boundary_type(boundary_type):
    """The method validates the choice of boundary type.

    Parameters
    ==========
    boundary_type : int

    Returns
    =======
    int:
        The given boundary_type, if valid

    Raises
    ======
    ValueError:
        If an incorrect boundary type is specified.
    """
    if (
        boundary_type == BOUNDARY_TYPE_SAT
        or boundary_type == BOUNDARY_TYPE_GHOST_POINTS
    ):
        return boundary_type
    raise ValueError("Invalid BOUNDARY_TYPE for SBP encountered.")


################################################################################
# Base class for SBP operators
################################################################################


class SBP(object):
    """The class that models summation by parts finite difference operators.

    The values of g_{00} and g_{NN} given here
    are based on the values given in theorem 2.1 of Strand's paper,
    "Summation by parts for finite difference approximations for d/dx."

    It is highly unlikely, therefore, that they will ever need to change.
    """

    name = "Dx"
    g00 = -1.0
    gNN = 1.0
    pbound = None

    def __init__(self, boundary_type=BOUNDARY_TYPE_SAT):
        """SBP operators can use different methods to handle internal boundaries
        caused by running the code via mpi.

        Parameters
        ==========
        boundary_type: int
        """
        self.bdyRegion = self.Ql.shape
        self.w = self.A.shape[0] // 2
        self.log = logging.getLogger("SBP")
        self.boundary_type = validate_boundary_type(boundary_type)

    def __call__(self, u, dx, boundary_ID=None):
        """Calculates the derivative.

        Parameters
        ==========
        u : tslice.TimeSlice
            The function data.
        dx : float
            The step size.
        boundary_ID: int, Optional
            One of None, grid.LEFT or grid.RIGHT. Specifies if only the
            right or left (or both) derivatives should be applied at the
            boundary.
        """
        r, c = self.bdyRegion
        if u.shape[0] + 1 <= 2 * c:
            self.log.error("Domain too small for application of operator")
            raise ValueError("Domain too small for application of operator")
        du = be.convolve(u, self.A, mode="same")
        if boundary_ID is None:
            du[0:r] = be.dot(self.Ql, u[0:c])
            du[-r:] = be.dot(self.Qr, u[-c:])
        elif boundary_ID == grid.LEFT:
            du[0:r] = be.dot(self.Ql, u[0:c])
        elif boundary_ID == grid.RIGHT:
            du[-r:] = be.dot(self.Qr, u[-c:])
        return du / (dx**self.order)

    def penalty_boundary(self, dx, vector_selection):
        """Returns the penalty for use with penalty boundaries.

        The vector self.pbound is, in the notation of CGA,
        given by P^(-1)H^(-1)e_i, i=0,1 where
        e_0 = (1,0,...,0)^T and
        e_1 = (0,...,0,1)^T.

        If u = (u_0, ..., u_n) is the discretised function then
        the vector e_0 should be selected (vector_selection = 0) when the
        boundary value is to be applied to u_0. Similarly, the vector
        e_1 (vector_selection = 1) should be selected when the boundary value
        is to be applied to u_n.

        To support the user vector_selection = "right" is aliased to
        vector_selection = 1 and vector_selection = "left" is aliased to
        vector_selection = 0

        When creating an SBP operator ensure that self.pbound is set up
        correctly in the init() method of the operators.

        We return only the portion of the array, to be applied to
        d/dt of u, that is non-zero. In the case of e_0 the vector should be
        applied as

        (d/dt u)[:n] = - tau * characteristic * result of this method.

        In the case of e_1 the vector should be applied as

        (d/dt u)[-n:] = - tau * characteristic * result of this method.

        The necessary change of orientation of the result of this method is
        taken care of here. For further details about the implementation of
        the penalty boundary method we recommend, "TIME-STABLE BOUNDARY
        CONDITIONS FOR FINITE-DIFFERENCE SCHEMES SOLVING HYPERBOLIC SYSTEMS:
        METHODOLOGY AND APPLICATION TO HIGH-ORDER COMPACT SCHEMES" by
        Carpenter, Gottlieb and Abarbanel. (Sorry for the caps).

        Note that this method will only work for first order SBP operators,
        currently.

        Parameters
        ==========
        dx : float
            The spatial step size
        vector_selection : int
        """
        if vector_selection == 1 or vector_selection == "right":
            return self.gNN * self.pbound[::-1] / dx
        elif vector_selection == 0 or vector_selection == "left":
            return self.g00 * self.pbound / dx
        else:
            return Exception(
                "vector_selection must be either 0 or 1. Please \
            see the penalty_boundary doc string for further details."
            )

    def ghost_points(self):
        """Ghost points required for the penalty boundary method.

        Returns
        =======
        two tuple of ints:
            The number of ghost points for each boundary.
        """
        if self.boundary_type == BOUNDARY_TYPE_GHOST_POINTS:
            r, _ = self.bdyRegion
            return r, r
        if self.boundary_type == BOUNDARY_TYPE_SAT:
            return 0, 1
        raise ValueError("Unknown boundary type encountered.")

    def internal_points(self):
        """Internal points required for the penalty boundary method.

        Returns
        =======
        two tuple of ints:
            The number of ghost points for each boundary.
        """
        if self.boundary_type == BOUNDARY_TYPE_GHOST_POINTS:
            r, _ = self.bdyRegion
            return r, r
        if self.boundary_type == BOUNDARY_TYPE_SAT:
            return 1, 1
        raise ValueError("Unknown boundary type encountered.")

    def __str__(self):
        return "Differential operator " % self.name

    # def save(self):
    #     """Outputs a textual representation of the operator in the users
    #     home directory.
    #     """
    #     filename = os.path.expanduser("~/" + self.name)
    #     print(filename)
    #     be.savetxt(filename + "_left.txt", self.Ql)
    #     be.savetxt(filename + "_right.txt", self.Qr)
    #     be.savetxt(filename + "_mid.txt", self.A)


################################################################################
# Second order
################################################################################


class D21_CNG(SBP):
    """An SBP operator that is second order accurate internal and first
    order accurate on the boundary.

    Taken from "A stable and conservative interface treatment of arbitrary
    spatial accuracy", Carpenter, Nordstrom, and Gottlieb.

    The inner product is the identity.
    """

    def __init__(self, *args, **kwargs):
        self.A = -be.array([-1.0 / 2, 0.0, 1.0 / 2])
        self.name = "D21_CNG"
        self.order = 1

        Q = be.zeros((2, 3))

        Q[0, 0] = -1
        Q[0, 1] = 1.0
        Q[0, 2] = 0.0

        Q[1, 0] = -1.0
        Q[1, 1] = 0.0
        Q[1, 2] = 1

        Q = 0.5 * Q

        P = be.zeros((2, 2))
        P[0, 0] = 0.5
        P[0, 1] = 0
        P[1, 0] = 0
        P[1, 1] = 1

        Pinv = be.linalg.inv(P)
        self.pbound = Pinv[:, 0]
        self.Ql = be.dot(Pinv, Q)
        self.Qr = -self.Ql[::-1, ::-1]

        super(D21_CNG, self).__init__(*args, **kwargs)


################################################################################
# Fourth order accurate differential operators
################################################################################


class D42(SBP):
    """This is an SBP operator which is fourth order accurate on the interior
    and second order accurate at the boundary.

    This operator is the D42 operator given in
    the paper, "Optimized high-order derivative and dissipation operators
    satisfying summation by parts, and applications in three-dimensional
    multi-block evolutions" by Diener, Dorband, Schnetter and Tiglio.

    More detail on this operator can be found on page 59, of Strand's paper,
    "Summation by parts for finite difference approximations for d/dx"
    under the heading "Second-order accuracy at the boundary". This
    paper also gives the norm used.

    The norm in this case is given as be.diag([17./48,59./48,43./48,49./48]).
    Note the additional factors included in the code for initialisation.
    """

    def __init__(self, *args, **kwargs):
        self.name = "D42"
        self.order = 1
        self.A = be.array([-1.0 / 12.0, 2.0 / 3.0, 0.0, -2.0 / 3.0, 1.0 / 12.0])
        self.Ql = be.array(
            [
                [-24.0 / 17.0, 59.0 / 34.0, -4.0 / 17.0, -3.0 / 34.0, 0, 0],
                [-1.0 / 2.0, 0, 1.0 / 2.0, 0, 0, 0],
                [4.0 / 43.0, -59.0 / 86.0, 0, 59.0 / 86.0, -4.0 / 43.0, 0],
                [3.0 / 98.0, 0, -59.0 / 98.0, 0, 32.0 / 49.0, -4.0 / 49.0],
            ]
        )
        self.Qr = -self.Ql[::-1, ::-1]
        # P is the identity, H is as given above
        self.pbound = be.array([48.0 / 17])
        self.bdyRegion = self.Ql.shape
        super(D42, self).__init__(*args, **kwargs)


class D43_Tiglioetal(SBP):
    """D43 is a finite difference operator which has the SBP property.

    It is 4th order accurate in the interior and 3rd order accurate at the
    boundaries. It is from the paper, "Optimized high-order derivative and
    dissipation operators
    satisfying summation by parts, and applications in three-dimensional
    multi-block evolutions" by Diener, Dorband, Schnetter and Tiglio (DDST).

    The operator corresponds to the operator with minimum error.

    From Strand, page 75, we know that the norm is restricted full, as is also
    mentioned in DDST. It appears that, without perhaps staring at Cactus
    code, the values of the three parameters are not given. This is
    needed as the numerical values of the norm used depends on these values.

    Fortunately we do have the explicit values for Q and as such can calculate
    the correct values for the norm from them. We only need
    to calculate h_00 as the norm is restricted full. Using the notation of
    page 62
    of Strand we get,
    h_00 = 4.186595269326998 = x_1.
    """

    def __init__(self, *args, **kwargs):
        self.name = "D43_Tiglioetal"
        self.order = 1
        self.A = be.array([-1.0 / 12.0, 2.0 / 3.0, 0.0, -2.0 / 3.0, 1.0 / 12.0])

        self.Ql = be.zeros((5, 7))
        self.Ql[0, 0] = -2.09329763466349871588733
        self.Ql[0, 1] = 4.0398572053206615302160
        self.Ql[0, 2] = -3.0597858079809922953240
        self.Ql[0, 3] = 1.37319053865399486354933
        self.Ql[0, 4] = -0.25996430133016538255400
        self.Ql[0, 5] = 0
        self.Ql[0, 6] = 0
        self.Ql[1, 0] = -0.31641585285940445272297
        self.Ql[1, 1] = -0.53930788973980422327388
        self.Ql[1, 2] = 0.98517732028644343383297
        self.Ql[1, 3] = -0.05264665989297578146709
        self.Ql[1, 4] = -0.113807251750624235013258
        self.Ql[1, 5] = 0.039879767889849911803103
        self.Ql[1, 6] = -0.0028794339334846531588787
        self.Ql[2, 0] = 0.13026916185021164524452
        self.Ql[2, 1] = -0.87966858995059249256890
        self.Ql[2, 2] = 0.38609640961100070000134
        self.Ql[2, 3] = 0.31358369072435588745988
        self.Ql[2, 4] = 0.085318941913678384633511
        self.Ql[2, 5] = -0.039046615792734640274641
        self.Ql[2, 6] = 0.0034470016440805155042908
        self.Ql[3, 0] = -0.01724512193824647912172
        self.Ql[3, 1] = 0.16272288227127504381134
        self.Ql[3, 2] = -0.81349810248648813029217
        self.Ql[3, 3] = 0.13833269266479833215645
        self.Ql[3, 4] = 0.59743854328548053399616
        self.Ql[3, 5] = -0.066026434346299887619324
        self.Ql[3, 6] = -0.0017244594505194129307249
        self.Ql[4, 0] = -0.00883569468552192965061
        self.Ql[4, 1] = 0.03056074759203203857284
        self.Ql[4, 2] = 0.05021168274530854232278
        self.Ql[4, 3] = -0.66307364652444929534068
        self.Ql[4, 4] = 0.014878787464005191116088
        self.Ql[4, 5] = 0.65882706381707471953820
        self.Ql[4, 6] = -0.082568940408449266558615

        self.Qr = -self.Ql[::-1, ::-1]
        # P is the identity.
        self.pbound = be.array([4.186595370392226897362216859769846226369])
        super(D43_Tiglioetal, self).__init__(*args, **kwargs)


class D43_CNG(SBP):
    """
    An SBP operator that is fourth order accurate in the interior and
    third order accurate on the boundary.

    Taken from "A stable and conservative interface treatment of arbitrary
    spatial accuracy", the matrix H is the identity.
    Note that the P[1,3] entry of the norm matrix given in the paper is wrong!
    """

    def __init__(self, *args, **kwargs):
        self.r1 = -(2177.0 * math.sqrt(295369.0) - 1166427.0) / (25488.0)
        self.r2 = (66195.0 * math.sqrt(53.0 * 5573.0) - 35909375.0) / 101952.0
        self.A = be.array([-1.0 / 12.0, 2.0 / 3.0, 0.0, -2.0 / 3.0, 1.0 / 12.0])
        self.name = "D43_CNG"

        self.order = 1
        a = self.r1
        b = self.r2
        Q = be.zeros((4, 7))

        Q[0, 0] = -0.5
        Q[0, 1] = -(864.0 * b + 6480 * a + 305) / 4320.0
        Q[0, 2] = (216 * b + 1620 * a + 725) / 540.0
        Q[0, 3] = -(864 * b + 6480 * a + 3335) / 4320

        self.g00 = -1  # Q[0,0]
        self.gnn = 1  # -Q[0,0]

        Q[1, 0] = -Q[0, 1]
        Q[1, 1] = 0.0
        Q[1, 2] = -(864.0 * b + 6480 * a + 2315) / 1440.0
        Q[1, 3] = (108 * b + 810 * a + 415) / 270

        Q[2, 0] = -Q[0, 2]
        Q[2, 1] = -Q[1, 2]
        Q[2, 2] = 0.0
        Q[2, 3] = -(864 * b + 6480 * a + 785) / 4320

        Q[3, 0] = -Q[0, 3]
        Q[3, 1] = -Q[1, 3]
        Q[3, 2] = -Q[2, 3]
        Q[3, 3] = 0.0

        Q[2, 4] = -1.0 / 12.0
        Q[3, 5] = -1.0 / 12.0
        Q[3, 4] = 8.0 / 12.0

        P = be.zeros((4, 4))
        P[0, 0] = -(216 * b + 2160 * a - 2125) / (12960)
        P[0, 1] = (81 * b + 675 * a + 415) / 540
        P[0, 2] = -(72 * b + 720 * a + 445) / (1440)
        P[0, 3] = -(108 * b + 756 * a + 421) / 1296

        P[1, 0] = P[0, 1]
        P[1, 1] = -(4104 * b + 32400 * a + 11225) / 4320
        P[1, 2] = (1836 * b + 14580 * a + 7295) / 2160
        P[1, 3] = -(216 * b + 2160 * a + 655) / (4320)

        P[2, 0] = P[0, 2]
        P[2, 1] = P[1, 2]
        P[2, 2] = -(4104 * b + 32400 * a + 12785) / 4320
        P[2, 3] = (81 * b + 675 * a + 335) / (540)

        P[3, 0] = P[0, 3]
        P[3, 1] = P[1, 3]
        P[3, 2] = P[2, 3]
        P[3, 3] = -(216 * b + 2160 * a - 12085) / (12960)

        Pinv = be.linalg.inv(P)
        self.pbound = Pinv[:, 0]
        self.Ql = be.dot(Pinv, Q)
        self.Qr = -self.Ql[::-1, ::-1]

        # Note that the operation Ql[::-1,::-1]
        # is not the transpose. The instructions in
        # the relevant paper are misleading.
        # >>> a = be.array([[1,2,3],[4,5,6],[7,8,9]])
        # >>> a
        # array([[1, 2, 3],
        #       [4, 5, 6],
        #       [7, 8, 9]])
        # >>> a[::-1,::-1]
        # array([[9, 8, 7],
        #       [6, 5, 4],
        #       [3, 2, 1]])

        super(D43_CNG, self).__init__(*args, **kwargs)


class D43_Strand(SBP):
    """An SBP operator that is fourth order accurate in the interior and
    third order accurate on the boundary.

    See page 66 of Strand's paper, "Summation by parts finite difference
    approximations for first derivatives". The norm in this case
    is restricted full and we have h00 = 3./11.

    Note the terms introduced in the initialisation code.
    """

    def __init__(self, *args, **kwargs):
        self.A = be.array([-1.0 / 12.0, 2.0 / 3.0, 0.0, -2.0 / 3.0, 1.0 / 12.0])
        self.name = "D43_Strand"

        self.order = 1
        Q = be.mat(be.zeros((5, 7)))
        Q[0, 0] = -11.0 / 6
        Q[0, 1] = 3.0
        Q[0, 2] = -3.0 / 2
        Q[0, 3] = 1.0 / 3
        Q[0, 4] = 0
        Q[0, 5] = 0
        Q[0, 6] = 0
        Q[1, 0] = -0.389422071485311842975177265599
        Q[1, 1] = -0.269537639034869460503559633378
        Q[1, 2] = 0.639037937659262938432677856167
        Q[1, 3] = 0.0943327360845463774750968877551
        Q[1, 4] = -0.0805183715808445133581024825052
        Q[1, 5] = 0.00610740835721650092906463755990
        Q[1, 6] = 0
        Q[2, 0] = 0.111249966676253227197631191911
        Q[2, 1] = -0.786153109432785509340645292042
        Q[2, 2] = 0.198779437635276432052935915726
        Q[2, 3] = 0.508080676928351487908752085966
        Q[2, 4] = -0.0241370624126563706018867104954
        Q[2, 5] = -0.00781990939443926721678719106507
        Q[2, 6] = 0
        Q[3, 0] = 0.0190512060948850190478223587421
        Q[3, 1] = 0.0269311042007326141816664674713
        Q[3, 2] = -0.633860292039252305642283500163
        Q[3, 3] = 0.0517726709186493664626888177616
        Q[3, 4] = 0.592764606048964306931634491846
        Q[3, 5] = -0.0543688142698406758774679261355
        Q[3, 6] = -0.00229048095413832510406070952285
        Q[4, 0] = -0.00249870649542362738624804675220
        Q[4, 1] = 0.00546392445304455008494236684036
        Q[4, 2] = 0.0870248056190193154450416111553
        Q[4, 3] = -0.686097670431383548237962511314
        Q[4, 4] = 0.0189855304809436619879348998899
        Q[4, 5] = 0.659895344563505072850627735853
        Q[4, 6] = -0.0827732281897054247443360556719

        self.Ql = Q
        self.Qr = -self.Ql[::-1, ::-1]
        self.pbound = be.array([11.0 / 3])

        super(D43_Strand, self).__init__(*args, **kwargs)


################################################################################
# Higher order accurate differential operators
################################################################################


class D65_min_err(SBP):
    """D65_min_err is a first order derivative according to Diener et al,
    which has minimised error. Coefficients taken from the source file
    of the paper in the arxive.org repository.

    As above the norm is not given. We can calculate h_00 by noting
    that h is restricted full and that {h^(-1)q}_{00} = -1/2.
    The result is, to 15 decimal places,
    h_{00} = 4.930709842221048

    """

    def __init__(self, *args, **kwargs):
        self.name = "D65"
        self.A = be.array(
            [
                1.0 / 60.0,
                -3.0 / 20.0,
                3.0 / 4.0,
                0.0,
                -3.0 / 4.0,
                3.0 / 20.0,
                -1.0 / 60.0,
            ]
        )
        self.Ql = be.mat(be.zeros((7, 10)))

        self.Ql[0, 0] = -2.465354921110524023660777656111276003457
        self.Ql[0, 1] = 6.092129526663144141964665936667656020742
        self.Ql[0, 2] = -7.730323816657860354911664841669140051855
        self.Ql[0, 3] = 6.973765088877147139882219788892186735807
        self.Ql[0, 4] = -3.980323816657860354911664841669140051855
        self.Ql[0, 5] = 1.292129526663144141964665936667656020742
        self.Ql[0, 6] = -0.1820215877771906903274443227779426701237
        self.Ql[0, 7] = 0
        self.Ql[0, 8] = 0
        self.Ql[0, 9] = 0

        self.Ql[1, 0] = -0.2234725650784319828746535134412736890421
        self.Ql[1, 1] = -0.9329308121107134563129925525068570679651
        self.Ql[1, 2] = 1.586820596545839371759081303802027231274
        self.Ql[1, 3] = -0.3647002340377160216914505558624668821400
        self.Ql[1, 4] = -0.2666957784872806143914117440166232718819
        self.Ql[1, 5] = 0.3112949048634705032101261273629794071371
        self.Ql[1, 6] = -0.1404504214762266650000768489896480092493
        self.Ql[1, 7] = 0.03488568514730479833596013512958238764128
        self.Ql[1, 8] = -0.004964021886392518344179263072091597647654
        self.Ql[1, 9] = 0.0002126465201465853095969115943714918742904

        self.Ql[2, 0] = 0.1582216737061633151406179477554921935333
        self.Ql[2, 1] = -1.137049298003377811733609086574457439398
        self.Ql[2, 2] = 1.212364522932578587741649981040340946798
        self.Ql[2, 3] = -0.9562288729513894906148167047868730813830
        self.Ql[2, 4] = 1.066548057336766350478498057851678826640
        self.Ql[2, 5] = -0.3478788551267041838265477441805600110467
        self.Ql[2, 6] = -0.03133923293520187620333693909408071632123
        self.Ql[2, 7] = 0.04098845955755862691072597869183962277781
        self.Ql[2, 8] = -0.005963188634687155197078928402509551508436
        self.Ql[2, 9] = 0.0003367341182936373038974376991292099082999

        self.Ql[3, 0] = 0.02915734641890708196910927068736798144670
        self.Ql[3, 1] = -0.1169665089768926152768236581512624861308
        self.Ql[3, 2] = -0.1112219092451476301503253995474190870412
        self.Ql[3, 3] = -0.7924486261248032107393766820001361351677
        self.Ql[3, 4] = 1.266650704820613624987450232358951199911
        self.Ql[3, 5] = -0.2899273290506621673153239836530375587273
        self.Ql[3, 6] = 0.002515684257201926199329020583484434062150
        self.Ql[3, 7] = 0.01329713961871764653006682056620518602804
        self.Ql[3, 8] = -0.001124464399630667352932212208930962568134
        self.Ql[3, 9] = 0.00006796268169601114882659136477742818715059

        self.Ql[4, 0] = -0.04582150000326981674750984653096293434777
        self.Ql[4, 1] = 0.2240986548857151482718685516611524323427
        self.Ql[4, 2] = -0.3246718493011818141660859125588209338018
        self.Ql[4, 3] = -0.3929792921782506986152017485694441380503
        self.Ql[4, 4] = 0.1166355818729375628072830916953646214341
        self.Ql[4, 5] = 0.3449626905957060254933930895775644438105
        self.Ql[4, 6] = 0.1430419813354607083034935179267283951745
        self.Ql[4, 7] = -0.07764802499372607792980458731991885121073
        self.Ql[4, 8] = 0.01332439335504217034559288889042994978834
        self.Ql[4, 9] = -0.0009426355684332077630290447720929851395193

        self.Ql[5, 0] = 0.003172814452954821196677290327889903944225
        self.Ql[5, 1] = 0.00001061446045061551877105554145609103530766
        self.Ql[5, 2] = -0.08747763580209736614983637747947172321794
        self.Ql[5, 3] = 0.3975827322299876034907453299884380895682
        self.Ql[5, 4] = -1.148835072393422871630425744497391344782
        self.Ql[5, 5] = 0.3583006649535242306065761818925080902380
        self.Ql[5, 6] = 0.5647665154270147564019144982190032455071
        self.Ql[5, 7] = -0.09698196887272109736153117076061707705561
        self.Ql[5, 8] = 0.008843905091972988427261446924164441884143
        self.Ql[5, 9] = 0.0006174304523363194998474898440202828786385

        self.Ql[6, 0] = -0.008639107540858839028043929986084287776394
        self.Ql[6, 1] = 0.04722773954485212324714352753530343274219
        self.Ql[6, 2] = -0.1008747537650261142294540111407681552350
        self.Ql[6, 3] = 0.08043834953845218736895768965086958762389
        self.Ql[6, 4] = 0.1295138674713300902982857323205417604553
        self.Ql[6, 5] = -0.7909424166489541737614153656634872155367
        self.Ql[6, 6] = 0.03807866847647628589685997987877954466259
        self.Ql[6, 7] = 0.7367055699548196242687865288427927434250
        self.Ql[6, 8] = -0.1480235854665196220062411065981933720158
        self.Ql[6, 9] = 0.01651566843542843794512095516024596165494

        self.Qr = -self.Ql[::-1, ::-1]
        self.order = 1
        self.pbound = be.array([1 / 4.930709842221048])
        super(D65_min_err, self).__init__(*args, **kwargs)


################################################################################
# Second Derivative - Second Order SBP operators
################################################################################


class D43_2_CNG(SBP):
    """
    A second derivative SBP operator with fourth order interal accuracy and
    third order accuracy on the boundary.

    Taken from "A stable and conservative interface treatment of arbitrary
    spatial accuracy", the matrix H is the identity.
    """

    def __init__(self, *args, **kwargs):
        self.r1 = -(2177.0 * math.sqrt(295369.0) - 1166427.0) / (25488.0)
        self.r2 = (66195.0 * math.sqrt(53.0 * 5573.0) - 35909375.0) / 101952.0
        self.A = be.array(
            [-1.0 / 12.0, 16.0 / 12.0, -30.0 / 12.0, 16.0 / 12.0, -1.0 / 12.0]
        )
        self.name = "D43_2_CNG"

        self.order = 2
        a = self.r1
        b = self.r2
        Ql = be.zeros((3, 5))

        Ql[0, 0] = 35.0 / 12
        Ql[0, 1] = -26.0 / 3
        Ql[0, 2] = 19.0 / 2
        Ql[0, 3] = -14.0 / 3
        Ql[0, 4] = 11.0 / 12

        self.g00 = 1
        self.gnn = 1

        Ql[1, 0] = 11.0 / 12
        Ql[1, 1] = -5.0 / 3
        Ql[1, 2] = 1.0 / 2
        Ql[1, 3] = 1.0 / 3
        Ql[1, 4] = -1.0 / 12

        Ql[2, 0] = -1.0 / 12
        Ql[2, 1] = 16.0 / 12
        Ql[2, 2] = -30.0 / 12
        Ql[2, 3] = 16.0 / 12
        Ql[2, 4] = -1.0 / 12

        P = be.zeros((4, 4))
        P[0, 0] = -(216 * b + 2160 * a - 2125) / (12960)
        P[0, 1] = (81 * b + 675 * a + 415) / 540
        P[0, 2] = -(72 * b + 720 * a + 445) / (1440)
        P[0, 3] = -(108 * b + 756 * a + 421) / 1296

        P[1, 0] = P[0, 1]
        P[1, 1] = -(4104 * b + 32400 * a + 11225) / 4320
        P[1, 2] = (1836 * b + 14580 * a + 7295) / 2160
        P[1, 3] = -(216 * b + 2160 * a + 665) / (4320)

        P[2, 0] = P[0, 2]
        P[2, 1] = P[1, 2]
        P[2, 2] = -(4104 * b + 32400 * a + 12785) / 4320
        P[2, 3] = (81 * b + 675 * a + 335) / (540)

        P[3, 0] = P[0, 3]
        P[3, 1] = P[1, 3]
        P[3, 2] = P[2, 3]
        P[3, 3] = -(216 * b + 2160 * a - 12085) / (12960)

        Pinv = be.linalg.inv(P)
        self.pbound = Pinv[:, 0]
        self.Ql = Ql
        self.Qr = self.Ql[::-1, ::-1]

        super(D43_2_CNG, self).__init__(*args, **kwargs)
