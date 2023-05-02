# import standard modules
import numpy as np

# import our modules
from coffee.tslices import tslices
from coffee.system import System
from coffee.backend import backend as be

be.set_backend("numpy")


class TwoDadvection(System):
    def timestep(self, tslice):
        ssizes = tslice.domain.step_sizes
        spatial_divisor = (1 / ssizes[0]) + (1 / ssizes[1])
        dt = self.CFL / spatial_divisor
        return dt

    ############################################################################
    # Constructor
    ############################################################################
    def __init__(self, xdirec, ydirec, Dx, Dy, CFL, tau=None, equation_coords="Polar"):
        super(TwoDadvection, self).__init__()
        self.CFL = CFL
        self.xcoef = xdirec
        self.ycoef = ydirec
        self.Dx = Dx
        self.Dy = Dy
        self.numvar = 1
        self.tau = tau
        self.equation_coords = equation_coords
        self.name = """<TwoDadvection xdirec = %f, ydirec = %f, Dx = %s, 
        Dy = %s, CLF = %f, tau = %s>""" % (
            xdirec,
            ydirec,
            Dx.name,
            Dy.name,
            CFL,
            repr(tau),
        )

    ############################################################################
    # Configuration for initial conditions and boundary conditions.
    ############################################################################
    def initial_data(self, t0, r):
        return self.centralBump(t0, r)

    def boundary(self, t, Psi):
        return be.zeros_like(Psi.data[0])

    def first_right(self, t, Psi):
        return be.zeros_like(Psi.domain.axes[0])

    def first_left(self, t, Psi):
        return (0.0, 0.0)

    ############################################################################
    # Evolution Routine
    ############################################################################
    def evaluate(self, t, Psi, intStep=None):
        # Define useful variables
        (f0,) = Psi.data

        x = Psi.domain.axes[0]
        y = Psi.domain.axes[1]
        dx = Psi.domain.step_sizes[0]
        dy = Psi.domain.step_sizes[1]
        tau = self.tau

        ########################################################################
        # Calculate derivatives and impose boundary conditions
        ########################################################################
        Dxf = be.apply_along_axis(lambda x: self.Dx(x, dx), 0, f0)

        Dyf = be.apply_along_axis(lambda y: self.Dy(y, dy), 1, f0)

        ########################################################################
        # Impose boundary conditions
        ########################################################################

        # implementation follows Carpenter et al.
        # using the SAT method
        # at the boundaries we need boundary conditions
        # implemented as penalty terms the objects in diffop.py know how to
        # do this.
        #
        # tau is the penalty parameter and will need to take on different
        # values depending on the operator.
        #

        pt_x_r = self.Dx.penalty_boundary(dx, "right")
        pt_x_r_shape = pt_x_r.size
        pt_x_l = self.Dx.penalty_boundary(dx, "left")
        pt_x_l_shape = pt_x_l.size

        pt_y_r = self.Dy.penalty_boundary(dy, "right")
        pt_y_r_shape = pt_y_r.size
        pt_y_l = self.Dy.penalty_boundary(dy, "left")
        pt_y_l_shape = pt_y_l.size

        # First do internal boundaries
        _, b_values = Psi.communicate()  # compare to OneDAdvection for an
        # alternative way to handle this.
        for d_slice, data in b_values:
            # the calculation of sigma constants is taken from Carpenter,
            # Nordstorm and Gottlieb. Note that in this paper the metric H is
            # always set to the identity matrix. Beware: in some presentations
            # of SBP operators it is not the identitiy. This is accounted for
            # in the calculation of pt below.
            # I think that this paper implicitly assumes that 'a' is positive
            # hence the difference for psi4 from the calculations given in
            # the paper. This change accounts for the negative eigenvalue
            # associated to psi4.
            # Note that sigma3 = sigma1 - eigenvalue_on_boundary, at least when
            # the eigenvalue is positive. For negative eigenvalue it seems to me
            # that the roles of sigma3 and sigma1 are reversed.
            x_chara = self.xcoef
            y_chara = self.ycoef
            if x_chara > 0:
                sigma3x = 0.25
                sigma1x = sigma3x - 1
            else:
                sigma1x = 0.25
                sigma3x = sigma1x - 1
            if y_chara > 0:
                sigma3y = 0.25
                sigma1y = sigma3y - 1
            else:
                sigma1y = 0.25
                sigma3y = sigma1y - 1

            if d_slice[1] == slice(-1, None, None):
                Dxf[-pt_x_r_shape:] += (
                    sigma1x * x_chara * pt_x_r * (f0[d_slice[1:]] - data[0])
                )
            elif d_slice[1] == slice(None, 1, None):
                Dxf[:pt_x_l_shape] += (
                    sigma3x * x_chara * pt_x_l * (f0[d_slice[1:]] - data[0])
                )
            elif d_slice[1] == slice(None, None, None):
                if d_slice[2] == slice(-1, None, None):
                    Dyf[:, -pt_y_r_shape:] += (
                        sigma1y * y_chara * pt_y_r * (f0[d_slice[1:]] - data[0])
                    )
                elif d_slice[2] == slice(None, 1, None):
                    Dyf[:, :pt_y_l_shape] += (
                        sigma3y * y_chara * pt_y_l * (f0[d_slice[1:]] - data[0])
                    )

        # Now do the external boundaries
        b_data = Psi.external_slices()
        for dim, direction, d_slice in b_data:
            d_slice = d_slice[1:]
            if dim == 0:
                if self.xcoef > 0 and direction == 1:
                    Dxf[-pt_x_r_shape:] -= (
                        tau
                        * self.xcoef
                        * (f0[d_slice] - self.boundary(t, Psi)[d_slice])
                        * pt_x_r
                    )
                if self.xcoef < 0 and direction == -1:
                    Dxf[:pt_x_l_shape] += (
                        tau
                        * self.xcoef
                        * (f0[d_slice] - self.boundary(t, Psi)[d_slice])
                        * pt_x_l
                    )
            elif dim == 1:
                if self.ycoef > 0 and direction == 1:
                    Dyf[:, -pt_y_r_shape:] -= (
                        tau
                        * self.ycoef
                        * (f0[d_slice] - self.boundary(t, Psi)[d_slice])
                        * pt_y_r
                    )
                if self.ycoef < 0 and direction == -1:
                    Dyf[:, :pt_y_l_shape] += (
                        tau
                        * self.ycoef
                        * (f0[d_slice] - self.boundary(t, Psi)[d_slice])
                        * pt_y_l
                    )

        Dtf = self.xcoef * Dxf + self.ycoef * Dyf

        print(("t =", t))

        # now all time derivatives are computed
        # package them into a time slice and return
        rtslice = tslices.TimeSlice([Dtf], Psi.domain, time=t)
        return rtslice

    ############################################################################
    # Boundary functions
    ############################################################################
    def dirichlet_boundary(self, u, intStep=None):
        u.fields[0][:, 0] = u.fields[0][:, -1]
        return u

    ############################################################################
    # Initial Value Routines
    ############################################################################
    def centralBump(self, t0, grid):
        from mpi4py import MPI

        rank = MPI.COMM_WORLD.rank
        r, phi = grid.axes[0], grid.axes[1]
        r3ind = int(r.shape[0] / 6)
        r3val = r[r3ind]
        r6val = r[6 * r3ind]
        phi3ind = int(phi.shape[0] / 5)
        phi3val = phi[phi3ind]
        phi6val = phi[4 * phi3ind]
        rmid = int(r.shape[0] / 2)
        phimid = int(phi.shape[0] / 2)
        r_mesh, phi_mesh = grid.meshes

        rv = be.exp(-20 * (r_mesh - r[rmid]) ** 2) * be.exp(
            -5 * (phi_mesh - phi[phimid]) ** 2
        )
        if rank == 0:
            rtslice = tslices.TimeSlice([rv], grid, t0)
        else:
            rtslice = tslices.TimeSlice([be.zeros_like(rv)], grid, t0)
        return rtslice
