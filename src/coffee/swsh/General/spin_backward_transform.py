from coffee.settings import be
import scipy.fft as fft

def Compute_Gmn(alm, Gmn, DD, s, lmax, smax, Sampling_on_torus_theta, Sampling_on_torus_phi):
    Tolerance = 1e-13

    def Sup(a, b):
        if a <= b:
            return b
        else:
            return a

    def pindex(l, i):
        if i >= 0:
            return i
        else:
            return 2 * l + 1 + i

    abss = abs(s)

    # case m=n=0
    for l in range(abss, lmax + 1, 2):
        if Tolerance < abs(DD[1][l][pindex(smax, s)][0][0]):
            Gmn[0] += DD[1][l][pindex(smax, s)][0][0] * alm[l * (l + 1)]

    # case m=0
    n = 0
    while n <= lmax:
        n += 1
        l = Sup(abss, n)
        while l <= lmax:
            l_centro = l * (l + 1)
            if Tolerance < abs(DD[1][l][pindex(smax, s)][0][n]):
                Gmn[n] += DD[1][l][pindex(smax, s)][0][n] * alm[l_centro + n]
                Gmn[Sampling_on_torus_phi - n] += DD[1][l][pindex(smax, s)][0][2 * l + 1 - n] * alm[l_centro - n]
            l += 2

    # case n=0
    m = 0
    Power1 = (-1.) ** (-s)
    while m <= lmax:
        m += 1
        l = m
        while l <= lmax:
            l_centro = l * (l + 1)
            if Tolerance < abs(DD[1][l][pindex(smax, s)][m][0]):
                Gmn[m * Sampling_on_torus_phi] += DD[1][l][pindex(smax, s)][m][0] * alm[l_centro]
                Gmn[(Sampling_on_torus_theta - m) * Sampling_on_torus_phi] = Power1 * Gmn[m * Sampling_on_torus_phi]
            l += 2

    # other cases
    m = 0
    while m <= lmax:
        m += 1
        n = 0
        while n <= lmax:
            n += 1
            l = Sup(Sup(abss, n), Sup(abss, m))
            while l <= lmax:
                print(m, n, l)
                if Tolerance < abs(DD[1][l][pindex(smax, s)][m][n]):
                    Gmn[m * Sampling_on_torus_phi + n] += DD[1][l][pindex(smax, s)][m][n] * alm[l * (l + 1) + n]
                    Gmn[m * Sampling_on_torus_phi + (Sampling_on_torus_phi - n)] += DD[1][l][pindex(smax, s)][m][
                        2 * l + 1 - n] * alm[l * (l + 1) - n]
                l += 1
            Gmn[(Sampling_on_torus_theta - m) * Sampling_on_torus_phi + n] = (-1.) ** (-s - n) * Gmn[m * Sampling_on_torus_phi + n]
            Gmn[(Sampling_on_torus_theta - m) * Sampling_on_torus_phi + (Sampling_on_torus_phi - n)] = (-1.) ** (-s + n) * Gmn[m * Sampling_on_torus_phi + (Sampling_on_torus_phi - n)]


def spin_backward_transform(f, alm, s, lmax, smax, Ntheta, Nphi, DD):
    Sampling_on_torus_theta = 2 * (Ntheta - 1)
    Sampling_on_torus_phi = Nphi

    if 2 * lmax + 1 > Sampling_on_torus_phi:
        print("Number of points in phi are less than 2*lmax+1.")
        print("Abort swsh-backward-transform!")
        return

    if 2 * lmax + 1 > Sampling_on_torus_theta:
        print("Number of points in theta are less than lmax+3/2.")
        print("Abort swsh-backward-transform!")
        return

    Gmn = be.zeros((Sampling_on_torus_theta * Sampling_on_torus_phi,), dtype=be.complex128)

    Compute_Gmn(alm, Gmn, DD, s, lmax, smax, Sampling_on_torus_theta, Sampling_on_torus_phi)

    Gmn_reshaped = Gmn.reshape((Sampling_on_torus_theta, \
                                Sampling_on_torus_phi))

    f_torus = fft.ifft2(Gmn_reshaped) * Sampling_on_torus_theta*Sampling_on_torus_phi

    f[:] = f_torus[:Ntheta, :Nphi]
