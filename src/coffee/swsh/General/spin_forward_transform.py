from coffee.settings import be
import scipy.fft as fft

def pindex(l, i):
    if i >= 0:
        return i
    else:
        return 2 * l + 1 + i

def Compute_Imn_by_extending_f_to_torus(f, Imn, W, s, Sampling_on_torus_theta, Sampling_on_torus_phi):
    Norm = 2 * be.pi / (Sampling_on_torus_theta * Sampling_on_torus_phi)

    F = be.zeros(Sampling_on_torus_theta*Sampling_on_torus_phi, dtype=be.complex128)

    Ntheta = Sampling_on_torus_theta // 2 + 1

    power_s = (-1.0)**s

    for itheta in range(Ntheta):
        for iphi in range(Sampling_on_torus_phi):
            opp_iphi = (iphi + Sampling_on_torus_phi // 2) % Sampling_on_torus_phi

            F[itheta * Sampling_on_torus_phi  + iphi] = be.real(W[itheta]) * f[itheta * Sampling_on_torus_phi + iphi] * Norm

            if itheta > 0 and itheta < Ntheta:
                F[(Sampling_on_torus_theta - itheta)*Sampling_on_torus_phi + opp_iphi] = power_s * be.real(W[Sampling_on_torus_theta - itheta]) * f[itheta*Sampling_on_torus_phi + iphi] * Norm

    F_reshaped = be.reshape(F, (Sampling_on_torus_theta, Sampling_on_torus_phi))
    F_out = fft.fft2(F_reshaped, norm='backward')
    F_out_reshaped = be.reshape(F_out,Sampling_on_torus_theta*Sampling_on_torus_phi)
    Imn[:] = F_out_reshaped[:]

def Compute_Jmn(f, Jmn, W, s, lmax, Sampling_on_torus_theta, Sampling_on_torus_phi):
    Imn = be.zeros((Sampling_on_torus_theta*Sampling_on_torus_phi), dtype=be.complex128)

    Compute_Imn_by_extending_f_to_torus(f, Imn, W, s, Sampling_on_torus_theta, Sampling_on_torus_phi)

    Nn = 2 * lmax + 1
    Nm = lmax

    power_s_plus_n = (-1.0)**(s + be.arange(Nm + 1))
    power_s_minus_n = (-1.0)**(s - be.arange(Nm + 1))

    for n in range(Nm + 1):
        for m in range(Nm + 1):
            if m == 0:
                Jmn[n] = Imn[n]
                if n > 0 and n < lmax:
                    Jmn[Nn - n] = Imn[Sampling_on_torus_phi - n]
            else:
                Jmn[m * Nn + n] = Imn[m * Sampling_on_torus_phi + n] + power_s_plus_n[n] * Imn[(Sampling_on_torus_theta - m) * Sampling_on_torus_phi + n]
                if n > 0 and n <= lmax:
                    Jmn[m * Nn + (Nn - n)] = Imn[m * Sampling_on_torus_phi + (Sampling_on_torus_phi - n)] + power_s_minus_n[n] * Imn[(Sampling_on_torus_theta - m) * Sampling_on_torus_phi + (Sampling_on_torus_phi - n)]

def spin_forward_transform(f, aln, s, lmax, smax, Ntheta, Nphi, DD, W):
    Sampling_on_torus_theta = 2 * (Ntheta - 1)
    Sampling_on_torus_phi = Nphi

    if 2 * lmax + 1 > Sampling_on_torus_phi:
        print("Number of points in phi are less than 2*lmax+1 (read the manual!)")
        print("Abort swsh-forward-transform")
        return

    if 2 * lmax + 1 > Sampling_on_torus_theta:
        print("Number of points in theta are less than lmax+3/2 (read the manual!)")
        print("Abort swsh-forward-transform")
        return

    Jmn = be.zeros((lmax + 1) * (2 * lmax + 1), dtype=be.complex128)

    Compute_Jmn(f, Jmn, W, s, lmax, Sampling_on_torus_theta, Sampling_on_torus_phi)

    s_index = pindex(smax, s)
    number_of_n = (2*lmax+1)
    Tolerance = 1e-13
    for l in range(abs(s), lmax + 1):
        l_centro = l * (l + 1)
    
        for m in range(l % 2, l + 1, 2):
            if Tolerance < abs(DD[0, l, s_index, m, 0]):
                aln[l_centro] += DD[0, l, s_index, m, 0] * Jmn[m * number_of_n]
    
        for n in range(1, l + 1):
            for m in range(l + 1):
                if Tolerance < abs(DD[0, l, s_index, m, n]):
                    aln[l_centro + n] += DD[0, l, s_index, m, n] * Jmn[m * number_of_n + n]
                    aln[l_centro - n] += DD[0, l, s_index, m, 2 * l + 1 - n] * Jmn[m * number_of_n + (number_of_n - n)]