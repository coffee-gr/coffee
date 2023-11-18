from coffee.settings import be

from coffee.swsh.General import spin_backward_transform as sb
from coffee.swsh.General import spin_forward_transform as sf

def spin_ethU(f, alm, s, lmax, smax, Ntheta, Nphi, DD, W):
    sf.spin_forward_transform(f, alm, s, lmax, smax, Ntheta, Nphi, DD, W)
    Nlm = lmax * (lmax + 2) + 1
    m = 0
    l = 0
    i = 0

    while i < Nlm:
        square_root_terms_up = -be.sqrt((l - s) * (l + s + 1))

        while i <= m:
            if l < abs(s):
                alm[i] = 0.0
            else:
                alm[i] = square_root_terms_up * alm[i]
            i += 1

        l += 1
        m = (2 * l + 1) + m

    s += 1
    f = f.reshape((Ntheta, Nphi))
    sb.spin_backward_transform(f, alm, s, lmax, smax, Ntheta, Nphi, DD)

def spin_ethD(f, alm, s, lmax, smax, Ntheta, Nphi, DD, W):
    sf.spin_forward_transform(f, alm, s, lmax, smax, Ntheta, Nphi, DD, W)
    Nlm = lmax * (lmax + 2) + 1
    m = 0
    l = 0
    i = 0

    while i < Nlm:
        square_root_terms_down = be.sqrt((l + s) * (l - s + 1))

        while i <= m:
            if l < abs(s):
                alm[i] = 0.0
            else:
                alm[i] = square_root_terms_down * alm[i]
            i += 1

        l += 1
        m = (2 * l + 1) + m

    s -= 1
    f = f.reshape((Ntheta, Nphi))
    sb.spin_backward_transform(f, alm, s, lmax, smax, Ntheta, Nphi, DD)