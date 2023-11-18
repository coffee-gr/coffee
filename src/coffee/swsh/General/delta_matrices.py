from coffee.settings import be

def Compute_Delta_Matrices(D, lmax):
    D[0][0][0] = 1.0

    for l in range(1, lmax + 1):
        for m2 in range(0, l + 1):
            if m2 == 0: 
                D[l][l][0] = -be.sqrt((2. * l - 1) / (2. * l)) * D[l - 1][l - 1][0]
            else: 
                D[l][l][m2] = be.sqrt(((l / 2.) * (2. * l - 1)) / ((l + m2) * (l + m2 - 1.))) * D[l - 1][l - 1][m2 - 1]
            for m1 in range(l - 1, -1, -1):
                if m1 == l - 1:
                    D[l][l - 1][m2] = ((2. * m2) / be.sqrt((l - m1) * (l + m1 + 1.))) * D[l][l][m2]
                else:
                    D[l][m1][m2] = ((2. * m2) / be.sqrt((l - m1) * (l + m1 + 1.))) * D[l][m1 + 1][m2] - be.sqrt(
                        ((l - m1 - 1.) * (l + m1 + 2.)) / ((l - m1) * (l + m1 + 1.))) * D[l][m1 + 2][m2]

def Precompute_Delta_Matrices(lmax):
    D = be.zeros((lmax + 1, lmax + 1, lmax + 1))
    Compute_Delta_Matrices(D, lmax)
    return D

def Dlmn(D, l, m, n):
    if m >= 0 and n >= 0:
        return D[l][m][n]
    elif m >= 0 and n < 0:
        return (-1) ** (l + m) * D[l][m][abs(n)]
    elif m < 0 and n >= 0:
        return (-1) ** (l + n) * D[l][abs(m)][n]
    else:
        return (-1) ** (2 * l - m + n) * D[l][abs(m)][abs(n)]

def FACTOR1(s, l, m, n):
    return (-1.0) ** (s + l + m) * (1j) ** (s + n) * be.sqrt((2. * l + 1) / (4. * be.pi))

def FACTOR2(s, l, m, n):
    return (-1.0) ** (2 * l - 2 * m + s + n)

def pindex(l, i):
    if i >= 0:
        return i
    else:
        return 2 * l + 1 + i

def Precompute_Delta_Delta_Matrices(lmax, smax):
    D = be.zeros((lmax + 1, lmax + 1, lmax + 1), dtype=be.float64)
    Compute_Delta_Matrices(D, lmax)

    DD = be.zeros((2, lmax + 1, 2 * smax + 1, 2 * lmax + 1, 2 * lmax + 1), dtype=be.complex128)

    for p in range(2):
        # DD[0] for the forward transform
        # DD[1] for the backward transform
        DD[p] = be.zeros((lmax + 1, 2 * smax + 1, 2 * lmax + 1, 2 * lmax + 1), dtype=be.complex128)
        for i in range(lmax + 1):
            DD[p][i] = be.zeros((2 * smax + 1, 2 * lmax + 1, 2 * lmax + 1), dtype=be.complex128)
            for j in range(2 * smax + 1):
                DD[p][i][j] = be.zeros((2 * lmax + 1, 2 * lmax + 1), dtype=be.complex128)

    l, s, m, n = 0, 0, 0, 0

    Tolerance = 1e-14

    def Inf(a, b):
        if a <= b:
            return a
        else:
            return b

    m_counter = 0
    sm_counter = 0

    for l in range(lmax + 1):
        for s in range(Inf(smax, l) + 1):
            for n in range(l + 1):
                for m in range(l + 1):
                    sm_counter += 1

                    if Tolerance < abs(Dlmn(D, l, m, n) * Dlmn(D, l, m, s)):
                        # (s+, m+, n+)
                        DD[0][l][pindex(smax, s)][pindex(l, m)][pindex(l, n)] = (
                            FACTOR1(s, l, m, n) * Dlmn(D, l, m, n) * Dlmn(D, l, m, s)
                        )
                        DD[1][l][pindex(smax, s)][pindex(l, m)][pindex(l, n)] = (
                            FACTOR2(s, l, m, n) * DD[0][l][pindex(smax, s)][pindex(l, m)][pindex(l, n)]
                        )

                        # (s+, m+, n-)
                        if n > 0:
                            DD[0][l][pindex(smax, s)][pindex(l, m)][pindex(l, -n)] = (
                                FACTOR1(s, l, m, -n) * Dlmn(D, l, m, -n) * Dlmn(D, l, m, s)
                            )
                            DD[1][l][pindex(smax, s)][pindex(l, m)][pindex(l, -n)] = (
                                FACTOR2(s, l, m, -n) * DD[0][l][pindex(smax, s)][pindex(l, m)][pindex(l, -n)]
                            )

                        # (s+, m-, n+)
                        if m > 0:
                            DD[0][l][pindex(smax, s)][pindex(l, -m)][pindex(l, n)] = (
                                FACTOR1(s, l, -m, n) * Dlmn(D, l, -m, n) * Dlmn(D, l, -m, s)
                            )
                            DD[1][l][pindex(smax, s)][pindex(l, -m)][pindex(l, n)] = (
                                FACTOR2(s, l, -m, n) * DD[0][l][pindex(smax, s)][pindex(l, -m)][pindex(l, n)]
                            )

                        # (s+, m-, n-)
                        if m > 0 and n > 0:
                            DD[0][l][pindex(smax, s)][pindex(l, -m)][pindex(l, -n)] = (
                                FACTOR1(s, l, -m, -n) * Dlmn(D, l, -m, -n) * Dlmn(D, l, -m, s)
                            )
                            DD[1][l][pindex(smax, s)][pindex(l, -m)][pindex(l, -n)] = (
                                FACTOR2(s, l, -m, -n) * DD[0][l][pindex(smax, s)][pindex(l, -m)][pindex(l, -n)]
                            )

                        # (s-, m+, n+)
                        if s > 0:
                            DD[0][l][pindex(smax, -s)][pindex(l, m)][pindex(l, n)] = FACTOR1(-s, l, m, n) * Dlmn(D, l, m, n) * Dlmn(D, l, m, -s)
                            DD[1][l][pindex(smax, -s)][pindex(l, m)][pindex(l, n)] = FACTOR2(-s, l, m, n) * DD[0][l][pindex(smax, -s)][pindex(l, m)][pindex(l, n)]

                        # (s-, m+, n-)
                        if s > 0 and n > 0:
                            DD[0][l][pindex(smax, -s)][pindex(l, m)][pindex(l, -n)] = FACTOR1(-s, l, m, -n) * Dlmn(D, l, m, -n) * Dlmn(D, l, m, -s)
                            DD[1][l][pindex(smax, -s)][pindex(l, m)][pindex(l, -n)] = FACTOR2(-s, l, m, -n) * DD[0][l][pindex(smax, -s)][pindex(l, m)][pindex(l, -n)]

                        # (s-, m-, n+)
                        if s > 0 and m > 0:
                            DD[0][l][pindex(smax, -s)][pindex(l, -m)][pindex(l, n)] = FACTOR1(-s, l, -m, n) * Dlmn(D, l, -m, n) * Dlmn(D, l, -m, -s)
                            DD[1][l][pindex(smax, -s)][pindex(l, -m)][pindex(l, n)] = FACTOR2(-s, l, -m, n) * DD[0][l][pindex(smax, -s)][pindex(l, -m)][pindex(l, n)]

                        # (s-, m-, n-)
                        if s > 0 and m > 0 and n > 0:
                            DD[0][l][pindex(smax, -s)][pindex(l, -m)][pindex(l, -n)] = FACTOR1(-s, l, -m, -n) * Dlmn(D, l, -m, -n) * Dlmn(D, l, -m, -s)
                            DD[1][l][pindex(smax, -s)][pindex(l, -m)][pindex(l, -n)] = FACTOR2(-s, l, -m, -n) * DD[0][l][pindex(smax, -s)][pindex(l, -m)][pindex(l, -n)]

                    else:
                        m_counter += 1                       

    return DD