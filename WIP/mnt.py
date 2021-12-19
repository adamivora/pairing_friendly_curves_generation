import argparse
from math import ceil, log2

from sage.all import ZZ, IntegerModRing, floor, is_prime, is_square, sqrt
from sage.arith.misc import is_squarefree

"""
Implementation follows the appendix of the paper 
'On prime-order elliptic curves with embedding degrees k = 3, 4 and 6' by Koray Karabina and Edlyn Teske
"""

def pell_solve_1(D, m):
    """
    Algorithm 1 Pell Equation Solver
    Input: D ∈ ZZ, m ∈ ZZ \ {0} : D > m^2, D is not a perfect square
    Output: all minimal positive solutions (x, y) : x^2 - Dy^2 = m
    """
    assert m != 0
    assert D > m*m
    assert not is_square(D)

    B = [0, 1]
    P = [None, 0]
    Q = [None, 1]
    a = [None, floor(sqrt(D))]
    G = [1, a[1]]

    i = 1
    while True:
        i += 1
        P.append(a[i - 1] * Q[i - 1] - P[i - 1])
        Q.append((D - P[i]**2) / Q[i - 1])
        a.append(floor((P[i] + sqrt(D)) / Q[i]))
        B.append(a[i] * B[i - 1] + B[i - 2])
        G.append(a[i] * G[i - 1] + G[i - 2])

        if Q[i] == 1 and i % 2 == 1:
            break

    sols = []
    for j in range(1, i):
        f2 = (G[j]**2 - D * B[j]**2) / m
        if is_square(f2):
            f = sqrt(f2)
            sols.append((f * G[j], f * B[j]))
    return sols

def pell_solve_2(D, m):
    """
    Algorithm 2 Pell Equation Solver 2
    Input: D ∈ ZZ, m ∈ ZZ \ {0} : D ≤ m^2, D is not a perfect square
    Output: all fundamental solutions (x, y) : x^2 - Dy^2 = m
    """
    assert m != 0
    assert D <= m**2
    assert not is_square(D)

    u, v = pell_solve_1(D, 1)[0]
    if m > 0:
        L1 = 0
        L2 = sqrt(m * (u - 1) / (2 * D))
    else:
        L1 = sqrt(-m / D)
        L2 = sqrt(-m * (v + 1) / (2 * D))

    sols = []
    for y in range(ceil(L1), floor(L2) + 1):
        x2 = m + D * y**2
        if is_square(x2):
            x = sqrt(x2)

            if not (int(-x*x - y*y*D) % m == 0 and int(2*y*x) % m == 0):
                sols.append((x, y))
                sols.append((-x, y))
            else:
                sols.append((x, y))
    return sols

def gen_curve(N, z):
    """
    Algorithm 3 Elliptic curve parameters, embedding degree k = 6
    Input: N, z
    Output: EC parameters (q, n, D) where q - 1 is an N-bit prime, q^6 ≡ 1 (mod n)
    but q^i !≡ 1 (mod n) for 1 ≤ i ≤ 5, and D ≤ z (where 4q - t^2 = DY^2)
    """

    for D in range(9, 3 * z + 1, 24):
        if not is_squarefree(ZZ(D) / 3) or is_square(D):
            continue
        if not is_square(IntegerModRing(D)(-2)):
            continue

        if D > 64:
            sols = pell_solve_1(D, -8)
            if not sols:
                continue
            x0, y0 = sols[0]
        else:
            sols = pell_solve_2(D, -8)
            if not sols:
                continue
            x0, y0 = sols[0]

        u, v = pell_solve_1(D, 1)[0]
        x, y = x0, y0

        remainder = int(x) % 6
        if remainder == 5:
            remainder = -1
        if remainder == 1 or remainder == -1:
            while abs(x) <= 2**ceil(N / 2):
                l = (x - remainder) / 6
                if (((N - 3) / 2) <= log2(l) < ((N - 2) / 2)):
                    q = 4 * l**2 + 1
                    n = 4 * l**2 - remainder * 2 * l + 1
                    if is_prime(q) and is_prime(n):
                        return (q, n, D / 3)
                
                old_x = x
                x = x * u + y * v * D
                y = old_x * v + u * y

        x, y = x0 * u - y0 * v * D, u * y0 - x0 * v
        remainder = int(x) % 6
        if remainder == 5:
            remainder = -1
        if remainder == 1 or remainder == -1:
            while abs(x) <= 2**ceil(N / 2):
                l = (x - remainder) / 6
                if l < 0:
                    break
                if (((N - 3) / 2) <= log2(l) < ((N - 2) / 2)):
                    q = 4 * l**2 + 1
                    n = 4 * l**2 - remainder * 2 * l + 1
                    if is_prime(q) and is_prime(n):
                        return (q, n, D / 3)
                
                old_x = x
                x = x * u - y * v * D
                y = u * y - old_x * v

    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Pairing friendly generator.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--N', metavar='p', type=int,
                        help='the q bitlength')
    parser.add_argument('--z', metavar='p', type=int,
                        help='the maximum value of the CM discriminant')

    args = parser.parse_args()
    q, n, D = gen_curve(args.N, args.z)
    print(q, n, D)
