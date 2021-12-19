import argparse
from math import sqrt

from sage.all import (GF, ZZ, EllipticCurve, Integer, Mod, PolynomialRing,
                      hilbert_class_polynomial, kronecker, random)
from sage.arith.misc import fundamental_discriminant, is_prime


def choose_random_fundamental_discriminant(r):
    D = -Integer(Mod(int(random() * 1000), r))
    i = 0
    while not kronecker(D, r) == 1:  # expected number of iterations of the while loop is 2
        D = -Integer(Mod(int(random() * 1000), r))
        i += 1
    D = fundamental_discriminant(D)
    if not kronecker(D, r) == 1:
        return 0
    return D


def find_sqrt(p, D):
    F = GF(p)
    D = F(D)
    x = PolynomialRing(F, 'x').gen()
    f = x ^ 2 - D
    return ZZ(f.roots()[0][0])


def cornacchia(p, d):
    r"""
    input - a prime p, an absolute value of discriminant d
    output - a primitive solution (x,y) in integers to the equation
    x^2 + d*y^2 = p, if there exist one, otherwise return none.
    example:
    sage: find_sol(5,1) 
    sage: (2,1)
    """
    if kronecker(-d, p) == 1:
        t = find_sqrt(p, -d)
        bound = p ^ (.5)
        n = p
        while True:
            n, t = t, n % t
            if t < bound:
                break
        if ((p - t^2) / d) ^ (.5) in ZZ:
            return (t, ZZ(((p - t^2) / d) ^ (.5)))
        else:
            return None
    else:
        return None


def gen_curve(p, n, r, D=None):
    assert is_prime(p), "p has to be prime!"
    assert is_prime(n), "n has to be prime!"

    N = r * n  # number of points on the curve
    assert abs(N - (p + 1)) <= 2 * sqrt(p), "Hasse bound is not satisfied!"

    # b) set t = p + 1 - N
    t = p + 1 - N

    # c) choose a pair of integers (D, V) such that 4p - t^2 = DV^2
    DV2 = 4 * p - (t * t)

    # d) construct the Hilbert class polynomial P_D(X).
    pd_x = hilbert_class_polynomial(-D)

    # e) find a solution j0 in F(p) of P_D(X) = 0 modulo p.
    j0 = pd_x.any_root(GF(p))

    while True:
        # f) choose c ∈ F(p)* and construct an elliptic curve over F(p) with the j-invariant j0.
        c = GF(p).random_element()
        if c == 0:
            continue

        if j0 == 0:  # 2) j0 = 0_F
            E = EllipticCurve(GF(p), [0, c])
        elif j0 == 1728:  # 3) j0 = 1728
            E = EllipticCurve(GF(p), [c, 0])
        else:  # 1) j0 ≠ 0_F, 1 728
            a4 = (3 * c**2 * j0) / (1728 - j0)
            a6 = (2 * c**3 * j0) / (1728 - j0)
            E = EllipticCurve(GF(p), [a4, a6])

        # g) construct a random point G on E_{D,j0,c}[F(p)] such that G ≠ O_E and r·G ≠ O_E
        G = E.random_point()
        while G == 0 or r * G == 0:
            G = E.random_point()

        # h) set G = r·G
        G = r * G

        # i) If n·G = O_E, output curve parameters of E_{D,j0,c} and the base point G.
        #    If n·G ≠ O_E, go to step f) to choose another c.
        n = G.order()
        if is_prime(n):
            return E, G


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Pairing friendly generator.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', metavar='p', type=int,
                        help='the definition prime field order')
    parser.add_argument('-n', metavar='n', type=int,
                        help='the largest prime divisor of N = rn, which is the number of points')
    parser.add_argument('-r', metavar='r', type=int,
                        help='the cofactor of the curve')

    args = parser.parse_args()
    E, G = gen_curve(args.p, args.n, args.r)
    print(E, G)
