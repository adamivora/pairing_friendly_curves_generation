import argparse

from sage.all import (GF, EllipticCurve, Integer, Mod,
                      fundamental_discriminant, hilbert_class_polynomial,
                      kronecker, power_mod, randint, random, random_prime,
                      sqrt)
from sage.arith.misc import is_prime

from utils import is_embedding_degree


def find_prime(num_bits, B):
    p = 0
    while not p % B == 1:
        p = random_prime(2 ** num_bits, lbound=2 ** (num_bits - 1))
    return p

def find_fundamental_discriminant(n):
    D = 0
    while not kronecker(D, n) == 1:
        D = -Integer(Mod(int(random() * 1000), n))
    D = fundamental_discriminant(D)
    assert kronecker(D, n) == 1
    return D

def find_element_of_order(B, n):
    assert n % B == 1
    h = 0
    def order(h, k, p):
        bool = True
        g = h
        for _ in range(1, k):
            bool = bool and g != 1
            g = Mod(g * h, p)
        bool = bool and g == 1
        return bool
    while not order(h, B, n):
        h = power_mod(randint(2, n - 1), (n - 1) // B, n)
    return h

def gen_curve(B, num_bits):
    p = 0
    # g) If p is not prime, then go to step a).
    while not is_prime(p):
        # a) Choose a small square-free positive integer D and n in R such that -D is a square modulo n.
        n = find_prime(num_bits, B)
        D = -find_fundamental_discriminant(n)

        # b) Find a B-th primitive root of unity z in F(n).
        gf = GF(n)
        z = gf(find_element_of_order(B, n))
        assert z.multiplicative_order() == B

        # c) t' = z + 1.
        t_prime = z + 1

        # d) y' = (t' - 2) / √ (-D) (mod n).
        y_prime = (t_prime - 2) // sqrt(gf(-D))

        # e) Let t be an integer such that t is equal to t' modulo n,
        # and let y be an integer such that y is equal to y' modulo n.
        t, y = t_prime, y_prime

        # f) p = (t^2 + Dy^2) / 4.
        p = ((t ** 2) + D * (y ** 2)) / 4

    # i) Construct the Hilbert class polynomial P_D(X).
    pd_x = hilbert_class_polynomial(D)

    # j) Find a solution j_0 in F(p) of P_D(X) = 0 modulo p.
    j0 = pd_x.any_root(GF(p))

    while True:
        # k) Choose c ∈ F(p)* and construct an elliptic curve over F(p) with the j-invariant j_0.
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
        
        # l) Set a cofactor r = (p + 1 - t) / n.
        r = (p + 1 - t) / n

        # m) construct a random point G on E_{D,j_0,c}[F(p)] such that G ≠ O_E and r·G ≠ O_E
        G = E.random_point()
        while G == 0 or r * G == 0:
            G = E.random_point()

        # n) set G = r·G
        G = r * G

        # o) If n·G = O_E, output n, G, and the elliptic curve E.
        n = G.order()
        if is_prime(n):
            break

        # p) Else, go to step k) to choose another c ∈ F(p)*.

    assert is_embedding_degree(E, B), "Embedding degree is not B!"
    return n, G, E

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Pairing friendly generator.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--num-curves', metavar='m', type=int,
                        help='how many curves to generate', default=10000)
    parser.add_argument('--B', metavar='m', type=int,
                        help='the desired embedding degree', default=12)
    parser.add_argument('--num-bits', metavar='p', type=int,
                        help='number of bits of the prime that specifies the underlying finite field', default=128)

    args = parser.parse_args()
    curve = gen_curve(args.B, args.num_bits)
    if curve is None:
        print('No curve found with given parameters!')
        exit(1)

    n, G, E = curve
    print(E)
