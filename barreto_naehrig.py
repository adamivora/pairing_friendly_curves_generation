import argparse
import itertools

from sage.all import GF, EllipticCurve, kronecker, sqrt
from sage.arith.misc import is_prime

from utils import is_embedding_degree


def gen_curve(m, p_max):
    # a) Let P(u) = 36u^4 + 36u^3 + 24u^2 + 6u + 1
    def P(u): return 36 * u**4 + 36 * u**3 + 24 * u**2 + 6 * u + 1

    # b) Compute the smallest u ≈ 2^(m/4) such that ⎾log_2 P(−u)⏋ = m.
    u = int(-(2**(m / 4)) / sqrt(6))

    p = P(-u)
    found_curve = False

    # c) While p ≤ p_max
    while True:
        while p.bit_length() <= p_max:
            # 1) t = 6u^2 + 1.
            t = 6 * (u*u) + 1

            # 2) p = P(−u) and n = p + 1 − t.
            p = P(-u)
            n = p + 1 - t

            # 3) If p and n are prime, then go to step e).
            if is_prime(p) and is_prime(n):
                found_curve = True
                break

            # 4) p = P(u) and n = p + 1 − t.
            p = P(u)
            n = p + 1 - t

            # 5) If p and n are prime, then go to step e).
            if is_prime(p) and is_prime(n):
                found_curve = True
                break

            # 6) u = u + 1 and go to step 1).
            u += 1

        # d) Stop and output "fail".
        if not found_curve:
            return None

        gf = GF(p)
        # f) b = 0.
        for b in itertools.count(1):
            # g) If b + 1 is not represented by b + 1 = y0^2 modulo p for an integer y_0, then b = b + 1 and go to step g).
            if kronecker(b + 1, p) != 1:
                continue

            # h) Set an elliptic curve E: y^2 = x^3 + b.
            E = EllipticCurve(gf, [0, b])

            # i) Compute a square root y_0 = √(b + 1) modulo p.
            y0 = gf(b + 1).sqrt()

            # j) Set the basepoint G = (1, y_0) ∈ E.
            G = E(1, y0)

            # k) If n⋅G ≠ O_E, then set b = b + 1 and go to step g).
            if n * G == 0:
                break

        assert is_embedding_degree(E, 12), "Embedding degree is not 12!"
        u += 1
        # l) Output p, E, n, and G.
        yield p, E, n, G


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Pairing friendly generator.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--num-curves', metavar='m', type=int,
                        help='how many curves to generate', default=10000)
    parser.add_argument('--m', metavar='m', type=int,
                        help='desired curve order in bits', default=100)
    parser.add_argument('--p-max', metavar='p', type=int,
                        help='the maximum p bitlength', default=1024)

    args = parser.parse_args()

    curves = gen_curve(args.m, args.p_max)
    remaining_curves = args.num_curves
    for curve in curves:
        if curve is None:
            print('No curve found with given parameters!')
            exit(1)
        p, E, n, G = curve
        print(E)
        remaining_curves -= 1
        if remaining_curves == 0:
            break
