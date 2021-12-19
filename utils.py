from sage.all import EllipticCurve

def is_embedding_degree(E: EllipticCurve, k):
    return (E.base_field().order() ** k - 1) % E.order() == 0
