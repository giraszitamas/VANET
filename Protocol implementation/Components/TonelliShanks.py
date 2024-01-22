def legendre(a, p):
    """
    Legendre symbol
    The Legendre symbol ( a | p) denotes the value of a ^ ((p-1)/2) (mod p)
    (a | p) ≡   1     if a is a square (mod p)
    (a | p) ≡  -1     if a is not a square (mod p)
    (a | p) ≡   0     if a ≡ 0
    :param a: integer;
    :param p: integer;
    :return: Legendre symbol;
    """
    return pow(a, (p - 1) // 2, p)


def tonelli(n, p):
    """
    A Tonelli-Shanks algoritmus implementációja, segítségével megoldható az r^2 ≡ n (mod p), ahol p prím; az n modulo p
     négyzetgyökét számolhatjuk ki vele.
    :param n: az alap;
    :param p: a modulus, prím;
    :return: r^2 ≡ n (mod p) --> r;
    """

    assert legendre(n, p) == 1, "not a square (mod p)"

    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r


if __name__ == '__main__':
    ttest = [(10, 13), (56, 101), (1030, 10009), (44402, 100049),
         (665820697, 1000000009), (881398088036, 1000000000039),
             (41660815127637347468140745042827704103445750172002, 10 ** 50 + 577)]
    for n, p in ttest:
        r = tonelli(n, p)
        assert (r * r - n) % p == 0
        print("n = %d p = %d" % (n, p))
        print("\t  roots : %d %d" % (r, p - r))
