def half_extended_gcd(aa, bb):
    """
    Kibővített euklideszi algoritmus;
    :param aa: egész szám;
    :param bb: egész szám;
    :return: LNKO(aa, bb);
    """
    lastrem, rem = abs(aa), abs(bb)
    x, lastx = 0, 1
    while rem:
        lastrem, (quotient, rem) = rem, divmod(lastrem, rem)
        x, lastx = lastx - quotient * x, x
    return lastrem, lastx


# Modular inverse: compute the multiplicative inverse i of a mod m:
#     i*a = a*i = 1 mod m
def modular_inverse(a, m):
    """
    Moduláris multiplikatív inverze számítása
    :param a: alap;
    :param m: modulus;
    :return: a^-1 (mod m);
    """
    g, x = half_extended_gcd(a, m)
    if g != 1:
        raise ValueError
    return x % m


def negate(P, q):
    """
    Projektív térbeli elem negáltja
    :param P: Projektív térbeli elem;
    :param q: prím, a véges test;
    :return: -P;
    """
    if P == 0: return P
    if P < 0: return -P
    if P < 0:
        P = P % q
    if q < P:
        P = P % q
    return (q - P) % q


def testBit(x, kth):
    """
    Adott helyi értékü bit ellenőrzés
    :return: igaz ha 1 és hamis 0;
    """
    return (x & 1 << kth) != 0
