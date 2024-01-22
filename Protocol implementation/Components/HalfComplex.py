from Components import Math


class HalfComplex(object):

    def __init__(self, field, real):
        """
        :param field: véges test;
        :param real: valós rész;
        """
        self.field = field
        self.real = real

    def HCpow(self, p, n, q):
        """
        Computation of Lucas sequence elements
        :param p: alap;
        :param n: exponens;
        :param q: prím, a véges test;
        :return: p^n (mod q);
        """
        p = p.real << 1

        v0 = 2
        v1 = p

        t = n.bit_length() - 1

        for j in range(t, -1, -1):
            if Math.testBit(n, j):
                v0 = (v0 * v1) % q - p
                v1 = (v1 ** 2) % q - 2
            else:
                v1 = (v0 * v1) % q - p
                v0 = (v0 ** 2) % q - 2

        if Math.testBit(v0, 0):
            return HalfComplex(q, (v0 * Math.modular_inverse(2, q)) % q)
        else:
            return HalfComplex(q, v0 >> 1)

    def toString(self):
        """
        :return: az elem string reprezentációja;
        """
        return str(self.real)
