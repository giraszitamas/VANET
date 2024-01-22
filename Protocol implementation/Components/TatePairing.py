from Components.Fp2Element import Fp2Element
from Components import Complex
from Components.HalfComplex import HalfComplex
from Components import Math


class TatePairing(object):

    def computeF(self, P, Q, ec):
        """
        This class implements the compressed Tate pairing as described in Michael Scott, Paulo S. L. M. Barreto:
        Compressed Pairings. CRYPTO 2004: 140-156.
        The output value is reduced to half length. Instead of keeping the full a+bi value of the Tate pairing, it may
        be possible for cryptographic purposes to discard b altogether, leaving the values defined only up to
        conjugation, which means one of the pairing arguments will only be defined up to a sign. If the output of the
        Tate pairing e(P,Q) = a+bi, then the output of the compressed Tata pairing ce(P,Q) is only a. Note that when
        using compressed pairing ce(P,Q) =ce(-P,Q)=ce(-P,-Q)=ce(P,-Q), so only use the compressed Tate Pairing when the
        sign is not important.
        @author Changyu Dong
        @version 1.0
        @see uk.ac.ic.doc.jpair.pairing.HalfComplex
        @see uk.ac.ic.doc.jpair.pairing.HalfComplexField
        :param P: elliptikus görbén lévő pont;
        :param Q: elliptikus görbén lévő pont;
        :param ec: elliptikus görbe példánya;
        :return: P és Q pontok bilineáris leképezése;
        """

        f = Fp2Element(1, 0)
        V = ec.aToJ(P)
        nP = ec.neg(P)
        b = naf(ec.r, 2)
        for i in range(len(b) - 2, -1, -1):
            u = self.encDouble(self, V, Q, ec)
            f = Fp2Element.multFp2(Fp2Element, Fp2Element.squareFp2(Fp2Element, f, ec.q), u, ec.q)
            if b[i] == 1:
                u = self.encAdd(self, V, P, Q, ec)
                f = Fp2Element.multFp2(Fp2Element, f, u, ec.q)
            if b[i] == -1:
                u = self.encAdd(self, V, nP, Q, ec)
                f = Fp2Element.multFp2(Fp2Element, f, u, ec.q)

        finalExp = (ec.q + 1) // ec.r
        conj = Complex.conjugate(f, ec.q)
        f = Complex.divide(conj, f, ec.q)
        f = Complex.toHalfComplex(f, ec.q)
        return HalfComplex.HCpow(HalfComplex, f, finalExp, ec.q)

    def compute(self, P, Q, ec):
        """
        A Tate pairing implementation.
        This implementation uses the pairing friendly curve {@code Y^2 = X^3 + aX + b} defined over GF(p)
        where {@code p = 3 mod 4}. G1 is taken as an order-q subgroup of the group formed by all the points on the curve.
        The curve has an embedding degree of 2. It has a corresponding twisted curve {@code Y^2 = X^3 + aX - b}.
        Points from the twisted curve are used in the computation as elements in G2 to avoid operations in the extension field.
        The algorithm is taken from "Efficient Computation of Tate Pairings in Projective Coordinates over General Characteristic Fields",
        Proc. 7th Int. Conference on Inf Security and Cryptology (ICISC 2004), Eds. C.Park and S. Chee,
        LNCS 3506, Springer 2005,  168-181.
        @author Changyu
        """
        f = Fp2Element(1, 0)
        V = ec.aToJ(P)
        nP = ec.neg(P)
        n = ec.r - 1
        b = naf(ec.r, 2)
        for i in range(len(b) - 2, -1, -1):
            u = self.encDouble(self, V, Q, ec)
            # f = f^2 * u % prime
            f = Fp2Element.multFp2(Fp2Element, Fp2Element.squareFp2(Fp2Element, f, ec.q), u, ec.q)
            if b[i] == 1:
                u = self.encAdd(self, V, P, Q, ec)
                f = Fp2Element.multFp2(Fp2Element, f, u, ec.q)
            if b[i] == -1:
                u = self.encAdd(self, V, nP, Q, ec)
                f = Fp2Element.multFp2(Fp2Element, f, u, ec.q)

        finalExp = (ec.q + 1) // ec.r
        conj = Complex.conjugate(f, ec.q)
        f = Complex.divide(conj, f, ec.q)
        return Complex.complexPow(f, finalExp, ec)

    def encDouble(self, P, Q, EC):
        """
        used by tate pairing, point doubling in Jacobian coordinates, and return the value of f
        :param P: Jacobian térbeli elliptikus görbén lévő pont duplázáshoz (JacobianPoint);
        :param Q: affin térbeli elliptikus görbén lévő pont, viszonyitási pont (ECPoint);
        :param EC: elliptikus görbe példánya;
        :return: Fp2Element(real, img);
        """
        x = P.x
        y = P.y
        z = P.z
        q = EC.q

        t1 = (y * y) % q

        t2 = (x * t1) % q
        t2 = (t2 + t2) % q
        t2 = (t2 + t2) % q
        #         t2 = (4 * t2) % q

        t3 = (t1 ** 2) % q
        t3 = (t3 + t3) % q
        t3 = (t3 + t3) % q
        t3 = (t3 + t3) % q
        #         t3 = (6 * t3) % q

        t4 = (z ** 2) % q

        t5 = (x ** 2) % q
        t5 = (3 * t5) % q
        t5 = (t5 + EC.a * (t4 ** 2)) % q

        x3 = (t5 ** 2) % q
        x3 = (x3 - (t2 + t2)) % q

        y3 = (t5 * (t2 - x3)) % q
        y3 = (y3 - t3) % q

        z3 = (y * z) % q
        z3 = (z3 + z3) % q

        P.x = x3
        P.y = y3
        P.z = z3

        real = (t4 * Q.x) % q
        real = (real + x) % q
        real = (t5 * real) % q
        real = (real - t1) % q
        real = (real - t1) % q

        img = (z3 * t4) % q
        img = (img * Q.y) % q

        return Fp2Element(real, img)

    def encAdd(self, A, P, Q, EC):
        """
        used by Tate paring, add two point, save result in the first argument, return the value of f
        :param A: Jacobian térbeli elliptikus görbepont, az összeadás egyik tagja
        :param P: Jacobian térbeli elliptikus görbepont, az összeadás másik tagja
        :param Q: Jacobian térbeli elliptikus görbepont, viszonyitási pont
        :param EC: elliptikus görbe példánya;
        :return: Fp2Element(real, img)
        """

        x1 = A.x
        y1 = A.y
        z1 = A.z

        q = EC.q

        x = P.x
        y = P.y

        #          t1=z1^2
        t1 = (z1 ** 2) % q
        #          t2=z1t1
        t2 = (z1 * t1) % q
        #          t3=xt1
        t3 = (x * t1) % q
        #          t4=Yt2
        t4 = (y * t2) % q
        #          t5=t3-x1
        t5 = (t3 - x1) % q
        #          t6=t4-y1
        t6 = (t4 - y1) % q
        #          t7=t5^2
        t7 = (t5 ** 2) % q
        #          t8=t5t7
        t8 = (t5 * t7) % q
        #          t9=x1t7
        t9 = (x1 * t7) % q

        x3 = (t6 ** 2) % q
        x3 = (x3 - (t8 + t9 + t9)) % q
        # x3 = this.field.subtract(x3, this.field.add(t8, this.field.add(t9, t9)));

        y3 = (t6 * (t9 - x3)) % q
        y3 = (y3 - y1 * t8) % q

        z3 = (z1 * t5) % q

        A.x = x3
        A.y = y3
        A.z = z3

        img = (z3 * Q.y) % q
        real = (Q.x + x) % q
        real = (real * t6) % q
        real = (real - (z3 * y)) % q

        return Fp2Element(real, img)


def naf(k, w):
    """
    k Window-Non-Adjacent - windowed naf form of BigInt k, w is the window size
    The window NAF is at most 1 element longer than the binary
    representation of the integer k. byte can be used instead of short or
    int unless the window width is larger than 8. For larger width use
    short or int. However, a width of more than 8 is not efficient for
    m = log2(q) smaller than 2305 Bits. Note: Values for m larger than
    1000 Bits are currently not used in practice.
    :param k: big integer k;
    :param w: ablak mérete;
    :return: k Window-Non-Adjacent alakja;
    """
    wnaf = [None] * (k.bit_length() + 1)
    pow2wB = 1 << w
    i = 0
    length = 0
    while k >= 1:
        # if k is odd
        if Math.testBit(k, 0):
            reminder = k % pow2wB
            # if reminder > 2^(width - 1) - 1
            if Math.testBit(reminder, w - 1):
                wnaf[i] = reminder - pow2wB
            else:
                wnaf[i] = reminder
            k = k - wnaf[i]
            length = i
        else:
            wnaf[i] = 0
        k = k >> 1
        i += 1
    length += 1
    wnafShort = wnaf[0:length]
    return wnafShort
