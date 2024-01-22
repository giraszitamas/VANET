import hashlib
import random
import time

from Components.ECPoint import ECPoint
from Components.Fp2Point import Fp2Point
from Components.JacobPoint import JacobPoint
from Components import Math
from Components import TonelliShanks


class EC(object):
    """System of Elliptic Curve"""

    def __init__(self):
        """
        - a, b: az elliptikus görbe paraméterei;
        - q: prím szám;
        formula y^2=x^3+x;
        q =  7313295762564678553220399414112155363840682896273128302543102778210584118101444624864132462285921835023839111762785054210425140241018649354445745491039387;
        rend =  730750818665451459101842416358141509827966402561;
        cofactor 10007920040268628970387373215664582404186858178692152430205359413268619141100079249246263148037326528074908;
        embadded degree k = 2;
        """
        self.a = 1
        self.b = 0
        self.q = 7313295762564678553220399414112155363840682896273128302543102778210584118101444624864132462285921835023839111762785054210425140241018649354445745491039387
        # just as unique ZERO value representation for "add": (not on curve)
        self.zero = ECPoint(self, 0, 0)
        self.r = 730750818665451459101842416358141509827966402561
        self.basepoint = ECPoint(self,
                                 709493439226547518671621749763636350732076745439920897858029934151951208424241756347352095017515948632575152914173127042095409296000312028409473049395679,
                                 2539945134587823181281630459392477399176359428987561216702875801069327651488451605997719345033752243529105120554302078221704343116029570188151122387224856)
        self.cofactor = 10007920040268628970387373215664582404186858178692152430205359413268619141100079249246263148037326528074908
        pass

    def neg(self, p):
        """Pont negáltja
        :param p: elliptikus görbepont;
        """
        return ECPoint(self, p.x, -p.y % self.q)

    def at(self, RID):
        """Az elliptikus görbe egy pontja az x helyen
        :param RID: string, valós azonosító, ebből lesz a x < q, vagy None, akkor véletlen pont generálása;
        :returns: ((x, y), (x,-y)) vagy "Not found exception";
        """
        while True:
            x, ysq = ecCalc(RID, self.q)
            while TonelliShanks.legendre(ysq, self.q) != 1:
                x, ysq = ecCalc(RID, self.q)
            y = TonelliShanks.tonelli(ysq, self.q)

            P = ECPoint(self, x, y)
            P = self.mulJ(P, self.cofactor)
            tmp = self.mulJ(P, self.r)

            if tmp.x is None and tmp.y is None:
                return P

    def at_gen(self, RID):
        """Az elliptikus görbe egy pontja az x helyen
        :param RID: string, valós azonosító, ebből lesz a x < q, vagy None, akkor véletlen pont generálása;
        :returns: ((x, y), (x,-y)) vagy "Not found exception";
        """
        if RID is not None:
            x = hashStringToInt(RID)
            while x > self.q:
                x = hashStringToInt(RID)
            return self.mulJ(self.basepoint, x)
        else:
            x = random.getrandbits(128)
            while x > self.q:
                x = hashStringToInt(RID)
            return self.mulJ(self.basepoint, x)

    # Return the (x,y) point where this line intersects our curve
    #  Q1 and Q2 are two points on the line of slope m
    def line_intersect(self, Q1, Q2, m):
        """
        Meghatározza egy adott egyenes és a görbe metszéspontját
        :param Q1: egy az egyenesen lévő pont;
        :param Q2: egy az egyenesen lévő pont;
        :param m: slope;
        :return: (x,y) metszéspont;
        """
        v = (Q1.y + self.q - (m * Q1.x) % self.q) % self.q
        x = (m * m + self.q - Q1.x + self.q - Q2.x) % self.q
        y = (self.q - (m * x) % self.q + self.q - v) % self.q
        return ECPoint(self, x, y)

    # Prime field division: return num/den mod p
    def field_div(self, num, den):
        """Prím test osztás: return num/den mod p"""
        inverse_den = Math.modular_inverse(den % self.q, self.q)
        return self.field_mul(num % self.q, inverse_den)

    # Prime field multiplication: return a*b mod p
    def field_mul(self, a, b):
        """Prím test szorzás: return a*b mod p"""
        return (a * b) % self.q

    # Prime field exponentiation: raise num to power mod p
    def field_exp(self, num, power):
        """Prím test hatványozás: raise num to power mod p"""
        return pow(num % self.q, power, self.q)

    def add2(self, Q1, Q2):
        """
        Az elliptikus görbe pontjainak összeadása
        :param Q1: az elliptikus görbe egy pontja;
        :param Q2: az elliptikus görbe egy pontja;
        :return: Q1 + Q2;
        """
        # Identity special cases
        if Q1.x == self.q:  # Q1 is identity
            return Q2
        if Q2.x == self.q:  # Q2 is identity
            return Q1

        # Equality special cases
        if Q1.x == Q2.x:
            if Q1.y == Q2.y:  # adding point to itself
                return self.double(Q1)
            else:  # vertical pair--result is the identity
                return Q1
        # Ordinary case
        m = self.field_div(Q1.y + self.q - Q2.y, Q1.x + self.q - Q2.x)
        return self.line_intersect(Q1, Q2, m)

    def double(self, P):
        """
        Az elliptikus görbe egy pontjának duplázása
        :param P: elliptikus görbepont;
        :return: P + P;
        """
        s = (3 * P.x ** 2 + self.a * Math.modular_inverse(2 * P.y, self.q) % self.q)
        if s <= 0:
            s += self.q
        t = (P.y - s * P.x) % self.q
        if t <= 0:
            t += self.q
        r1 = (s ** 2 - 2 * P.x) % self.q
        if r1 <= 0:
            r1 += self.q
        r2 = (s * r1 + -t) % self.q
        if r2 <= 0:
            r2 += self.q
        return ECPoint(self, r1, r2)

    def mulJ(self, P, x):
        """
        Affin térbeli pontok szorzása
        :param P: affin térbeli elliptikus görbén lévő pont;
        :param x: konstans érték
        :return: Jacobian térbeli pont, melynek értéke Q1 + Q2;
        """
        return self.jToA(self.jMultiply(P, x))

    def aToJ(self, P):
        """
        Affin térbeli pont átalakítása Jacobian térbeli ponttá
        :param P: affin térbeli elliptikus görbepont P (x, y);
        :return: Jacobian térbeli pont P(x, y, 1);
        """
        return JacobPoint(P.x, P.y, 1)

    def jToA(self, P):
        """
        Jacobian térbeli pont átalakítása affin térbeli ponttá
        :param P: Jacobian térbeli pont (x, y, z);
        :return: Affint térbeli pont P (x, y);
        """
        if P.isInfinity():
            return Fp2Point(None, None)
        zInverse = Math.modular_inverse(P.z, self.q)
        square = (zInverse ** 2) % self.q
        # x =X/Z^2
        x = (P.x * square) % self.q
        # y=Y/Z^3
        y = (P.y * (square * zInverse)) % self.q
        return Fp2Point(x, y)

    def jMultiply(self, P, x):
        """
        Affin térbeli pont szorzása
        :param P: affin térbeli elliptikus görbén lévő pont;
        :param x: konstans
        :return: Jacobian térbeli pont P (x, y, z) = Q1 * Q2;
        """
        result = self.aToJ(P)
        degree = x.bit_length() - 2
        for i in range(degree, -1, -1):
            self.jDbl(result)
            if Math.testBit(x, i):
                self.jAdd(result, P)
        return result

    def jDbl(self, P):
        """
        Jacobian térbeli pont duplázása
        :param P: Jacobian térbeli pont;
        :return: P + P;
        """
        x = P.x
        y = P.y
        z = P.z
        q = self.q

        #         t1 = y ** 2
        t1 = (y * y) % q
        #         t2 = 4x * t1
        t2 = (x * t1) % q
        t2 = (t2 + t2) % q
        t2 = (t2 + t2) % q
        #         t3 = (8 * t1) * 2
        t3 = (t1 ** 2) % q
        t3 = (t3 + t3) % q
        t3 = (t3 + t3) % q
        t3 = (t3 + t3) % q
        #         t4 = z ** 2
        t4 = (z ** 2) % q
        #         t5 = 3x ** 2 + (a * t4) ** 2
        t5 = (x ** 2) % q
        t5 = (3 * t5) % q
        t5 = (t5 + self.a * (t4 ** 2)) % q
        #         x3 = t5 ** 2 - 2 * t2
        x3 = (t5 ** 2) % q
        x3 = (x3 - (t2 + t2)) % q
        #         y3 = t5 * (t2 - x3) - t3
        y3 = (t5 * (t2 - x3)) % q
        y3 = (y3 - t3) % q
        #         z3 = 2 * y * z
        z3 = (y * z) % q
        z3 = (z3 + z3) % q

        P.x = x3
        P.y = y3
        P.z = z3
        return

    def jAdd(self, P, Q):
        """
        Jacobian térbeli pontok összeadása
        :param P: Jacobian térbeli pont;
        :param Q: Jacobian térbeli pont;
        :return: P + Q;
        """

        x1 = P.x
        y1 = P.y
        z1 = P.z

        q = self.q

        x = Q.x
        y = Q.y

        # t1=z1^2
        t1 = (z1 ** 2) % q
        # t2=z1t1
        t2 = (z1 * t1) % q
        # t3=xt1
        t3 = (x * t1) % q
        # t4=Yt2
        t4 = (y * t2) % q
        # t5=t3-x1
        t5 = (t3 - x1) % q
        # t6=t4-y1
        t6 = (t4 - y1) % q
        # t7=t5^2
        t7 = (t5 ** 2) % q
        # t8=t5t7
        t8 = (t5 * t7) % q
        # t9=x1t7
        t9 = (x1 * t7) % q
        # x3=t6^2-(t8+2t9)
        x3 = (t6 ** 2) % q
        x3 = (x3 - (t8 + t9 + t9)) % q
        # y3=t6(t9-x3)-y1t8
        y3 = (t6 * (t9 - x3)) % q
        y3 = (y3 - y1 * t8) % q
        # z3=z1t5
        z3 = (z1 * t5) % q

        P.x = x3
        P.y = y3
        P.z = z3

        return

    pass


def hashStringToInt(x_str):
    """
    Stringből állít elő egész számot hash függvény segítségével
    :param x_str: string;
    :return: integer;
    """
    RID = x_str + str(time.time())
    hash_def = hashlib.sha256()
    hash_def.update(RID.encode())
    return int(hash_def.hexdigest(), base=16)


def hashTimeToInt(x_str):
    """
    String időbélyegből állít elő egész számot hash függvény segítségével
    :param x_str: string;
    :return: integer;
    """
    hash_def = hashlib.sha256()
    hash_def.update(x_str.encode())
    return int(hash_def.hexdigest(), base=16)


def hashMsgAndTime(x_str):
    """
    Stringből állít elő egész számot hash függvény segítségével
    :param x_str: string;
    :return: integer és időbélyeg;
    """
    t = time.time()
    RID = x_str + str(t)
    hash_def = hashlib.sha256()
    hash_def.update(RID.encode())
    return int(hash_def.hexdigest(), base=16), t


def hashMsgAndTimeStr(x_str, timeStamp):
    """
    Stringből és időbélyegből állít elő egész számot hash függvény segítségével
    :param x_str: string;
    :return: integer;
    """
    RID = x_str + str(timeStamp)
    hash_def = hashlib.sha256()
    hash_def.update(RID.encode())
    return int(hash_def.hexdigest(), base=16)


def ecCalc(RID, q):
    """
    Véletlen pont megtalálásához szükséges megfelelő érték keresése, x < q
    :param RID: valós azonosító;
    :param q: prím, a véges test;
    :return: megfelelő x véletlen;
    """
    if RID is not None:
        x = hashStringToInt(RID)
        while x > q:
            x = hashStringToInt(RID)
    else:
        x = random.getrandbits(128)
        while x > q:
            x = random.getrandbits(128)
    assert x < q
    ysq = (x ** 3 + x) % q
    return [x, ysq]
