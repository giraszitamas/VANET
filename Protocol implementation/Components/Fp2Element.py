class Fp2Element(object):

    def __init__(self, real, img):
        """
        Projektív térbeli elem inicializálása
        :param real: valós rész;
        :param img: képzetes rész;
        """
        self.real = real
        self.img = img

    def multFp2(self, e1, e2, q):
        """
        Projektív térbeli elemek szorzása
        :param e1: szorzandó;
        :param e2: szorzó;
        :param q: prím;
        :return: e1 * e2 (mod q);
        """
        r = ((e1.real * e2.real) - (e1.img * e2.img)) % q
        if r < 0: r += q
        i = ((e1.real * e2.img) + (e1.img * e2.real)) % q
        if i < 0: i += q
        return Fp2Element(r, i)

    def squareFp2(self, e, q):
        """
        Projektív térbeli elem négyzetre emelése
        :param e: az alap;
        :param q: prím;
        :return: e^2 (mod q);
        """
        r = ((e.real * e.real) - (e.img * e.img)) % q
        if r < 0: r += q
        i = (e.real * e.img * 2) % q
        if i < 0: i += q
        return Fp2Element(r, i)

    def conj(self, e, q):
        """
        Konjugált meghatározása
        :param e: alap;
        :param q: prím;
        :return: e konjugáltja;
        """
        return Fp2Element(e.real, q - e.img)

    def isOne(self):
        """
        :return: Ha a valós rész 1 és a képzetes 0 akkor igaz, egyébként hamis;
        """
        return self.real == 1 and self.img == 0

    def isZero(self):
        """
        :return: Ha a valós rész és a képzetes rész is 0 akkor igaz, egyébként hamis;
        """
        return self.real == 0 and self.img == 0

    def toString(self):
        """
        :return: az elem string reprezentációja;
        """
        return str(self.real) + " " + str(self.img) + "i"
