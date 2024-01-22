import hashlib
import random
from Components.TatePairing import TatePairing


class BonehFranklin(object):

    def H2(self, bp, n):
        """
        MIN_Length = 256
        SHA-256 hash, melynek visszatérési értékének hosszát a második paraméter adja meg
        :param bp: hashelni kívánt pont;
        :param n: hossz;
        :return: adott hosszúságú hash;
        """
        bp = str(bp)
        BPhash = hashlib.sha256()
        BPhash.update(bp.encode())
        BPhash = int(BPhash.hexdigest(), base=16)
        BPhash = str(BPhash)
        while len(BPhash) < n:
            BPhash = "a" + BPhash
        return BPhash

    def encryp(self, msg, Qr, gammaP, P, EC):
        """
        Titkosítás
        :param msg: titkosítandó üzenet;
        :param Qr: küldő nyilvános kulcs;
        :param gammaP: rendszer nyilvános paramétere;
        :param P: rendszer nyilvános paramétere;
        :param EC: elliptikus görbe osztály példánya;
        :return: titkosíott üzenet egy értékpár formájában (rP, v);
        """
        r = random.getrandbits(128)
        rP = EC.mulJ(P, r)
        rgammaP = EC.mulJ(gammaP, r)
        gIDr = TatePairing.computeF(TatePairing, Qr, rgammaP, EC)
        gIDr = self.H2(self, gIDr.toString(), len(msg))
        v = sxor(msg, gIDr)
        return rP, v

    def decrypt(self, rP, v, gammaQr, EC):
        """
        Visszafejtés
        :param rP: titkosított üzenet megfelelő része (rP);
        :param v: titkosított üzenet megfelelő része (v);
        :param gammaQr: címzett titkos kulcsa;
        :param EC: elliptikus görbe osztály paramétere;
        :return: visszafejtett üzenet;
        """
        h2 = self.H2(self, TatePairing.computeF(TatePairing, gammaQr, rP, EC).toString(), len(v))
        return sxor(v, h2)


def sxor(s1, s2):
    """
    Stringek közötti kizáró vagy művelet. A stringekből először karaktertömböket hoz létre, majd a karakterek ASCII
    értékeit XOR-olja és visszaalakítja azokat karakterré ezután a kapott értékeket összefűzi stringgé.
    :param s1: 1# string;
    :param s2: 2# string;
    :return: 1# string XOR 2# string;
    """
    # convert strings to a list of character pair tuples
    # go through each tuple, converting them to ASCII code (ord)
    # perform exclusive or on the ASCII code
    # then convert the result back to ASCII (chr)
    # merge the resulting array of characters as a string
    return ''.join(chr(ord(a) ^ ord(b)) for a, b in zip(s1, s2))
