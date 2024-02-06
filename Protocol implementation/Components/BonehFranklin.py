import hashlib
import random
from Components.ECPoint import ECPoint
from Components.TatePairing import TatePairing
import math
def point_to_byte(P):
    x_byte = P.x.to_bytes(64)
    y_byte = P.y.to_bytes(64)
    return x_byte + y_byte

def byte_to_point(B, ec):
    x_byte = int.from_bytes(B[0:64])
    y_byte = int.from_bytes(B[64:])

    return ECPoint(ec, x_byte, y_byte)


class BonehFranklin_stream(object):

    def H2(self, bp):
        bp = bp.to_bytes(64)
        BPhash = hashlib.sha256()
        BPhash.update(bp)
        BPhash = int(BPhash.hexdigest(), base=16).to_bytes(32)

        return BPhash

    def encryp(self, msg, Qr, gammaP, P, EC):
        r = random.getrandbits(128)
        rP = EC.mulJ(P, r)
        rgammaP = EC.mulJ(gammaP, r)
        gIDr = TatePairing.computeF(TatePairing, Qr, rgammaP, EC)
        gIDr = self.H2(self, gIDr.real)
        v = self.xor(self,msg, gIDr)
        return rP, v

    def decrypt(self, rP, v, gammaQr, EC):
        h2 = self.H2(self, TatePairing.computeF(TatePairing, gammaQr, rP, EC).real)
        return self.xor(self,v, h2)


    def xor(self,s1, s2):
        ret = (s1[0] ^ s2[0]).to_bytes()
        for i in range(1,len(s1)):
            ret+=(s1[i] ^ s2[i%len(s2)]).to_bytes()
        return ret


class BonehFranklin_block(object):

    def H2(self, bp):
        bp = bp.to_bytes(64)
        BPhash = hashlib.sha256()
        BPhash.update(bp)
        BPhash = int(BPhash.hexdigest(), base=16).to_bytes(32)

        return BPhash

    def encryp(self, msg, Qr, gammaP, P, EC):


        r = random.getrandbits(128)
        rP = EC.mulJ(P, r)
        rgammaP = EC.mulJ(gammaP, r)
        gIDr = TatePairing.computeF(TatePairing, Qr, rgammaP, EC)
        gIDr = self.H2(self, gIDr.real)

        to_enc_msg = msg[:32]


        ret = point_to_byte(rP) + self.xor(self, to_enc_msg, gIDr)
        for i in range(1, int(len(msg)/32)):
            r = random.getrandbits(128)
            rP = EC.mulJ(P, r)
            rgammaP = EC.mulJ(gammaP, r)
            gIDr = TatePairing.computeF(TatePairing, Qr, rgammaP, EC)
            gIDr = self.H2(self, gIDr.real)

            to_enc_msg = msg[i*32:(i+1)*32]


            v = self.xor(self, to_enc_msg, gIDr)

            ret += point_to_byte(rP) + v

        return ret

    def decrypt(self, M1_block, gammaQr, EC):

        c = M1_block[:160]
        rP = byte_to_point(c[:128],EC)
        v = c[128:]


        h2 = self.H2(self, TatePairing.computeF(TatePairing, gammaQr, rP, EC).real)


        ret = self.xor(self,v, h2)

        for i in range(1,int(len(M1_block)/160)):
            c = M1_block[i*160:(i+1)*160]
            rP = byte_to_point(c[:128],EC)
            v = c[128:]

            h2 = self.H2(self, TatePairing.computeF(TatePairing, gammaQr, rP, EC).real)
            ret += self.xor(self,v, h2)
        return ret


    def xor(self,s1, s2):
        ret = (s1[0] ^ s2[0]).to_bytes()
        for i in range(1,len(s1)):
            ret+=(s1[i] ^ s2[i]).to_bytes()
        return ret

