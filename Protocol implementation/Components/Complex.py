from Components.Fp2Element import Fp2Element
from Components import Math
from builtins import range
from Components.HalfComplex import HalfComplex
import math


def square(P, r):
    """
    Komplex elem négyzetre emelése
    :param P: a komplex elem;
    :param r: a rend;
    :return: P négyzete;
    """
    #     0 or 1
    if P.isOne() or P.isZero():
        return P
    #     real number
    if P.img == 0:
        newReal = (P.real * P.real) % r
        return Fp2Element(newReal, 0)

    #     imaginary
    if P.real == 0:
        newReal = P.img * P.img
        newReal = Math.negate(newReal, r)
        return Fp2Element(newReal, 0)

    #     (a+bi)^2=(a+b)(a-b)+2abi
    newReal = ((P.real + P.img) * (P.real - P.img)) % r
    newImg = ((P.real + P.real) * P.img) % r
    return Fp2Element(newReal, newImg)


def conjugate(P, q):
    """
    Konjugált
    :param P: projektív térbeli pont
    :param q: prím, a véges test
    :return: Konjugált
    """
    return Fp2Element(P.real, Math.negate(P.img, q))


def multiply(Q1, Q2, q):
    """
    Projektív elemek szorzása
    :param Q1: projekív térbeli elem
    :param Q2: projektív térbeli elem
    :param q: prím, a véges test
    :return: p * q
    """
    if Q1.real == Q2.real and Q1.img == Q2.img:
        return square(Q1, q)

    ac = Q1.real * Q2.real
    bd = Q1.img * Q2.img
    real = (ac - bd) % q
    img = Q1.real + Q1.img
    img = img * (Q2.real + Q2.img)
    img = img - ac
    img = (img - bd) % q
    return Fp2Element(real, img)


def divide(P, Q, q):
    """
    Projekív térbeli elemek osztása
    :param P: projektív térbeli elem, osztó;
    :param Q: projektív térbeli elem, osztandó;
    :param q: prím, a véges test;
    :return: q/p;
    """
    conj = conjugate(Q, q)
    top = multiply(P, conj, q)
    bottom = (P.real ** 2 + P.img ** 2) % q
    real = (top.real * Math.modular_inverse(bottom, q)) % q
    img = (top.img * Math.modular_inverse(bottom, q)) % q
    return Fp2Element(real, img)


def complexPow(P, k, ec):
    """
    Complex pow
    :param P: Fp2Element
    :param k: exponent
    :param ec: EC object
    """
    inv = False
    if k < 0:
        k = abs(k)
        inv = True
    windowSize = 3
    tLength = 2 ** (windowSize - 1)
    u = P
    u2 = square(u, ec.q)
    t = [None] * tLength
    t[0] = u
    for i in range(1, tLength):
        t[i] = multiply(u2, t[i - 1], ec.q)
    windows = getWindows(k, windowSize)
    i = len(windows) - 1
    if windows[i][0] % 2 == 0:
        u = t[windows[i][0] // 2]
        u = multiply(u, t[0], ec.q)
    else:
        u = t[windows[i][0] // 2]
    i = i - 1
    while i >= 0:
        actual = windows[i]
        for j in range(0, actual[1]):
            u = square(u, ec.q)
        if actual[0] > 0:
            if windows[i][0] % 2 == 0:
                u = multiply(u, t[windows[i][0] // 2], ec.q)
                u = multiply(u, t[0], ec.q)
            else:
                u = multiply(u, t[windows[i][0] // 2], ec.q)
        i = i - 1
    if inv:
        return inverse(u)
    return u


def getWindows(exponent, windowSize):
    expBin = str(bin(exponent))[2:][::-1]
    currentWindow = ""
    elementSizeList = list()
    elementValueList = list()
    i = 0
    zw = False
    while i < len(expBin):
        if expBin[i] == "0":  # ZW
            currentWindow = currentWindow + expBin[i]
            zw = True
            i = i + 1
        else:  # NZW
            if zw:
                elementValueList.append(int(currentWindow[::-1], 2))
                elementSizeList.append(len(currentWindow))
                currentWindow = ""
                zw = False
            while len(currentWindow) < windowSize and i < len(expBin):
                currentWindow = currentWindow + expBin[i]
                i = i + 1
        if len(currentWindow) == windowSize or i == (len(expBin)):
            elementValueList.append(int(currentWindow[::-1], 2))
            elementSizeList.append(len(currentWindow))
            currentWindow = ""
            zw = False
    return list(zip(elementValueList, elementSizeList))


def inverse(P):
    return divide(Fp2Element(1, 0), P)
    # return divide(P, Fp2Element(1, 0))


def toHalfComplex(P, q):
    """
    Komplex átalakítása félkomplexxé
    :param P: a projektív térbeli pont
    :param q: prím, a véges test
    :return: HalfComplex(q, f.real)
    """
    return HalfComplex(q, P.real)
