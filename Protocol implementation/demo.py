# echo-server.py

from builtins import str
import hashlib
import hmac
import random
import timeit
import time

from Components import Complex
from Components import EC, HalfComplex
from Components import Math
from Components.BonehFranklin import BonehFranklin
from Components.ECPoint import ECPoint
from Components.Fp2Point import Fp2Point
from Components.Fp2Element import Fp2Element
from Components.TatePairing import TatePairing
from Components.HalfComplex import HalfComplex
from Crypto.Cipher import AES

from ecdsa import SigningKey, VerifyingKey, NIST384p


gamma = None
# Global parameters
P = None
gammaP = None
ec = None
# OBU
Qv = None
gammaQv = None
ownToken = None
RSUToken = None
t = 0
# RSU
Qr = None
gammaQr = None
xi = None
xiQr = None
# Local user list and global blacklist
users = []
blacklist = []


ec = EC.EC()
ec.basepoint = ec.at(None)

################################# gamma #################################
with open('params/gamma', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    gamma = int(content[0])
################################# P #################################
with open('params/P', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    coord = []
    for line in content:
        coord.append(int(line))
    P = ECPoint(ec, coord[0], coord[1])
################################# gammaP #################################
with open('params/gammaP', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    coord = []
    for line in content:
        coord.append(int(line))
    gammaP = ECPoint(ec, coord[0], coord[1])
################################# Qv #################################
with open('params/Qv', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    coord = []
    for line in content:
        coord.append(int(line))
    Qv = ECPoint(ec, coord[0], coord[1])
################################# gammaQv #################################
with open('params/gammaQv', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    coord = []
    for line in content:
        coord.append(int(line))
    gammaQv = ECPoint(ec, coord[0], coord[1])
################################# ownToken #################################
with open('params/ownToken', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    coord = []
    for line in content:
        coord.append(int(line))
    ownToken = ECPoint(ec, coord[0], coord[1])
################################# RSUToken #################################
with open('params/RSUToken', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    coord = []
    for line in content:
        coord.append(int(line))
    RSUToken = ECPoint(ec, coord[0], coord[1])
################################# Qr #################################
with open('params/Qr', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    coord = []
    for line in content:
        coord.append(int(line))
    Qr = ECPoint(ec, coord[0], coord[1])
################################# gammaQr #################################
with open('params/gammaQr', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    coord = []
    for line in content:
        coord.append(int(line))
    gammaQr = ECPoint(ec, coord[0], coord[1])
################################# xi #################################
with open('params/xi', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    xi = int(content[0])
################################# xiQr #################################
xiQr = ec.mulJ(Qr, xi)
################################# users #################################
with open('params/users', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    for user in content:
        users.append(HalfComplex(ec.q, int(user)))
################################# blacklist #################################
with open('params/blacklist', 'r') as filehandle:
    content = [current_place.rstrip() for current_place in filehandle.readlines()]
if content:
    for user in content:
        coord = user.split(" ")
        blacklist.append(Fp2Point(int(coord[0]), int(coord[1])))

def bilinear_map(gammaQv, Qr, ec):
    return TatePairing.computeF(TatePairing, gammaQv, Qr, ec)

def point_from_string(s):
    s = s.split(" ")
    return Fp2Point(int(s[0]), int(s[1]))

def registry_rnd_car(lp,gamma,v = 1):
    Qw = ec.at(lp)
    gQw = ec.mulJ(Qw, gamma)
    if v:
        print("Qw          = x:", Qw.toString().split()[0])
        print("              y:", Qw.toString().split()[1])
        print("---------------------------------------------------------------------------------------------------------------")
        print("gQw         = x:", gQw.toString().split()[0])
        print("              y:", gQw.toString().split()[1])


    return Qw,gQw

def obu_compute_cs_1(Qv, gammaQv, Qr, ec, v = 1):
    s = random.getrandbits(128)
    A1 = bilinear_map(gammaQv, Qr, ec)
    sgammaQv = ec.mulJ(gammaQv, s)
    msg = Qv.toString() + ";" + A1.toString() + ";"  + sgammaQv.toString()
    M1 = BonehFranklin.encryp(BonehFranklin, msg, Qr, gammaP, P, ec)
    if v:
        print("s           =",s)
        print("---------------------------------------------------------------------------------------------------------------")
        print("A1          =",A1.toString())
        print("---------------------------------------------------------------------------------------------------------------")
        print("sgammaQv    = x:",sgammaQv.toString().split()[0])
        print("              y:",sgammaQv.toString().split()[1])
        print("---------------------------------------------------------------------------------------------------------------")
        print("msg         =",msg)

    return s, M1

def rsu_compute_cs(xi, xiQr, gammaQr, M1, v = 1):

    decM1 = BonehFranklin.decrypt(BonehFranklin, M1[0], M1[1], gammaQr, ec)
    Qv_tmp, A1_obu, sgQv_tmp = decM1.split(";")
    Qv = point_from_string(Qv_tmp)
    sgQv = point_from_string(sgQv_tmp)

    try:
        blacklist.index(Qv.toString())
        RL = True
    except ValueError:
        RL = False


    A1_rsu = bilinear_map(Qv, gammaQr, ec).toString()


    validity = A1_rsu == A1_obu
    xi_s_gammaQv = ec.mulJ(sgQv, xi)
    if v:
        print("Dec M1      = ", decM1)
        print("---------------------------------------------------------------------------------------------------------------")
        print("Rec Qv      = x:", Qv.toString().split()[0])
        print("              y:", Qv.toString().split()[1])
        print("---------------------------------------------------------------------------------------------------------------")
        print("Rec A1      = ", A1_obu)
        print("")
        print("On revoc    =", RL)
        print("RSU A1      = ", A1_rsu)
        print("Validity    =", validity)
        print("---------------------------------------------------------------------------------------------------------------")
        print("xisgQv      = x:", xi_s_gammaQv.toString().split()[0])
        print("              y:", xi_s_gammaQv.toString().split()[1])

    return xi_s_gammaQv,RL, validity

def obu_compute_cs_2(Qv, gammaQv, xi_s_gammaQv, s, ec, Qr, xiQr, v = 1):
    s_inv = Math.modular_inverse(s, ec.r)
    xi_gammaQv = ec.mulJ(xi_s_gammaQv, s_inv)
    chk_1 = bilinear_map(xi_gammaQv, Qr, ec).toString()
    chk_2 = bilinear_map(gammaQv, xiQr, ec).toString()

    validity = chk_1 == chk_2
    if v:
        print("s            =",s)
        print("s-1          =",s_inv)

        print("xigQv        = x:", xi_gammaQv.toString().split()[0])
        print("               y:", xi_gammaQv.toString().split()[1])
        print("b(xigQv ,Qr) =",chk_1)
        print("b(gQv ,xiQr) =",chk_2)
        print("validity     =",validity)
    return xi_gammaQv, validity

def incident_report(Qv,xigQv, M, v = 1):
    a = random.getrandbits(128)
    Aid = ec.mulJ(Qv, a)
    A1 = ec.mulJ(xigQv, a)
    if v:
        print("Aid          = x:", Aid.toString().split()[0])
        print("               y:", Aid.toString().split()[1])
        print("A1           = x:", A1.toString().split()[0])
        print("               y:", A1.toString().split()[1])
        print("MSG          = x_pos:", M[0])
        print("               y_pos:", M[1])
        print("               val:", M[2])
    return Aid.toString(), A1.toString(), M

def verify_incident_report(MSG,Qw, gQw, xigQw, ec, v = 1):
    Aid = point_from_string(MSG[0])
    A1 = point_from_string(MSG[1])
    M = MSG[2]
    chk_1 = bilinear_map(Aid, xigQw, ec).toString()
    chk_2 = bilinear_map(A1, Qw, ec).toString()
    validity = chk_1 == chk_2
    if v:
        print("b(Aid, xigQw)     =", chk_1)
        print("b(A1, Qw)         =", chk_2)
        print("Validity          =",validity)
    return validity


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

print("###############################################################################################################")
print("############################################# Random Car Registry #############################################")
print("###############################################################################################################\n")

Qw,gQw = registry_rnd_car("IDX411",gamma)

print("###############################################################################################################")
print("#############################################        Params       #############################################")
print("###############################################################################################################\n")
print("Qv          = x:", Qv.toString().split()[0])
print("              y:", Qv.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("gQv         = x:", gammaQv.toString().split()[0])
print("              y:", gammaQv.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("Qr          = x:", Qr.toString().split()[0])
print("              y:", Qr.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("gQr         = x:", gammaQr.toString().split()[0])
print("              y:", gammaQr.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("xi          = ", xi)
print("xiQr        = x:", xiQr.toString().split()[0])
print("              y:", xiQr.toString().split()[1])
print("\n###############################################################################################################")
print("############################################# Communication setup #############################################")
print("###############################################################################################################\n")
print("############################################# Step 1: OBU compute #############################################")
s, M1 = obu_compute_cs_1(Qv, gammaQv, Qr, ec)
#print(s, M1)

print("\n############################################# Step 2: RSU compute #############################################")
xi_s_gammaQv,RL, validity_obu = rsu_compute_cs(xi, xiQr, gammaQr, M1)

print("\n############################################ Step 3: OBU Finalize ############################################")
xigQv,validity_rsu = obu_compute_cs_2(Qv, gammaQv, xi_s_gammaQv, s, ec, Qr, xiQr)


print("\n###############################################################################################################")
print("#############################################   Incident report   #############################################")
print("###############################################################################################################\n")
print("#############################################     OBU compute     #############################################")

MSG = incident_report(Qv,xigQv, (10.0,10.0,1))

print("\n###############################################################################################################")
print("############################################# Verify Incident report ############################################")
print("###############################################################################################################\n")
print("#############################################   OBU/RSU compute   #############################################")
print("The car:")
print("    Qw          = x:", Qw.toString().split()[0])
print("                  y:", Qw.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("    gQw         = x:", gQw.toString().split()[0])
print("                  y:", gQw.toString().split()[1])
print("The car enters RSU square, and regist:")

s_w, M1 = obu_compute_cs_1(Qw, gQw, Qr, ec,0)
xi_s_gammaQw,RL, validity_obu = rsu_compute_cs(xi, xiQr, gammaQr, M1,0)
xigQw,validity_key = obu_compute_cs_2(Qw, gQw, xi_s_gammaQw, s_w, ec, Qr, xiQr,0)

print("    OBU is on the revoc list:", RL)
print("    OBU has valid key       :",validity_obu)
print("    Xi and key has integrity:",validity_rsu)
print("    ")
print("The car reads the message:")

print("    Rec Aid      = x:", MSG[0].split()[0])
print("                   y:", MSG[0].split()[1])
print("    Rec A1       = x:", MSG[1].split()[0])
print("                   y:", MSG[1].split()[1])
print("    Rec MSG      = x_pos:", MSG[2][0])
print("                   y_pos:", MSG[2][1])
print("                   val:", MSG[2][2])
print("The car compute:")
validity = verify_incident_report(MSG,Qw, gQw, xigQw, ec)

