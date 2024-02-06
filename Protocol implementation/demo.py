# echo-server.py

from builtins import str
import hashlib
import hmac
import random
import timeit
import time
import socket
import struct
import random
import datetime
from Components import Complex
from Components import EC, HalfComplex
from Components import Math
from Components.BonehFranklin import BonehFranklin_stream, BonehFranklin_block, point_to_byte, byte_to_point
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
# Local user list and global blacklist
users = []
blacklist = []
pool = []


ec = EC.EC()
ec.basepoint = ec.at(None)


now=datetime.datetime.now()
logs_time = open("log_time_"+str(now)+".txt","w")
logs_message = open("log_message_"+str(now)+".txt","w")


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

def registry_rnd_car(lp,gamma,v = 1,car_id = None):
    start = timeit.default_timer()
    Qw = ec.at(lp)
    gQw = ec.mulJ(Qw, gamma)
    stop = timeit.default_timer()
    logs_time.write(str(car_id)+"#registry_rnd_car#"+str(stop-start)+"\n")
    if v:
        print("LP             :", lp)
        print("Qw            x:", Qw.toString().split()[0])
        print("              y:", Qw.toString().split()[1])
        print("---------------------------------------------------------------------------------------------------------------")
        print("gQw           x:", gQw.toString().split()[0])
        print("              y:", gQw.toString().split()[1])


    return Qw,gQw



def obu_compute_cs_1(Qv, gammaQv, Qr, ec, v = 1,car_id = None):
    start = timeit.default_timer()

    s = random.getrandbits(128)
    A1 = bilinear_map(gammaQv, Qr, ec)
    sgammaQv = ec.mulJ(gammaQv, s)

    Qv_byte = point_to_byte(Qv)
    sQgv_byte = point_to_byte(sgammaQv)

    msg_byte = Qv_byte+A1.real.to_bytes(64)+sQgv_byte

    M1 = BonehFranklin_block.encryp(BonehFranklin_block, msg_byte, Qr, gammaP, P, ec)

    stop = timeit.default_timer()

    logs_time.write(str(car_id)+"#obu_compute_cs_1#"+str(stop-start)+"\n")
    msg = Qv.toString() + "/" + A1.toString() + "/"  + sgammaQv.toString()
    if v:
        print("s           =",s)
        print("---------------------------------------------------------------------------------------------------------------")
        print("A1          =",A1.toString())
        print("---------------------------------------------------------------------------------------------------------------")
        print("sgammaQv      x:",sgammaQv.toString().split()[0])
        print("              y:",sgammaQv.toString().split()[1])
        print("---------------------------------------------------------------------------------------------------------------")
        print("msg         =",msg)

    return s, M1

def rsu_compute_cs(xi, gammaQr, M1, v = 1,car_id = None):
    start = timeit.default_timer()

    decM1 = BonehFranklin_block.decrypt(BonehFranklin_block, M1, gammaQr, ec)
    Qv = byte_to_point(decM1[:128],ec)
    A1_obu = int.from_bytes(decM1[128:192])
    sgQv = byte_to_point(decM1[192:],ec)


    try:
        blacklist.index(Qv.toString())
        RL = True
    except ValueError:
        RL = False


    A1_rsu = bilinear_map(Qv, gammaQr, ec)


    validity = A1_rsu.real == A1_obu

    xi_s_gammaQv = ec.mulJ(sgQv, xi)

    stop = timeit.default_timer()
    logs_time.write(str(car_id)+"#rsu_compute_cs#"+str(stop-start)+"\n")
    if v:
        print("Dec M1      = ", Qv.toString() + "/" + str(A1_obu) + "/"  + sgQv.toString())
        print("---------------------------------------------------------------------------------------------------------------")
        print("Rec Qv      = x:", Qv.toString().split()[0])
        print("              y:", Qv.toString().split()[1])
        print("---------------------------------------------------------------------------------------------------------------")
        print("Rec A1      = ", A1_obu)
        print("")
        print("On revoc    =", RL)
        print("RSU A1      = ", A1_rsu.real)
        print("Validity    =", validity)
        print("---------------------------------------------------------------------------------------------------------------")
        print("xisgQv      = x:", xi_s_gammaQv.toString().split()[0])
        print("              y:", xi_s_gammaQv.toString().split()[1])

    return xi_s_gammaQv,RL, validity

def obu_compute_cs_2(Qv, gammaQv, xi_s_gammaQv, s, ec, Qr, xiQr, v = 1,car_id = None):
    start = timeit.default_timer()

    s_inv = Math.modular_inverse(s, ec.r)

    xi_gammaQv = ec.mulJ(xi_s_gammaQv, s_inv)

    chk_1 = bilinear_map(xi_gammaQv, Qr, ec).real
    chk_2 = bilinear_map(gammaQv, xiQr, ec).real

    validity = chk_1 == chk_2

    stop = timeit.default_timer()

    logs_time.write(str(car_id)+"#obu_compute_cs_2#"+str(stop-start)+"\n")
    if v:
        print("s            =",s)
        print("s-1          =",s_inv)

        print("xigQv        = x:", xi_gammaQv.toString().split()[0])
        print("               y:", xi_gammaQv.toString().split()[1])
        print("b(xigQv ,Qr) =",str(chk_1.real))
        print("b(gQv ,xiQr) =",str(chk_2.real))
        print("validity     =",validity)
    return xi_gammaQv, validity

def incident_report(Qv,xigQv, v = 1, M = None,car_id = None):
    start = timeit.default_timer()
    a = random.getrandbits(128)
    Aid = ec.mulJ(Qv, a)
    A1 = ec.mulJ(xigQv, a)
    stop = timeit.default_timer()
    logs_time.write(str(car_id)+"#incident_report"+"#"+str(stop-start)+"\n")
    if v:
        print("Aid          = x:", Aid.toString().split()[0])
        print("               y:", Aid.toString().split()[1])
        print("A1           = x:", A1.toString().split()[0])
        print("               y:", A1.toString().split()[1])

    return Aid, A1, M


def verify_incident_report(Aid, A1,Qw, xigQw, ec, v = 1,car_id = None):
    start = timeit.default_timer()
    chk_1 = bilinear_map(Aid, xigQw, ec).toString()
    chk_2 = bilinear_map(A1, Qw, ec).toString()
    validity = chk_1 == chk_2
    stop = timeit.default_timer()
    logs_time.write(str(car_id)+"#verify_incident_report_2"+"#"+str(stop-start)+"\n")
    if v:
        print("Qw           = x:", Qw.toString().split()[0])
        print("               y:", Qw.toString().split()[1])
        print("xigQw           = x:", xigQw.toString().split()[0])
        print("               y:", xigQw.toString().split()[1])
        print("Aid          = x:", Aid.toString().split()[0])
        print("               y:", Aid.toString().split()[1])
        print("A1           = x:", A1.toString().split()[0])
        print("               y:", A1.toString().split()[1])
        print("b(Aid, xigQw)     =", chk_1)
        print("b(A1, Qw)         =", chk_2)
        print("Validity          =",validity)
    return validity



print("###############################################################################################################")
print("############################################# Random Car Registry #############################################")
print("###############################################################################################################\n")

Qv,gQv = registry_rnd_car("idx411",gamma, v = 1, car_id= 1)
print("\n###############################################################################################################\n")
Qw,gQw = registry_rnd_car("idx412",gamma, v = 1, car_id= 2)
print("\n###############################################################################################################\n")
Qr,gQr = registry_rnd_car("rsu412",gamma, v = 1, car_id= 0)
xi = random.getrandbits(128)
xiQr =ec.mulJ(Qr, xi)
print("\n###############################################################################################################\n")

print("###############################################################################################################")
print("#############################################        Params       #############################################")
print("###############################################################################################################\n")
print("Qv          = x:", Qv.toString().split()[0])
print("              y:", Qv.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("gQv         = x:", gQv.toString().split()[0])
print("              y:", gQv.toString().split()[1])
print("\n###############################################################################################################\n")
print("Qw          = x:", Qv.toString().split()[0])
print("              y:", Qv.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("gQw         = x:", gQv.toString().split()[0])
print("              y:", gQv.toString().split()[1])
print("\n###############################################################################################################\n")
print("Qr          = x:", Qr.toString().split()[0])
print("              y:", Qr.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("gQr         = x:", gQr.toString().split()[0])
print("              y:", gQr.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("xi          = ", xi)
print("xiQr        = x:", xiQr.toString().split()[0])
print("              y:", xiQr.toString().split()[1])
print("\n###############################################################################################################\n")

print("\n###############################################################################################################")
print("############################################# Communication setup #############################################")
print("###############################################################################################################\n")
print("############################################# Step 1: OBU compute #############################################")

s, M1 = obu_compute_cs_1(Qv, gQv, Qr, ec, v = 1,car_id = 1)

print("\n############################################# Step 2: RSU compute #############################################")
xi_s_gammaQv,RL, validity = rsu_compute_cs(xi, gQr, M1, v = 1,car_id = 0)
print("\n############################################ Step 3: OBU Finalize ############################################")
xigQv, validity = obu_compute_cs_2(Qv, gQv, xi_s_gammaQv, s, ec, Qr, xiQr, v = 1,car_id = 1)
print("\n###############################################################################################################")
print("#############################################   Incident report   #############################################")
print("###############################################################################################################\n")
print("#############################################     OBU compute     #############################################")

print("Report an incident")
MSG = incident_report(Qv,xigQv, M = (10.0,10.0,1), car_id = 1)

print("A car enters RSU square, and regist")

s_w, M1_w = obu_compute_cs_1(Qw, gQw, Qr, ec, v = 0,car_id = 2)
xi_s_gammaQw,RL, validity_obu = rsu_compute_cs(xi, gQr, M1_w, v = 0,car_id = 2)
xigQw, validity_key = obu_compute_cs_2(Qw, gQw, xi_s_gammaQw, s_w, ec, Qr, xiQr, v = 0,car_id = 2)



print("The car:")
print("    Qw          = x:", Qw.toString().split()[0])
print("                  y:", Qw.toString().split()[1])
print("---------------------------------------------------------------------------------------------------------------")
print("    gQw         = x:", gQw.toString().split()[0])
print("                  y:", gQw.toString().split()[1])

print("    OBU is on the revoc list:", RL)
print("    OBU has valid key       :",validity_obu)
print("    Xi and key has integrity:",validity_key)
print("    ")



print("The car reads the message:")

print("    Rec Aid      = x:", MSG[0].toString().split()[0])
print("                   y:", MSG[0].toString().split()[1])
print("    Rec A1       = x:", MSG[1].toString().split()[0])
print("                   y:", MSG[1].toString().split()[1])
print("    Rec MSG      = x_pos:", MSG[2][0])
print("                   y_pos:", MSG[2][1])
print("                   val:", MSG[2][2])
print("The car compute:")
validity = verify_incident_report(MSG[0], MSG[1],Qw, xigQw, ec, v = 1,car_id = 2)
print(f"The validity of message is: {validity}")







