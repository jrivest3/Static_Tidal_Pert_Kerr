#KerrFrequencies.py

import math
from KerrB import *

def Frequencies(a, p, e, x, En0, Lz0, Q0, r10, r20, r30, r40, zm, zp, M=1):
    z1, z2, kr2, kTheta2, rp, rm, En = En0, Lz = Lz0, Q = Q0, r1 = r10, r2 = r20, r3 = r30, r4 = r40
    hr, hp, hm, Kkr, KkTheta, Ekr, EkTheta, Pihr, Pihp, Pihm, Piz1 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    CapitalUpsilonr, CapitalUpsilonTheta, CapitalUpsilonCurlyPhi, CapitalGamma, CapitalOmegar, CapitalOmegaTheta, CapitalOmegaCurlyPhi = 0, 0, 0, 0, 0, 0, 0

    if Q0 == 0:
        En, Lz, Q = KerrGeoConstantsOfMotion(a, p, e, x)
    if r40 == 0:
        r1, r2, r3, r4 = KerrGeoRadialRoots(a, p, e, x, En, Q)

    z1 = 1 - x * x
    if z1 == 0:
        z2 = a * a * (1 - En * En) + Lz * Lz
    else:
        z2 = Q / z1

    kr2 = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
    kTheta2 = z1 / z2 * a * a * (1 - En * En)
    rRadical = math.sqrt((M + a) * (M - a))
    rp = M + rRadical
    rm = M - rRadical
    hr = (r1 - r2) / (r1 - r3)
    hp = (r1 - r2) * (r3 - rp) / ((r1 - r3) * (r2 - rp))
    hm = (r1 - r2) * (r3 - rm) / ((r1 - r3) * (r2 - rm))
    Kkr = math.ellipk(kr2)
    KkTheta = math.ellipk(kTheta2)
    Ekr = math.ellipe(kr2)
    EkTheta = math.ellipe(kTheta2)
    Pihr = math.ellippi(hr, kr2)
    Pihp = math.ellippi(hp, kr2)
    Pihm = math.ellippi(hm, kr2)
    Piz1 = math.ellippi(z1, kTheta2)

    CapitalUpsilonr = (math.pi / 2) * math.sqrt((1 - En * En) * (r1 - r3) * (r2 - r4)) / Kkr
    CapitalUpsilonTheta = (math.pi / 2) * math.sqrt(z2) / KkTheta

    if a < M:
        CapitalGamma = 4 * M * M * En + En / (1 - En * En) * z2 * (1 - EkTheta / KkTheta) +(1 / Kkr) \
            * (En / 2 \
                * ((r3 * (r1 + r2 + r3) - r1 * r2) * Kkr + (r2 - r3) \
                * (r1 + r2 + r3 + r4) * Pihr + (r1 - r3) * (r2 - r4) * Ekr) \
                + 2 * M * En * (r3 * Kkr + (r2 - r3) * Pihr) \
                + 2 * M / (rp - rm) \
                * ( \
                ((4 * M * M * En - a * Lz) * rp - 2 * M * a * a * En) \
                / (r3 - rp) \
                * (Kkr - (r2 - r3) / (r2 - rp) * Pihp) \
                - ((4 * M * M * En - a * Lz) * rm - 2 * M * a * a * En) \
                / (r3 - rm) \
                * (Kkr - (r2 - r3) / (r2 - rm) * Pihm) \
                )
            )
        CapitalUpsilonCurlyPhi = Lz * (Piz1 / KkTheta) + a / (Kkr * (rp - rm)) * ((2 * M * En * rp - a * Lz) / (r3 - rp) * (Kkr - (r2 - r3) / (r2 - rp) * Pihp) - (2 * M * En * rm - a * Lz) / (r3 - rm) * (Kkr - (r2 - r3) / (r2 - rm) * Pihm))
    else:
        raise ValueError('a=1')

    CapitalOmegar = CapitalUpsilonr / CapitalGamma
    CapitalOmegaTheta = CapitalUpsilonTheta / CapitalGamma
    CapitalOmegaCurlyPhi = CapitalUpsilonCurlyPhi / CapitalGamma

    return [CapitalOmegar, CapitalOmegaTheta, CapitalOmegaCurlyPhi]


