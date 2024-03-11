#KerrFrequencies.py

# import math
from KerrB import *
from mpmath import ellipk,ellipe,ellippi

def Frequencies(a, p, e, x, En=0, Lz=0, Q=0, r1=0, r2=0, r3=0, r4=0, zm=0, zp=0, M=1):

    if Q==0:
        En, Lz, Q = KerrGeoConstantsOfMotion(a, p, e, x)
    if r4==0:
        r1, r2, r3, r4 = KerrGeoRadialRoots(a, p, e, x, En, Q)

    z1 = 1 - x * x
    if z1 == 0:
        z2 = a * a * (1 - En * En) + Lz * Lz
    else:
        z2 = Q / z1

    kr2 = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
    kTheta2 = z1 / z2 * a * a * (1 - En * En)
    rRadical = np.sqrt((M + a) * (M - a))
    rp = M + rRadical
    rm = M - rRadical
    hr = (r1 - r2) / (r1 - r3)
    hp = (r1 - r2) * (r3 - rp) / ((r1 - r3) * (r2 - rp))
    hm = (r1 - r2) * (r3 - rm) / ((r1 - r3) * (r2 - rm))
    Kkr = ellipk(kr2)
    KkTheta = ellipk(kTheta2)
    Ekr = ellipe(kr2)
    EkTheta = ellipe(kTheta2)
    Pihr = ellippi(hr, kr2)
    Pihp = ellippi(hp, kr2)
    Pihm = ellippi(hm, kr2)
    Piz1 = ellippi(z1, kTheta2)

    CapitalUpsilonr = (np.pi / 2) * np.sqrt((1 - En * En) * (r1 - r3) * (r2 - r4)) / Kkr
    CapitalUpsilonTheta = (np.pi / 2) * np.sqrt(z2) / KkTheta

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


