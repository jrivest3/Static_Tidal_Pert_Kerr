#GeodesicController
import numpy as np
# import scipy
from KerrB import *
from KerrFrequencies import *

class Geodesic:
    def __init__(self,M,a,p,e,x) -> None:
        this=self
        self.M=M
        M=self.M
        self.params=[a,p,e,x]
        self.a=a
        self.p=p
        self.e=e
        self.x=x

        def Delta(r):return r * r - 2 * r + a * a
        def Sigma(r,usq):return r * r + a * a * usq  #usq=cos(theta)^2

        self.zm=np.sqrt(1-x*x)
        self.r1=p/(1-e)
        self.r2=p/(1+e)
        zm,r1,r2=self.zm,self.r1,self.r2

        # Special case functions defined in KerrB
        this.En = KerrGeoEnergy(a, p, e, x)
        this.Lz = KerrGeoAngularMomentum(a, p, e, x, this.En)
        this.Q = KerrGeoCarterConstant(a, p, e, x, this.En, this.L)
        # these will return 'undefined' if no special cases apply.
        if self.Q is None:
            # Set Constants of Motion
            Del1 = Delta(r1)
            Del2 = Delta(r2)
            d1 = (r1 ** 2 + a ** 2 * zm ** 2) * Del1
            d2 = (r2 ** 2 + a ** 2 * zm ** 2) * Del2
            f1 = r1 ** 4 + a ** 2 * (r1 * (r1 + 2) + zm ** 2 * Del1)
            f2 = r2 ** 4 + a ** 2 * (r2 * (r2 + 2) + zm ** 2 * Del2)
            g1 = 2 * a * r1
            g2 = 2 * a * r2
            h1 = r1 * (r1 - 2) + zm ** 2 / (1 - zm ** 2) * Del1
            h2 = r2 * (r2 - 2) + zm ** 2 / (1 - zm ** 2) * Del2
            Kappa = d1 * h2 - h1 * d2
            Epsilon = d1 * g2 - g1 * d2
            Rho = f1 * h2 - h1 * f2
            Eta = f1 * g2 - g1 * f2
            Zeta = g1 * h2 - h1 * g2

            # En, Lz, and Q belong to Geodesic, not the controller.
            if self.En is None:
                self.En = np.sqrt((Kappa * Rho + 2 * Epsilon * Zeta - x * 2 * np.sqrt(Zeta * (Zeta * Epsilon ** 2 + Rho * Epsilon * Kappa - Eta * Kappa ** 2) / x ** 2)) / (Rho ** 2 + 4 * Eta * Zeta))
            if self.Lz is None:
                self.Lz = (0 - self.En * g1 + x * np.sqrt((-d1 * h1 + self.En ** 2 * (g1 ** 2 + f1 * h1)) / x ** 2)) / h1
            self.Q = zm ** 2 * (a ** 2 * (1 - self.En ** 2) + self.Lz ** 2 / (1 - zm ** 2))
            this.CoMs = [this.En, this.Lz, this.Q]
            AplusB = (2 * M) / (1 - this.En ** 2) - (this.r1 + this.r2) 
            AB = (a ** 2 * this.Q) / ((1 - this.En ** 2) * this.r1 * this.r2)
            this.r3 = (AplusB + np.sqrt(AplusB ** 2 - 4 * AB)) / 2 #(*Eq.(11)*)
            this.r4 = AB / this.r3
            this.zp = np.sqrt(this.Q / (1 - x ** 2))
            this.RadialRoots = [this.r1, this.r2, this.r3, this.r4]
            this.PolarRoots = [this.zm, this.zp]
            this.Frequencies = Frequencies(a, p, e, x, this.En, this.Lz, this.Q, this.r1, this.r2, this.r3, this.r4)

