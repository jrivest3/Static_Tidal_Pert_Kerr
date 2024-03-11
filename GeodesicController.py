#GeodesicController
import numpy as np
# import scipy
# from KerrB import *
from KerrFrequencies import *

class Geodesic:
    def __init__(self,M,a,p,e,x,r_init=None,th_init=None,rdot_init=None,thdot_init=None) -> None: #t_init and phi_init will be set to zero without loss of generality
        this=self
        self.M=M
        self.params=(a,p,e,x)
        self.name=f"Geodesic: {self.params}"
        self.a=a
        self.p=p
        self.e=e
        self.x=x
        
        

        def Delta(r):return r * r - 2 * r * M + a * a
        def Sigma(r,usq):return r * r + a * a * usq  #usq=cos(theta)^2

        self.zm=np.sqrt(1-x*x)
        self.r1=p/(1-e)
        self.r2=p/(1+e)
        zm,r1,r2=self.zm,self.r1,self.r2

        # Special case functions defined in KerrB
        this.En = KerrGeoEnergy(a, p, e, x)
        this.Lz = KerrGeoAngularMomentum(a, p, e, x, this.En)
        this.Q = KerrGeoCarterConstant(a, p, e, x, this.En, this.Lz)
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
        
        self.CoMs = [this.En, this.Lz, this.Q]
        AplusB = (2 * M) / (1 - this.En ** 2) - (this.r1 + this.r2) 
        AB = (a ** 2 * this.Q) / ((1 - this.En ** 2) * this.r1 * this.r2)
        self.r3 = (AplusB + np.sqrt(AplusB ** 2 - 4 * AB)) / 2 #(*Eq.(11)*)
        self.r4 = AB / this.r3
        self.zp = np.sqrt(this.Q / (1 - x * x)) if x*x<1 else self.zm
        self.RadialRoots = [this.r1, this.r2, this.r3, this.r4]
        self.PolarRoots = [this.zm, this.zp]
        self.Frequencies = Frequencies(a, p, e, x, this.En, this.Lz, this.Q, this.r1, this.r2, this.r3, this.r4)
        # self.isResonant = False 
        #I may wan't a dictionary of Geodesics with resonances as Keys 
        # perhaps Resonant or closed Geodesics should be a subclass

        En,Lz,Q=this.En, this.Lz, this.Q
        Phisq = Lz *Lz / (1 - En *En)
        q = Q / (1 - En *En) # is mu^2 1 or -1?
        U0sq = 1-x*x # zm^2
        asqU1sq = ((a *a + q + Phisq) + np.sqrt((a *a + q + Phisq) ** 2 - 4 * a *a * q)) / (2) #* a**2); #u0 and u1 just need to be swapped. now u1^2>1
        
        self.r0=r1 if r_init is None else r_init
        self.th0=np.arccos(zm) if th_init is None else th_init
        r0,th0=self.r0,self.th0
        self.y0=np.array([0,r0,th0,0]) 
        
        usq=np.cos(th0)*np.cos(th0) # double check if this is current theta or theta_min
        Delt0=Delta(r0)
        Sigm0=Sigma(r0,usq)
        
        self.rdot0=0 if rdot_init is None else rdot_init
        self.thdot0=0 if thdot_init is None else thdot_init
        self.ydot0=np.array([(a * (Lz - a * En * (1 - usq)) + (r0 *r0 + a *a) * ((r0 *r0 + a *a) * En - Lz * a) / Delt0) / Sigm0,self.rdot0,self.thdot0,(2 * M * r0 * a * En + (Sigm0 - 2 * M * r0) * Lz / (1 - usq)) / (Delt0 * Sigm0)])

        self.ICs=np.array([self.y0,self.ydot0])
        
        self.ProperTime=np.array([0])
        self.PhaseSpaceTrajectory={"positions":np.array([self.y0]),"velocities":np.array([self.ydot0])}

        self.info={}

    # def __repr__(self) -> str:
    #     pass
    def __str__(self) -> str:
        return self.name

    def UpdateTrajOneStep(self,t:float,Pos_array:np.ndarray,Vel_array:np.ndarray)->None:
        # try:
        # if array.shape is one step
        np.append(self.ProperTime,t,axis=0)
        np.append(self.PhaseSpaceTrajectory["positions"],Pos_array,axis=0)
        np.append(self.PhaseSpaceTrajectory["velocities"],Vel_array,axis=0)
        # except ValueError:
        # else: .join instead of append, or append in loop?

    # def findResonance():
    #     if this geodesic has a resonance, return the resonance and set isResonant to True

    def Perturb(self,p_Model:str,p_name:str,z_array:np.ndarray)->None: # Should perturbed Geodesics be a subclass?
        self.name=self.name+" Perturbed by " + p_name + f": {z_array}"
        self.perturber=p_name
        self.perturbation_type=p_Model
        self.zs=z_array


