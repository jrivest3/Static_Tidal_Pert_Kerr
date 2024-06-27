#GeodesicController
import numpy as np
# import scipy
# from KerrB import *
from KerrFrequencies import *

# from CristhValuesArrays import TotalGamma




class Geodesic_Set:
    def __init__(self,M,a,p,e,x,r_psi_Nsteps=10,th_psi_Nsteps=10,phi_Nsteps=10) -> None: #,rdot_init=0,thdot_init=0 t_init and phi_init will be set to zero without loss of generality
        this=self
        self.M=M
        self.params=(a,p,e,x)
        self.name=f"Geodesics: {self.params}"
        self.a=a
        self.p=p
        self.e=e
        self.x=x
        
        self.xi0s= np.linspace(0,2*np.pi,r_psi_Nsteps,endpoint=False)
        self.chi0s= np.linspace(0,2*np.pi,th_psi_Nsteps,endpoint=False)
        self.phi0s= np.linspace(0,2*np.pi,phi_Nsteps,endpoint=False)
        # self.init_phases=np.meshgrid()


        
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
        self.Base_Freqs = Frequencies(a, p, e, x, this.En, this.Lz, this.Q, this.r1, this.r2, this.r3, this.r4)
        
        # self.isResonant = False 
        #I may want a dictionary of Geodesics with resonances as Keys 
        # perhaps Resonant or closed Geodesics should be a subclass

        class Geodesic:
            tau_array=None
            Trajectory=None # shape=(size(t_eval),4,dim)
            # will come from a SpaceTime method call
            def __init__(self,ICs) -> None:
                self.name=this.name+f", {ICs=}"
                self.ICs= ICs
                # self.tau_array=
                # self.Trajectory=None
                
                # self.ProperTime=np.array([0])
                # self.PhaseSpaceTrajectory={"positions":np.array([self.y0]),"velocities":np.array([self.ydot0]),"accelerations":np.array([]),"jerks":np.array([])}

                # self.info={}

            def __str__(self) -> str:return self.name 

            def check_Traj(self):
                if self.Trajectory is None: 
                    print('No Trajectory exists for '+self.name) 
                    return -1
                else: return 0
            # def cal_Approx_Freqs(cls): # Need phase trajectories
            #     # Only makes sense if perturbed
            #     if cls.check_Traj()!=0: pass
            #     tList=cls.tau_array
            #     Traj=cls.Trajectory # shape=?
            #     assert tList is not None and Traj is not None
            #     cls.Approx_Freqs= np.polyfit(tList, Traj[0], 1)[0], np.polyfit(tList, etaList, 1)[0], np.polyfit(tList, phiList, 1)[0]
            # def rootFindParametersFromFrequencies(self):#,Omega_xi, Omega_eta, Omega_phi, init_apex):
            #     try:
            #         return rootFindParametersFromFrequencies(*self.Approx_Freqs, self.params)
            #     # def f(apex):
            #     #     omega_xi, omega_eta, omega_phi = Frequencies(*apex)
            #     #     return (omega_xi - OmegaR)**2 + (omega_eta - OmegaT)**2 + (omega_phi - OmegaP)**2

            #     # return opt.minimize(f, np.array([a,p,e,x]),
            #     #         bounds = [(0,.999), (1.6, float('inf')),(0, .999), (-1, 1)],
            #     #         tol=1e-16)
            # # @classmethod
            # def find_Unperturbed_Counterpart_Trajectory(cls):
            #     # Only makes sense if perturbed
            #     if cls.check_Traj()!=0: pass
            #     cls.Nearest_Matching_Traj={'Freq.s':this.,'params':}
            #     # call rootFindParametersFromFrequencies with 
            #     # new Traj has same ICs, not the same params

        En,Lz,Q=this.En, this.Lz, this.Q
        Phisq = Lz *Lz / (1 - En *En)
        q = Q / (1 - En *En) # is mu^2 1 or -1?
        # U0sq = 1-x*x # zm^2
        asqU1sq = ((a *a + q + Phisq) + np.sqrt((a *a + q + Phisq) ** 2 - 4 * a *a * q)) / (2) #* a**2); #u0 and u1 just need to be swapped. now u1^2>1
        


        self.Trajectories={}
        for xi0 in self.xi0s:
            for chi0 in self.chi0s:
                for phi0 in self.phi0s:
                    r0= p/(1+e*np.cos(xi0)) # r2 +(r1-r2)*np.cos(xi0)
                    th0= np.arccos(zm)*np.cos(chi0)
                    y0=np.array([0,r0,th0,phi0]) 
                    usq=np.cos(th0)*np.cos(th0) # double check if this is current theta or theta_min
                    Delt0=Delta(r0)
                    Sigm0=Sigma(r0,usq)

                    # Phisq = Lz *Lz / (1 - En *En)
                    # q = Q / (1 - En *En) # is mu^2 1 or -1?
                    # U0sq = zm *zm; #((a**2 + q + Phisq) - np.sqrt((a**2 + q + Phisq)**2 - 4* a**2* q))/(2.* a**2);// this is zm^2, right?
                    # asqU1sq = ((a *a + q + Phisq) + np.sqrt((a *a + q + Phisq) *(a *a + q + Phisq) - 4 * a *a * q)) / (2.)#* a**2); #u0 and u1 just need to be swapped. now u1^2>1

                    sinchi = np.sin(chi0)
                    # sinchisq = sinchi * sinchi
                    # usq = U0sq * sinchisq
                    tdot0 = (a * (Lz - a * En * (1 - usq)) + (r0 *r0 + a *a) * ((r0 *r0 + a *a) * En - Lz * a) / Delt0) / Sigm0
                    # xidot0 = np.sqrt((1 - En *En) * r1 * r2 * (r0 - this.r3) * (r0 - this.r4)) / Sigm0 / r0 #np.sqrt((mu**2 - En**2)*np.cos[xi]*(r - r3)*(r - r4))/Sig;
                    chidot0 = np.sqrt((1 - En *En) * (asqU1sq - usq * a *a)) / Sigm0
                    # phidot0 = (2. * M * r0 * a * En + (Sig - 2. * M * r0) * Lz / (1. - usq)) / (Delt0 * Sigm0)
                    
                    
                    rdot0= r0*e*np.sin(xi0)*np.sqrt((1-En*En)/(1-e*e)*(r0-self.r3)*(r0-self.r4))/ Sigm0
                    thdot0=-sinchi/np.tan(th0)*chidot0
                    phdot0=(2 * M * r0 * a * En + (Sigm0 - 2 * M * r0) * Lz / (1 - usq)) / (Delt0 * Sigm0)
                    ydot0=np.array([tdot0,rdot0,thdot0,phdot0])
                    # self.yd0=

                    ICs=np.array([y0,ydot0])
                    self.Trajectories[(xi0,chi0,phi0)]=Geodesic(ICs)
                    
                    # self.Geodesics[(xi0,chi0,phi0)].ProperTime=np.array([0])
                    # self.Geodesics[(xi0,chi0,phi0)].PhaseSpaceTrajectory={"positions":np.array([y0]),"velocities":np.array([ydot0]),"accelerations":np.array([]),"jerks":np.array([])}

                    # self.Geodesics[(xi0,chi0,phi0)].info={}

        

    # def __repr__(self) -> str:
    #     pass
    def __str__(self) -> str: # should also give some information about the Trajectories
        return self.name #+'(xi_0,chi_o,phi_0)='+str(self.init_phase)

    # def UpdateTrajOneStep(self,t:float,Pos_array:np.ndarray,Vel_array:np.ndarray,Acc_array:np.ndarray,Jerk_array:np.ndarray)->None:
    #     # try:
    #     # if array.shape is one step
    #     np.append(self.ProperTime,t,axis=0)
    #     np.append(self.PhaseSpaceTrajectory["positions"],Pos_array,axis=0)
    #     np.append(self.PhaseSpaceTrajectory["velocities"],Vel_array,axis=0)
    #     np.append(self.PhaseSpaceTrajectory["accelerations"],Acc_array,axis=0)
    #     np.append(self.PhaseSpaceTrajectory["jerks"],Jerk_array,axis=0)
    #     # if Acc_array is not None: 
    #     #     if "accelerations" not in self.PhaseSpaceTrajectory: self.PhaseSpaceTrajectory
    #     #     np.append(self.PhaseSpaceTrajectory["accelerations"],Vel_array,axis=0)
    #     # except ValueError:
    #     # else: .join instead of append, or append in loop?

    # # def findResonance():
    # #     if this geodesic has a resonance, return the resonance and set isResonant to True

    # # def Perturb(self,p_Model:str,p_name:str,z_array:np.ndarray)->None: # Should perturbed Geodesics be a subclass?
    # #     self.name=self.name+" Perturbed by " + p_name + f": {z_array}"
    # #     self.perturber=p_name
    # #     self.perturbation_type=p_Model
    # #     self.zs=z_array

    