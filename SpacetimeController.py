#SpacetimeController

import time

import numpy as np
import spherical
import quaternionic

from CristhValuesArrays import TotalGamma
import GeodesicController 

from scipy.integrate import solve_ivp as ODE
from scipy.special import ellipk,ellipe

sqrtof2=np.sqrt(2)
sqrtof3=np.sqrt(3)
sqrtof5=np.sqrt(5)
sqrtof6=sqrtof2*sqrtof3

start=time.perf_counter()
restart=start

class Perturber:
    def __init__(self,zs=None,Eij=None,name=None,Model='Pre-set Perturbation') -> None:
        self.Model=Model
        self.name='Pre-set Perturbation' if name is None else name
        self.z_array=zs #if zs is not None else np.zeros(5)
        self.Eij= Eij #if Eij is not None else np.zeros((3,3))

    def set_name(self,name):
        #try:
        if isinstance(name,int): self.name=self.name+' '+str(name)
        elif isinstance(name,str): self.name=name
        else: self.name=str(name)

    def make_zs_from_Eij(self)->None:
        if self.zs is None and self.Eij is not None:
            # try: exception for shape of Eij
            Eij=self.Eij
            z0,z1,z2=-2*(Eij[0,0]+Eij[1,1]),-2*(Eij[0,2]-complex(0,1)*Eij[1,2]),-2*(Eij[0,0]-Eij[1,1])+complex(0,4)*Eij[0,1]
            self.z_array=np.array([np.conj(z2),np.conj(z1),z0,z1,z2])
        else: print(f"These z's are already {self.zs}")

    def add(self,a,b):
        new_zs=a.z_array+b.z_array
        new_Eij=a.Eij+b.Eij
        return Perturber(new_zs,new_Eij,'Net Perturbation','Net Perturbation')
    
    def rotate_zs(self,rot_theta,rot_phi)->None:
        pass

    def make_zs_from_params(self)->None:
        pass

class PointMass(Perturber):
    def __init__(self,theta_p,phi_p,Epsilon=None,Mass=None,Distance=None, zs=None, Eij=None,name='Point Mass',Model='Point Mass') -> None:
        super().__init__(zs, Eij,name,Model)
        self.mass=Mass
        self.distance=Distance
        self.epsilon= Epsilon if Epsilon is not None else Mass/(Distance*Distance*Distance) # central mass M=1, M^2*Mp/Rp^3
        self.theta=theta_p
        self.phi=phi_p

        if self.z_array is None:
            if self.Eij is None:
                z0,z1,z2=( -1 + -3 * np.cos( 2 * theta_p ) ),6 * np.cos( theta_p ) * ( np.cos( \
                    phi_p ) + complex( 0,-1 ) * np.sin( phi_p ) ) * np.sin( theta_p ),6 * ( \
                    np.cos( 2 * phi_p ) + complex( 0,-1 ) * np.sin( 2 * phi_p ) ) * ( \
                    np.sin( theta_p )*np.sin( theta_p ) ) #May be faster to write in terms of zm and convert trig to exp
                self.z_companion= np.array([np.conj(z2),np.conj(z1),z0,z1,z2]) # The Norm of z_comp ranges [4, 2*sqrt(19) ~8.7] for theta_p =[0,pi/2]
                self.z_array=self.epsilon*self.z_companion
            else: self.make_zs_from_Eij()
        
class AccretionRing(Perturber):
    def __init__(self,zm,psi,Epsilon=None,Mass=None,Distance=None, zs=None, Eij=None,name='Accretion Ring',Model='Accretion Ring') -> None:
        super().__init__(zs, Eij,name,Model)
        # self.Model=
        # self.name=
        self.mass=Mass
        self.distance=Distance
        self.epsilon= Epsilon if Epsilon is not None else Mass/(Distance*Distance*Distance) # central mass M=1, M^2*Mp/Rp^3
        self.inclination=zm
        self.rotation=psi

        if self.z_array is None:
            if self.Eij is None:
                z0,z1,z2=   2 + -3 * ( zm*zm ), \
                            4/( np.pi )/( zm ) * \
                                ( ( -1 + 2 * ( zm*zm ) ) * ellipe( ( zm*zm ) ) + -1 * ( -1 + ( zm*zm ) ) \
                                * ellipk( ( zm*zm ) ) ) * ( np.cos( psi ) + complex( 0,-1 ) * np.sin( psi ) ), \
                            -3/2 * ( zm*zm ) * ( np.cos( 2 * psi ) + complex( 0,-1 ) * np.sin( 2 * psi ) )
                self.z_ring= np.array([np.conj(z2),np.conj(z1),z0,z1,z2])
                self.z_array=self.epsilon*self.z_ring
            else: self.make_zs_from_Eij()
        

class Spacetime:
    
    def __init__(self,a,M=1) -> None:
        thisST=self
        # Relavent sYlm's
        self.M=M
        self.a=a
        self.NetPerturbation=Perturber(Model='Flat Background',name='Flat Background')
        self.Perturbers=[]
        self.PointMasses=[]
        self.Rings=[]
        self.OtherSources=[]
        
    def GeodesicReconstructor(self,p,e,x,r_init,theta_init,phi_init):
        self.Geodesic=GeodesicController.Geodesic(self.a,p,e,x)

    def Add_Perturbation_Source(self,Model_Or_zs_Or_Eij,*kwargs)->None:
        NewSource=0
        if type(Model_Or_zs_Or_Eij) is str:
            if Model_Or_zs_Or_Eij=='point':
                NewSource=PointMass(kwargs)
                self.PointMasses.append(NewSource)
                self.Perturbers.append(NewSource)
                if NewSource.name=='Point Mass':NewSource.set_name(len(self.PointMasses))
            elif Model_Or_zs_Or_Eij=='ring':
                NewSource=AccretionRing(kwargs)
                self.Rings.append(NewSource)
                self.Perturbers.append(NewSource)
                if NewSource.name=='Accretion Ring':NewSource.set_name(len(self.Rings))
        elif isinstance(Model_Or_zs_Or_Eij,type(np.ndarray(5,dtype=np.complex128))):
            NewSource=Perturber(zs=Model_Or_zs_Or_Eij,*kwargs)
            self.OtherSources.append(NewSource)
            self.Perturbers.append(NewSource)
            if NewSource.name=='Pre-set Perturbation':NewSource.set_name(len(self.OtherSources))
        elif isinstance(Model_Or_zs_Or_Eij,type(np.ndarray((3,3)))):
            NewSource=Perturber(Eij=Model_Or_zs_Or_Eij,*kwargs)
            self.OtherSources.append(NewSource)
            self.Perturbers.append(NewSource)
            if NewSource.name=='Pre-set Perturbation':NewSource.set_name(len(self.OtherSources))
        
        #try:
        if NewSource==0:print(f"{Model_Or_zs_Or_Eij} is not a valid argument. Try 'point','ring', np.array(5,dtype=np.complex128), or np.array((3,3)).")
        else:
            self.NetPerturbation=self.NetPerturbation+NewSource
            if len(self.Perturbers)==1: self.NetPerturbation.set_name(self.Perturbers[0].name)
        #except:print(f"{Model_Or_zs_Or_Eij} is not a valid argument. Try 'point','ring', np.array(5,dtype=np.complex128), or np.array((3,3)).")
    
    def Christoffels(self,radius,theta,phi): return TotalGamma(self.a,radius,theta,phi,self.NetPerturbation.z_array)

    