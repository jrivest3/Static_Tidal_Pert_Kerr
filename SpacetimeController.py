#SpacetimeController

# import time

#import numpy as np
# import spherical
# import quaternionic

from GeodesicController import *

from scipy.integrate import solve_ivp as ODE
# from scipy.special import ellipk,ellipe
from CristhValuesArrays import TotalGamma

sqrtof2=np.sqrt(2)
sqrtof3=np.sqrt(3)
sqrtof5=np.sqrt(5)
sqrtof6=sqrtof2*sqrtof3
sqrtpiover5=np.sqrt(np.pi/5)

# start=time.perf_counter()
# restart=start

class Perturber:
    def __init__(self,zs=None,Eij=None,name='Pre-set Perturbation',Model='Pre-set Perturbation',Epsilon=None) -> None:
        self.Model=Model #Perhaps assume a flat backgroud if zs are not given
        self.name=name
        self.z_array=zs #if zs is not None else np.zeros(5)
        self.Eij= Eij #if Eij is not None else np.zeros((3,3))
        if zs is None and Eij is not None: self.make_zs_from_Eij()
        self.Epsilon=Epsilon #For 

    def set_name(self,name):
        #try:
        if isinstance(name,int): self.name=self.name+' '+str(name)
        elif isinstance(name,str): self.name=name
        else: self.name=str(name)

    def make_zs_from_Eij(self)->None:
        if self.z_array is None and self.Eij is not None:
            # try: exception for shape of Eij
            sqrtpiover5=np.sqrt(np.pi/5)
            Eij=self.Eij
            z0,z1,z2=-2*np.sqrt(6)*sqrtpiover5*(Eij[0,0]+Eij[1,1]),-4*sqrtpiover5*(Eij[0,2]-complex(0,1)*Eij[1,2]),2*sqrtpiover5*(Eij[0,0]-Eij[1,1]-complex(0,2)*Eij[0,1])
            self.z_array=np.array([np.conj(z2),-np.conj(z1),z0,z1,z2])
            # # Harmonic coordinates
            # z0,z1,z2=-2*(Eij[0,0]+Eij[1,1]),-2*(Eij[0,2]-complex(0,1)*Eij[1,2]),-2*(Eij[0,0]-Eij[1,1])+complex(0,4)*Eij[0,1]
            # self.z_array=np.array([np.conj(z2),np.conj(z1),z0,z1,z2])
        elif self.Eij is None: print("No Eij's have been given for this perturber.")
        else: print(f"These z's are already {self.z_array}")

    def __add__(self,a):
        new_zs=self.z_array+a.z_array if self.z_array is not None else a.z_array
        new_Eij=self.Eij+a.Eij if self.Eij is not None and a.Eij is not None else None
        assert new_zs is not None
        return Perturber(zs=new_zs,Eij=new_Eij)#,name='Net Perturbation',Model='Net Perturbation')

    # def _outer_decorator(y):
    #     def _decorator(foo):
    #         def magic(self, *args, **kwargs) :
    #             print("start magic")
    #             if self.z_array is None: # will there be a scope problem with self?
    #                 if self.Eij is None:
    #                     return foo(self, *args, **kwargs)
    #             else:
    #                 raise ValueError("x ({}) != y ({})".format(self.x, y))
    #             print("end magic")
    #         return magic

    #     return _decorator

    # @_outer_decorator(y=3)
    # def bar(self, *args, **kwargs) :
    #     print("normal call")
    #     print("args: {}".format(args))
    #     print("kwargs: {}".format(kwargs))

    #     return 27
    # def epsilon_decorator(self,z_func):
    #     def wrapper():
    #         if self.z_array is None: # will there be a scope problem with self?
    #             if self.Eij is None:
    #                 self.z_characteristic=z_func()         
    #                 if Epsilon is None:
    #                     if Mass is None or Distance is None:
    #                         raise ValueError("Need to give value for either Epsilon or both Mass and Distance.")
    #                         #Error
    #                     # try: 
    #                     else:
    #                         self.epsilon=Mass/(Distance*Distance*Distance) # central mass M=1, M^2*Mp/Rp^3
    #                     # except ValueError:
    #                 # try:       
    #                 self.z_array=self.epsilon*self.z_characteristic
    #                 # except:
    #             else: self.make_zs_from_Eij()
    #     return wrapper

    def add_and_normalize_zs(self,z2_array=None,epsilon=1)->None:
        '''
        Used to aid in fixing the the strength of the perturbation. 
        If no second array is given, z_array is normalized. 
        The value of espilon becomes the new magnitude of the array.
        The option epsilon='eps' can be given to set the magnitude to the value of the Perturber's Epsilon attribute, if it has one.
        '''
        try:
            assert isinstance(self.z_array,type(np.ndarray(5,dtype=np.complex128))) #type(self.z_array) is np.ndarray
            if z2_array is not None: self.z_array+=z2_array
            eps=epsilon 
            if eps=='eps': 
                if self.Epsilon is None: 
                    print('This Perturber has no Epsilon given. Norm reset to 1')
                    eps=1
                else: eps=self.Epsilon
            self.z_array=eps*self.z_array/np.linalg.norm(self.z_array)
        except:
            print('Error in add_and_normalize_zs: starting z_array=',self.z_array,' secondary array is z2_array=',z2_array, 'epsilon=',eps)
            if type(self.z_array) is np.ndarray: print('norm is ',np.linalg.norm(self.z_array))
        

    def rotate_zs(self,rot_theta,rot_phi)->None:
        pass

    def make_zs_from_params(self)->None:
        pass

# class RequestParams(TypedDict):
#     url: str
#     allow_redirects: bool


# def request(**kwargs: Unpack[RequestParams]) -> None:

class PointMass(Perturber):
    def __init__(self,theta_p,phi_p,Mass=None,Distance=None,Epsilon=None, zs=None, Eij=None,name='Point Mass',Model='Point Mass') -> None:
        super().__init__(zs, Eij,name,Model)
        self.mass=Mass
        self.distance=Distance
        self.epsilon= Epsilon 

        self.theta=theta_p
        self.phi=phi_p

        if self.z_array is None: # Maybe use a decorator?
            if self.Eij is None:
                CosTHp=np.cos(theta_p)
                SinTHp=np.sin(theta_p)
                expPHp=np.exp(1.j*phi_p)
                z0,z1,z2=2*sqrtof6*sqrtpiover5* ( 1 - 3 * CosTHp*CosTHp ), 12*sqrtpiover5 /expPHp *CosTHp*SinTHp, -6*sqrtpiover5 /expPHp/expPHp *SinTHp*SinTHp#Cartesion (BL?)
                Norm = 4*sqrtof6*sqrtpiover5
                # We want to fix the norm to be epsilon.
                self.z_characteristic= np.array([np.conj(z2),-np.conj(z1),z0,z1,z2])/Norm # Norm = 4*sqrt(6*pi/5)~7.8 for all theta and phi.
                # # Harmonic Coordinates 
                # z0,z1,z2=( -1 + -3 * np.cos( 2 * theta_p ) ),6 * np.cos( theta_p ) * ( np.cos( \ 
                # phi_p ) + complex( 0,-1 ) * np.sin( phi_p ) ) * np.sin( theta_p ),6 * ( \
                # np.cos( 2 * phi_p ) + complex( 0,-1 ) * np.sin( 2 * phi_p ) ) * ( \
                # np.sin( theta_p )*np.sin( theta_p ) ) 
                # z_companion= np.array([np.conj(z2),np.conj(z1),z0,z1,z2]) # The Norm of z_comp ranges [4, 2*sqrt(19) ~8.7] for theta_p =[0,pi/2].
                
                if Epsilon is None:
                    if Mass is None or Distance is None:
                        raise ValueError("Need to give value for either Epsilon or both Mass and Distance.")
                        #Error
                    # try: 
                    else:
                        # put the Norm back in when specific mass and distance are given
                        self.epsilon=Norm*Mass/(Distance*Distance*Distance) # central mass M=1, M^2*Mp/Rp^3
                    # except ValueError:
                # try:
                if isinstance(self.epsilon,float):      
                    self.z_array=self.epsilon*self.z_characteristic
                else:
                    print(self.epsilon,type(self.epsilon))
                    assert isinstance(self.epsilon,float)
                # except:
            else: self.make_zs_from_Eij()
        assert self.z_array is not None
        
class AccretionRing(Perturber):
    def __init__(self,inclination,Long_of_Asc_Node,Epsilon=None,Mass=None,Distance=None, zs=None, Eij=None,name='Accretion Ring',Model='Accretion Ring') -> None:
        super().__init__(zs, Eij,name,Model)
        # self.Model=
        # self.name=
        self.mass=Mass
        self.distance=Distance
        self.epsilon= Epsilon

        if inclination<0 or inclination>np.pi/2:raise ValueError("inclination must be in the range [0,pi/2).")
        if inclination==np.pi/2:raise ValueError("Polar rings are not yet implemented. Choose inclination<pi/2")
        self.inclination=inclination
        # self.theta_min=np.pi/2-inclination
        zm=np.cos(np.pi/2-inclination)
        zm= 0 if np.isclose(zm,0,.00001,.00001) else zm
        self.Long_of_Asc_Node=Long_of_Asc_Node # angle from phi=0
        psi=Long_of_Asc_Node+np.pi/2 # value of phi at max inclination(zm)


        if self.z_array is None:
            if self.Eij is None:
                expPsi=np.exp(1.j*psi)
                if zm*zm==1:raise ValueError("Polar rings are not yet implemented. ")
                elif zm==0: 
                    z0,z1,z2= 2*sqrtof6*sqrtpiover5,0,0
                    zchar=np.array([np.conj(z2),-np.conj(z1),z0,z1,z2])
                    Norm=2*sqrtof6*sqrtpiover5 
                else: 
                    z0,z1,z2=   sqrtof6*sqrtpiover5*( 2 + -3 * zm*zm ), \
                            8*sqrtpiover5/( np.pi )/( zm ) * \
                                ( ( -1 + 2 * ( zm*zm ) ) * float(ellipe( ( zm*zm ) )) - ( -1 + ( zm*zm ) ) * float(ellipk( ( zm*zm ) )) ) \
                                /expPsi, \
                            3/2 * sqrtpiover5 * ( zm*zm ) /expPsi/expPsi
                    zchar=np.array([np.conj(z2),-np.conj(z1),z0,z1,z2])
                    Norm=np.linalg.norm(zchar) # ~3.84 near pole, peaks ~4.01 at theta~0.46276 (just over 7*pi/48, 7.070452*pi/48), 2*sqrtof6*sqrtpiover5 ~3.88 at the equator
                # We want to fix the norm to be epsilon.
                self.z_characteristic= np.array([np.conj(z2),-np.conj(z1),z0,z1,z2])/Norm
                # # Harmonic Coordinates
                # z0,z1,z2=   2 + -3 * ( zm*zm ), \
                #             4/( np.pi )/( zm ) * \
                #                 ( ( -1 + 2 * ( zm*zm ) ) * ellipe( ( zm*zm ) ) + -1 * ( -1 + ( zm*zm ) ) \
                #                 * ellipk( ( zm*zm ) ) ) * ( np.cos( psi ) + complex( 0,-1 ) * np.sin( psi ) ), \
                #             -3/2 * ( zm*zm ) * ( np.cos( 2 * psi ) + complex( 0,-1 ) * np.sin( 2 * psi ) )
                # self.z_ring= np.array([np.conj(z2),np.conj(z1),z0,z1,z2])
                if Epsilon is None:
                    if Mass is None or Distance is None:
                        raise ValueError("Need to give value for either Epsilon or both Mass and Distance.")
                        #Error
                    # try: 
                    else:
                        # put the Norm back in when specific mass and distance are given
                        self.epsilon=Norm*Mass/(Distance*Distance*Distance) # central mass M=1, M^2*Mp/Rp^3
                    # except ValueError:
                # try:
                assert isinstance(self.epsilon,float)       
                self.z_array=self.epsilon*self.z_characteristic
                # except ValueError:
                #     assert isinstance(self.epsilon,float)
                #     print(f"epsilon is still {self.epsilon}")
            else: self.make_zs_from_Eij()
        assert self.z_array is not None
        

class Spacetime:
    
    def __init__(self,a,M=1) -> None:
        thisST=self
        # Relavent sYlm's
        self.M=M
        self.a=a
        self.NetPerturbation=Perturber(Model='Flat Background',name='Flat Background')
        self.PointMasses=[]
        self.Rings=[]
        self.OtherSources=[]
        self.Perturbers={"NetPerturbation":self.NetPerturbation,"PointMasses":self.PointMasses,"Rings":self.Rings,"OtherSources":self.OtherSources}
        self.Geodesics={}# maybe use memoization?
        self.num_of_Geods={'unperturbed':0,'perturbed':0}

    # def __str__(self) -> str:
    #     pass
    
    #getter for zs
    
    def Add_Perturbation_Source(self,Model_Or_zs_Or_Eij,*args:float,**kwargs)->None:
        NewSource=0
        if type(Model_Or_zs_Or_Eij) is str:
            if Model_Or_zs_Or_Eij=='point':
                NewSource=PointMass(*args,**kwargs)
                self.PointMasses.append(NewSource)
                # self.Perturbers.append(NewSource)
                if NewSource.name=='Point Mass':NewSource.set_name('P:'+str(len(self.PointMasses)))
            elif Model_Or_zs_Or_Eij=='ring':
                NewSource=AccretionRing(*args,**kwargs)
                self.Rings.append(NewSource)
                # self.Perturbers.append(NewSource)
                if NewSource.name=='Accretion Ring':NewSource.set_name('R:'+str(len(self.Rings)))
        elif isinstance(Model_Or_zs_Or_Eij,type(np.ndarray(5,dtype=np.complex128))):
            NewSource=Perturber(zs=Model_Or_zs_Or_Eij,*kwargs)
            self.OtherSources.append(NewSource)
            # self.Perturbers.append(NewSource)
            if NewSource.name=='Pre-set Perturbation':NewSource.set_name('S:'+str(len(self.OtherSources)))
        elif isinstance(Model_Or_zs_Or_Eij,type(np.ndarray((3,3)))):
            NewSource=Perturber(Eij=Model_Or_zs_Or_Eij,*kwargs)
            self.OtherSources.append(NewSource)
            # self.Perturbers.append(NewSource)
            if NewSource.name=='Pre-set Perturbation':NewSource.set_name('S:'+str(len(self.OtherSources)))
        
        #try:
        if NewSource==0:print(f"{Model_Or_zs_Or_Eij} is not a valid argument. Try 'point','ring', np.array(5,dtype=np.complex128), or np.array((3,3)).")
        else:
            checkname=self.NetPerturbation.name
            self.NetPerturbation+=NewSource #Memory leak issue?
            if checkname=='Flat Background': self.NetPerturbation.set_name(NewSource.name)
            else: self.NetPerturbation.set_name('NetPerturbation')
            #self.Perturbers['NetPerturbation'] deletion is not necessary in Python, only overwriting the reference
            self.Perturbers['NetPerturbation']=self.NetPerturbation
        #except:print(f"{Model_Or_zs_Or_Eij} is not a valid argument. Try 'point','ring', np.array(5,dtype=np.complex128), or np.array((3,3)).")
    
    def Clear_Pertubations(self):
        delattr(self,'NetPerturbation')
        self.NetPerturbation=Perturber(Model='Flat Background',name='Flat Background')
        self.PointMasses=[]
        self.Rings=[]
        self.OtherSources=[]
    
    def Christoffels(self,radius,theta,phi): return TotalGamma(self.a,radius,theta,phi,z_array=self.NetPerturbation.z_array,M=self.M)

    def IntRHS(self):
        if self.NetPerturbation.name=='Flat Background':
            def base_RHS(t, y):
                time,radius,theta,phi=y[0:4]
                u=y[4:8]
                
                udot=-1*np.einsum('ijk,j,k->i',TotalGamma(self.a,radius,theta,phi),u,u)
                udd=-2*np.einsum('ijk,j,k->i',TotalGamma(self.a,radius,theta,phi),u,udot)
                uddd=-2*(np.einsum('ijk,j,k->i',TotalGamma(self.a,radius,theta,phi),udot,udot)+np.einsum('ijk,j,k->i',TotalGamma(self.a,radius,theta,phi),u,udd))
                return np.concatenate((u,udot,udd,uddd))
            return base_RHS
        else:
            def pert_RHS(t, y):
                time,radius,theta,phi=y[0:4]
                u=y[4:8]
                #xdot=u
                udot=-1*np.einsum('ijk,j,k->i',self.Christoffels(radius,theta,phi),u,u)
                udd=-2*np.einsum('ijk,j,k->i',self.Christoffels(radius,theta,phi),u,udot) # dGamma/dtau=0 for a static perturbation on a stationary background
                uddd=-2*(np.einsum('ijk,j,k->i',self.Christoffels(radius,theta,phi),udot,udot)+np.einsum('ijk,j,k->i',self.Christoffels(radius,theta,phi),u,udd))
                return np.concatenate((u,udot,udd,uddd))
            return pert_RHS
        
    def GeodesicConstructor(self,p,e,x,xi_steps=10,chi_steps=10,phi_steps=10):
        if (p,e,x) not in self.Geodesics: # use frequency or resonance as key?
            NewGeodesicSet=Geodesic_Set(self.M,self.a,p,e,x,xi_steps,chi_steps,phi_steps)
            self.Geodesics[(p,e,x)]=NewGeodesicSet
            for key in self.Geodesics[(p,e,x)].Trajectories:
                func=self.IntRHS()
                uududd=func(t=None,y=self.Geodesics[(p,e,x)].Trajectories[key].ICs.flatten())
                # vel0,=uududd[0:4],
                acc0,jerk0=uududd[4:8],uududd[8:12]
                # uddd=,uududd[12:]
                self.Geodesics[(p,e,x)].Trajectories[key].ICs=np.append(self.Geodesics[(p,e,x)].Trajectories[key].ICs,[acc0,jerk0],axis=0)


            #phase_keys=func(params and init pos and vel)
            # if self.NetPerturbation.z_array is None: self.num_of_Geods['unperturbed']+=1
            # else:
            # if self.NetPerturbation.z_array is not None:
            #     try:
            #         NewGeodesic.Perturb(self.NetPerturbation.Model,self.NetPerturbation.name,self.NetPerturbation.z_array)
            #         self.num_of_Geods['perturbed']+=1
            #     except ValueError:
            #         print("Perturb function needs dtypes (str,str,ndarray).")
            #nested memoization into self.Geodesics


        return self.Geodesics[(p,e,x)] 
    

    # def Integrator(self,tau_end:int,)
        

if __name__ == '__main__':

    # from typing import Dict, List, Tuple
    # import collections
    # import re
    import sys

    # USAGE = (f"Usage: {sys.argv[0]} "
    #          "[--help] | [-s <sep>] [first [incr]] last")

    # def seq(operands: List[int], sep: str = "\n") -> str:
    #     first, increment, last = 1, 1, 1
    #     if len(operands) == 1:
    #         last = operands[0]
    #     if len(operands) == 2:
    #         first, last = operands
    #         if first > last:
    #             increment = -1
    #     if len(operands) == 3:
    #         first, increment, last = operands
    #     last = last + 1 if increment > 0 else last - 1
    #     return sep.join(str(i) for i in range(first, last, increment))

    # def parse(args: List[str]) -> Tuple[str, List[int]]:
    #     arguments = collections.deque(args)
    #     separator = "\n"
    #     operands: List[int] = []
    #     while arguments:
    #         arg = arguments.popleft()
    #         if not len(operands):
    #             if arg == "--help":
    #                 print(USAGE)
    #                 sys.exit(0)
    #             if arg in ("-s", "--separator"):
    #                 separator = arguments.popleft() # if arguments else None
    #                 continue
    #         try:
    #             operands.append(int(arg))
    #         except ValueError:
    #             raise SystemExit(USAGE)
    #         if len(operands) > 3:
    #             raise SystemExit(USAGE)

    #     return separator, operands

    # def main() -> None:
    #     sep, operands = parse(sys.argv[1:])
    #     if not operands:
    #         raise SystemExit(USAGE)
    #     print(seq(operands, sep))

    # if __name__ == "__main__":
    #     main()

    opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
    args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]



