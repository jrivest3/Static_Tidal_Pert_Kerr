#SpacetimeController

# import time

#import numpy as np
# import spherical
# import quaternionic

from GeodesicController import *
from CristhValuesArrays import TotalGamma


from scipy.integrate import solve_ivp as ODE
# from scipy.special import ellipk,ellipe

sqrtof2=np.sqrt(2)
sqrtof3=np.sqrt(3)
sqrtof5=np.sqrt(5)
sqrtof6=sqrtof2*sqrtof3
sqrtpiover5=np.sqrt(np.pi/5)

# start=time.perf_counter()
# restart=start

class Perturber:
    def __init__(self,zs=None,Eij=None,name=None,Model='Pre-set Perturbation') -> None:
        self.Model=Model #Perhaps assume a flat backgroud if zs are not given
        self.name='Pre-set Perturbation' if name is None else name
        self.z_array=zs #if zs is not None else np.zeros(5)
        self.Eij= Eij #if Eij is not None else np.zeros((3,3))

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
            self.z_array=np.array([np.conj(z2),-z1,z0,z1,z2])
            # # Harmonic coordinates
            # z0,z1,z2=-2*(Eij[0,0]+Eij[1,1]),-2*(Eij[0,2]-complex(0,1)*Eij[1,2]),-2*(Eij[0,0]-Eij[1,1])+complex(0,4)*Eij[0,1]
            # self.z_array=np.array([np.conj(z2),np.conj(z1),z0,z1,z2])
        elif self.Eij is None: print("No Eij's have been given for this perturber.")
        else: print(f"These z's are already {self.z_array}")

    def add(self,a,b):
        new_zs=a.z_array+b.z_array if a.z_array is not None else b.z_array
        new_Eij=a.Eij+b.Eij if a.Eij is not None and b.Eij is not None else None
        return Perturber(zs=new_zs,Eij=new_Eij,name='Net Perturbation',Model='Net Perturbation')
    

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
                z0,z1,z2=2*sqrtof6*sqrtpiover5* ( 1 - 3 * CosTHp*CosTHp ), 12*sqrtpiover5 *expPHp *CosTHp*SinTHp, -6*sqrtpiover5 /expPHp/expPHp *SinTHp*SinTHp#Cartesion (BL?)
                Norm = 4*sqrtof6*sqrtpiover5
                self.z_characteristic= np.array([np.conj(z2),-z1,z0,z1,z2])/Norm # Norm = 4*sqrt(6*pi/5)~7.8 for all theta and phi.
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
        self.Long_of_Asc_Node=Long_of_Asc_Node # angle from phi=0
        psi=Long_of_Asc_Node+np.pi/2 # value of phi at max inclination(zm)


        if self.z_array is None:
            if self.Eij is None:
                expPsi=np.exp(1.j*psi)
                if zm*zm==1:raise ValueError("Polar rings are not yet implemented. ")
                z0,z1,z2=   sqrtof6*sqrtpiover5*( 2 + -3 * zm*zm ), \
                            8*sqrtpiover5/( np.pi )/( zm ) * \
                                ( ( -1 + 2 * ( zm*zm ) ) * ellipe( ( zm*zm ) ) - ( -1 + ( zm*zm ) ) * ellipk( ( zm*zm ) ) ) \
                                *expPsi, \
                            3/2 * sqrtpiover5 * ( zm*zm ) /expPsi/expPsi
                zchar=np.array([np.conj(z2),-z1,z0,z1,z2])
                Norm=np.linalg.norm(zchar) # ~3.84 near pole, peaks ~4.01 at theta~0.46276 (just over 7*pi/48), 2*sqrtof6*sqrtpiover5 ~3.88 at the equator
                self.z_characteristic= np.array([np.conj(z2),-z1,z0,z1,z2])/Norm
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

    #getter for zs
    
    def Add_Perturbation_Source(self,Model_Or_zs_Or_Eij,*args:float,**kwargs)->None:
        NewSource=0
        if type(Model_Or_zs_Or_Eij) is str:
            if Model_Or_zs_Or_Eij=='point':
                NewSource=PointMass(*args,**kwargs)
                self.PointMasses.append(NewSource)
                # self.Perturbers.append(NewSource)
                if NewSource.name=='Point Mass':NewSource.set_name(len(self.PointMasses))
            elif Model_Or_zs_Or_Eij=='ring':
                NewSource=AccretionRing(*args,**kwargs)
                self.Rings.append(NewSource)
                # self.Perturbers.append(NewSource)
                if NewSource.name=='Accretion Ring':NewSource.set_name(len(self.Rings))
        elif isinstance(Model_Or_zs_Or_Eij,type(np.ndarray(5,dtype=np.complex128))):
            NewSource=Perturber(zs=Model_Or_zs_Or_Eij,*kwargs)
            self.OtherSources.append(NewSource)
            # self.Perturbers.append(NewSource)
            if NewSource.name=='Pre-set Perturbation':NewSource.set_name(len(self.OtherSources))
        elif isinstance(Model_Or_zs_Or_Eij,type(np.ndarray((3,3)))):
            NewSource=Perturber(Eij=Model_Or_zs_Or_Eij,*kwargs)
            self.OtherSources.append(NewSource)
            # self.Perturbers.append(NewSource)
            if NewSource.name=='Pre-set Perturbation':NewSource.set_name(len(self.OtherSources))
        
        #try:
        if NewSource==0:print(f"{Model_Or_zs_Or_Eij} is not a valid argument. Try 'point','ring', np.array(5,dtype=np.complex128), or np.array((3,3)).")
        else:
            if self.NetPerturbation.name=='Flat Background': self.NetPerturbation.set_name(NewSource.name)
            self.NetPerturbation=self.NetPerturbation.add(self.NetPerturbation,NewSource) #Memory leak issue?
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
                #xdot=u
                udot=-1*np.einsum('ijk,j,k->i',TotalGamma(self.a,radius,theta,phi),u,u)
                return np.concatenate((u,udot))
            return base_RHS
        else:
            def pert_RHS(t, y):
                time,radius,theta,phi=y[0:4]
                u=y[4:8]
                #xdot=u
                udot=-1*np.einsum('ijk,j,k->i',self.Christoffels(radius,theta,phi),u,u)
                return np.concatenate((u,udot))
            return pert_RHS
        
    def GeodesicConstructor(self,p,e,x,**kwargs):
        if (self.a,p,e,x) not in self.Geodesics:
            NewGeodesic=Geodesic(self.M,self.a,p,e,x,**kwargs)
            self.Geodesics[NewGeodesic.params]=NewGeodesic
            #phase_keys=func(params and init pos and vel)
            if self.NetPerturbation.z_array is not None:
                try:
                    NewGeodesic.Perturb(self.NetPerturbation.Model,self.NetPerturbation.name,self.NetPerturbation.z_array)

                except ValueError:
                    print("Perturb function needs dtypes (str,str,ndarray).")
            #nested memoization into self.Geodesics
        return self.Geodesics[(self.a,p,e,x)]
    # def Integrator(self,)
        


