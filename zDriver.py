#zDriver.py
#
# from os import open
import time
from SpacetimeController import *
import pandas as pd
start_time=time.perf_counter()
USAGE=""
description="Print or return a list of lists of z's that can be looped through"
#, begins by constructing a dictionary of Spacetimes"

# import pandas as pd
# import matplotlib.pyplot as plt

# @pd.api.extensions.register_dataframe_accessor("spacetime")
# class STAccessor:
#     def __init__(self, pandas_obj):
#         self._validate(pandas_obj)
#         self._obj=pandas_obj
        
#     @staticmethod
#     def _validate(obj):
#         # verify are columns for date, tmin, tmax, and tavg
#         if not all(col in obj.columns for col in ['DATE','TMAX','TMIN','TAVG']):
#             raise AttributeError("Columns must include 'DATE','TMAX','TMIN', and 'TAVG'")
#         if not all(nm == obj['NAME'][0] for nm in obj['NAME']):
#             raise AttributeError("All values in NAME column must be the same")
  
#     @property
#     def start_date(self):
#         # return the time series start date
#         return pd.to_datetime(self._obj.sort_values('DATE',axis=0)['DATE'].iloc[0]).strftime('%Y-%m-%d')
    
#     @property
#     def end_date(self):
#         # return the time series end date
#         return pd.to_datetime(self._obj.sort_values('DATE',axis=0)['DATE'].iloc[-1]).strftime('%Y-%m-%d')
  
#     @property
#     def station_name(self):
#         # return the station name
#         return self._obj['NAME'][0]
    
#     def plot_temperature(self):
#         fig,ax = plt.subplots()
#         ax.plot(pd.to_datetime(self._obj['DATE']),self._obj['TAVG'],marker='o',markersize=4,markerfacecolor='w',lw=1,markevery=2,label='Average Daily Temperature')
#         ax.fill_between(pd.to_datetime(self._obj['DATE']),self._obj['TMIN'],self._obj['TMAX'],alpha=.4,color='yellow')
#         ax.axhline(1.0, linestyle=':', lw=1)
#         title_str = 'Daily Temperature Range: {start} to {end}\n{station_name}'.format(
#             start=self.start_date,
#             end=self.end_date,
#             station_name=self.station_name)
#         ax.set_title(title_str)
#         ax.xaxis.set_major_formatter(dates.DateFormatter('%m/%d'))
#         ax.set_ylabel('Temperature (F)')
#         ax.set_xlabel('Date')
#         ax.legend()


# from dataclasses import dataclass
# from dataclasses import field
# from dataclasses import InitVar
# # import pandas as pd
# import random 
# import math

# @dataclass
# class DataClass_Modern(object):

#     # Invisible attribute (init-only)
#     attr0:InitVar[int] = 81

#     # Initialized attribute
#     attr1:int   =0
#     attr2:float =0.
#     attr3:str   ='undefined'
#     attr4:list  = field(default_factory=list)

#     # Generated attribute
#     attr5:float = field(init=False)

#     # Generated attribute - read property
#     @property
#     def attr5(self)->float:
#         return math.sqrt(abs(self._attrHidden))

#     # Generated attr - set property (required by dataclasses)
#     @attr5.setter
#     def attr5(self,_):pass # Do nothing, this is a read-only attribute

#     def __post_init__(self,attr0):
#         # Make a copy of the init-only attribute to a local attribute that
#         # is used for the generated attribute (via a property)
#         self._attrHidden = attr0 # This attribute should remain hidden from pandas
#     @classmethod

#     def rand_factory(cls):
#         '''
#         Returns an object of the calling class with randomized initialization attributess
#         '''
#         return cls(
#             attr0=random.randint(-1000,1000),
#             attr1=random.randint(-10,10),
#             attr2=random.random(),
#             attr3=random.choice(
#                 [
#                     'Tool',
#                     'PinkFloyd',
#                     'Soundgarden',
#                     'FaithNoMore',
#                     'aPerfectCircle',
#                     'KingCrimson',
#                     'PearlJam',
#                     'ChildrenOfBodom']),
#             attr4=random.choices(range(100,999),k=3)
#         )

# if __name__ == '__main__':

#     rand_objects = [DataClass_Modern.rand_factory() for _ in range(100)]
#     df = pd.DataFrame(rand_objects)

    # print(df)

# class SpacetimeDataBase:
#     def __init__(self) -> None:
#         self.Spacetimes={}

#     def SpacetimeConstructor(self,a):# *args,**kwargs):
#         if tuple([a,z1,z2,z3,z4,z5]) not in self.Spacetimes:
#             NewST=Spacetime(a)

PertST={}
def SpacetimeConstructor(a,eps,ptype,theta,phi):# *args,**kwargs):
        if tuple([a,eps,ptype,theta,phi]) not in PertST:
            PertST[(a,eps,ptype,theta,phi)]=Spacetime(a)
            PertST[(a,eps,ptype,theta,phi)].Add_Perturbation_Source(ptype,theta,phi,Epsilon=eps)

    # theta_p,phi_p=np.pi/2,0
    # PertST[(theta_p,phi_p)]={}

eps_arr=[1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6] # 10**(-6)
a_arr=[.001,.1,.2,.3,.4,.5,.75,.9,.99]
#The trajectory will be rotated, so we don't need to rotate phi of the source.
thetas=np.linspace(0, np.pi/2, 5,endpoint=False)
[SpacetimeConstructor(a,eps,'point',theta,0) for theta in thetas for eps in eps_arr for a in a_arr]
[SpacetimeConstructor(a,eps,'point',np.pi/2,0) for eps in eps_arr for a in a_arr]
[SpacetimeConstructor(a,eps,'ring',theta,0) for theta in thetas for eps in eps_arr for a in a_arr]

# PertST_keys=[key for key in PertST]
PertST_zs=[PertST[key].NetPerturbation.z_array for key in PertST] 
print('Generated',len(PertST_zs),'sources/spacetimes')
print(time.perf_counter()-start_time)

ST_IDs=pd.DataFrame([key for key in PertST])
ST_ref_dict={
      'model': [PertST[key].NetPerturbation.model for key in PertST],
      'spin': [PertST[key].a for key in PertST],
      'strength': [PertST[key].NetPerturbation.model for key in PertST],
      'polar angle': [PertST[key].NetPerturbation.theta if PertST[key].NetPerturbation.theta else None for key in PertST],
      'azimuthal angle': [PertST[key].NetPerturbation.theta if PertST[key].NetPerturbation.theta else None for key in PertST]
}

# Resonance Controller

# with open('./spacetimes.txt', 'x') as f:
#     f.write(f'{}')
# print('wrote to ./spacetimes.txt')
# Compound/Complex Sources
# compPertST={}
# skipped_list=[]
# for keyi in PertST:
#     for keyj in PertST:
#         if keyj[0]=='ring':continue
#         elif keyi==keyj or keyi+keyj in compPertST or keyj+keyi in compPertST: continue
#         elif keyi[0]==keyj[0] and np.isclose(keyi[1],np.pi/2-keyj[1]) and np.isclose(keyi[2],keyj[2]+np.pi): continue
#         else:
#             new_charzs=PertST[keyi].NetPerturbation.z_array+PertST[keyj].NetPerturbation.z_array
#             new_zs=new_charzs/np.linalg.norm(new_charzs)
#             if any([np.allclose(new_zs,skipped,1e-3,1e-2) for skipped in skipped_list]): continue
#             elif any([np.allclose(new_zs,compPertST[key].NetPerturbation.z_array,1e-3,1e-2) for key in compPertST]): 
#                 skipped_list.append(new_zs)
#                 print(keyi,keyj,'skipped',time.perf_counter()-start_time)
#             else:
#                 compPertST[keyi+keyj]=Spacetime(a)#PertST[keyi].Add_Perturbation_Source(PertST[keyj].NetPerturbation.z_array)
#                 compPertST[keyi+keyj].Add_Perturbation_Source(new_zs) # compPertST[keyi+keyj].NetPerturbation.add(PertST[keyi].NetPerturbation.z_array,PertST[keyj].NetPerturbation.z_array))
#                     # compPertST[keyi+keyj].NetPerturbation.add_and_normalize_zs()

# print(time.perf_counter()-start_time)
# compPertST_zs=[compPertST[key].NetPerturbation.z_array for key in compPertST]#[PertST[keyi].NetPerturbation.z_array+PertST[keyj].NetPerturbation.z_array for keyi in PertST for keyj in PertST])
# print(len(PertST_zs),len(compPertST_zs))
# for keyi in compPertST:
#     rowstring=''
#     for keyj in compPertST:
#         if keyi!=keyj and np.allclose(compPertST[keyj].NetPerturbation.z_array,compPertST[keyi].NetPerturbation.z_array,1e-3,1e-3):
#             rowstring+=str(keyi)+' sim '+str(keyj)+'  '
#     if rowstring!='':print(rowstring)

# compPertST_keys=[key for key in compPertST]

# for i in range(len(PertST_keys)):
#     rowstring=''
#     for j in range(len(compPertST_keys)):
#         if np.allclose(compPertST[compPertST_keys[j]].NetPerturbation.z_array,PertST[PertST_keys[i]].NetPerturbation.z_array,1e-3,1e-2):
#             rowstring+=str(keyi)+' sim '+str(keyj)+'  '
#     if rowstring!='':print(rowstring)

# print(time.perf_counter()-start_time)
# We can control the norm of the net z's by renormalizing after adding sources, but the ratio of contributions from each source will remain; 
# e.g. if a point mass with eps1=1e-6 is added to another with eps2=2e-6, the second still contributes more after dividing the net z-array by the new norm
# for ptype in ['point','ring']:
#     for theta in thetas:
#         for phi in phis:

#             PertST[('combo',PertST[(ptype,theta,phi)].Add_Perturbation_Source('point',theta,phi,Epsilon=eps)
