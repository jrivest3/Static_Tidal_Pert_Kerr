#zDriver.py

import os,sys
# from timebudget import timebudget
import time
from numpy import pi,linspace
# from SpacetimeController import *
import pandas as pd
# import h5py
from typing import List,Tuple
from Create_DataBase import tuple_to_filename_string

start_time=time.perf_counter()

USAGE=(f"Usage: {sys.argv[0]} "
         "[--help] | [-s | --submit] <batch_size> <Nsteps> <num_epsilons> <num_spins>")
description="Print or return a list of lists of z's that can be looped through"
#, begins by constructing a dictionary of Spacetimes"

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]# # fortunately, spin is not allowed to be negative

if "--help" in opts:
    print(USAGE)
    print(description)
    sys.exit(0)

 # maybe activate my conda env for convenience
submit_flag=False
if ('-s' in opts) or ('--submit' in opts): submit_flag=True

max_eps_arr=[1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6]
max_a_arr=[.001,.1,.5,.9,.99]#,.2,.3,.4,.75
batch_size=int(args[0]) if len(args)>0 else 3
Int_Nsteps=int(args[1]) if len(args)>1 else 10
eps_num=int(args[2]) if len(args)>2 else 1#len(max_eps_arr)
a_num=int(args[3]) if len(args)>3 else 1#len(max_a_arr)
try:
    assert eps_num<=len(max_eps_arr)
    assert a_num<=len(max_a_arr)
except:
    print(f'{max_eps_arr=}, so eps_num<={len(max_eps_arr)}\n{max_a_arr=}, so a_num<={len(max_a_arr)}')
    raise SystemExit(USAGE)


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

class SpacetimeDataBase:
    eps_arr=max_eps_arr[:eps_num]
    a_arr=max_a_arr[:a_num]
    PertST_keys=[]
    
    def __init__(self) -> None:
        # self.Spacetimes={}

        # def SpacetimeConstructor(a,ptype,eps,theta,phi):# *args,**kwargs):
        #         if tuple([a,ptype,eps,theta,phi]) not in PertST:
        #             PertST[(a,ptype,eps,theta,phi)]=Spacetime(a)
        #             PertST[(a,ptype,eps,theta,phi)].Add_Perturbation_Source(ptype,theta,phi,Epsilon=eps)

        # theta_p,phi_p=np.pi/2,0
        # PertST[(theta_p,phi_p)]={}
        
        #The trajectory will be rotated, so we don't need to rotate phi of the source.
        self.thetas=linspace(0, pi/2, 5,endpoint=False)
        [self.PertST_keys.append((a,'point',eps,theta,0)) for theta in self.thetas for eps in self.eps_arr for a in self.a_arr]
        [self.PertST_keys.append((a,'point',eps,pi/2,0)) for eps in self.eps_arr for a in self.a_arr]
        [self.PertST_keys.append((a,'ring',eps,theta,0)) for theta in self.thetas for eps in self.eps_arr for a in self.a_arr]
        
        self.batches=[self.PertST_keys[i:i+batch_size] for i in range(0,len(self.PertST_keys),batch_size)]

        # PertST_keys=[key for key in PertST]
        # PertST_zs=[PertST[key].NetPerturbation.z_array for key in PertST] 
        print('Generated',len(self.PertST_keys),'sources/spacetimes, in',len(self.batches),'batches and ',time.perf_counter()-start_time,'sec')

        self.ST_IDs_df=pd.DataFrame([key for key in self.PertST_keys],columns=['spin','model','strength','polar angle','azimuthal angle'])
        # ST_ref_dict={
        #       'spin': [PertST[key].a for key in PertST],
        #       'model': [PertST[key].NetPerturbation.model for key in PertST],
        #       'strength': [PertST[key].NetPerturbation.model for key in PertST],
        #       'polar angle': [PertST[key].NetPerturbation.theta if PertST[key].NetPerturbation.theta else None for key in PertST],
        #       'azimuthal angle': [PertST[key].NetPerturbation.phi if PertST[key].NetPerturbation.phi else None for key in PertST]
        # }
        # Qthesame_df=pd.DataFrame(ST_ref_dict)
        # print(ST_IDs_df)
        # print(Qthesame_df)
        # Qthesame_df.to_csv('Spacetimes', sep='\t') # chunks?
        self.ST_IDs_df.to_csv('ST_IDs.csv')

        
    def Spacetime_pbs_script(self,batch_num:int, walltime: str = "01:00:00", conda_env_name:str="r2840_env"):#, job_name: str, python_script_path: str, args_batch: List, output_dir: str, nodes: int = 1, ppn: int = 1
        '''
        Function to generate a PBS script for submitting a Python function to a local cluster.

        Parameters:

        - batch_num: int
            Which batch in the self.batches list.

        - walltime: str (default: "01:00:00")
            Maximum time the job is allowed to run in the format "hh:mm:ss".

        Returns:
        - str:
            The generated PBS script as a string.

        '''

        # - job_name: str
        #     Name of the job to be submitted.
        # - python_script_path: str
        #     Path to the Python script that needs to be executed.
        # - args_batch: List
        #     List of arguments to be passed to the Python script.
        # - output_dir: str
        #     Directory where the output files will be stored.
        # - nodes: int (default: 1)
        #     Number of nodes to be used for the job.
        # - ppn: int (default: 1)
        # #     Number of processors per node.
        # Raises:
        # - ValueError:
        #     If the nodes or ppn values are not positive integers.
        # One CPU per Spacetime
        #for fparams in fnames:
        job_name = 'batch_'+str(batch_num)
        python_script_path = "Create_DataBase.py"
        # args_batch=(batch_num,batch,Int_Nsteps)
        output_dir = "./outputs"
        # nodes = 2 #not used
        # ppn = 4 #not used
        # walltime = "01:00:00"


        file_path = python_script_path #'path/to/your/file_or_directory'
        if not os.path.exists(file_path):
        #     # print(f"{file_path} exists.")
            
        #     if os.path.isfile(file_path):
        #         print(f"{file_path} is a regular file.")
        #     elif os.path.isdir(file_path):
        #         print(f"{file_path} is a directory.")
        # else:
            print(f"{file_path} does not exist.")
            raise FileNotFoundError


        
        # Validating nodes and ppn values
        # if nodes <= 0 or ppn <= 0:
        #     raise ValueError("Nodes and ppn values should be positive integers.")
        #PBS -l nodes={nodes}:ppn={ppn}
        # Constructing the PBS script content# print(pbs_script)
        #Making a PHS file for each simulation
        pbsfile = open('./pbsfiles/pbsfile'+str(batch_num)+'.pbs', 'x')
        # pbsfile.write(pbs_script)
        pbsfile.write('#!/bin/bash\n')
        # pbsfile.write('#PBS -N '+ 'AHjobfile'+ str(i) + '.py\n')
        # pbsfile.write('#PBS -j oe\n')
        # pbsfile.write('#PBS -l select=1:ncpus=1:mem=1gb -l place=free\n \n')
        # pbsfile.write('module load python\n \n')
        # # the main python file that calculates coeffs is passed along with the path to the simulation folder
        # pbsfile.write('python /ddn/home8/r2571/run_AH/AH_7_main.py {0}'.format(path) )
        #pbs_script = f'''#!/bin/bash
        pbsfile.write('#PBS -N %s\n'%(job_name))
        pbsfile.write('#PBS -l ncpus=1\n')
        pbsfile.write('#PBS -l mem=1gb\n')
        pbsfile.write(f'#PBS -l walltime={walltime}\n')
        pbsfile.write(f'#PBS -l cput={walltime}\n')
        pbsfile.write('#PBS -j oe\n')
        pbsfile.write(f'#PBS -o {output_dir}/{job_name}.out\n')
        pbsfile.write(f'#PBS -e {output_dir}/{job_name}.err\n')
        pbsfile.write('#PBS -m a\n')
        pbsfile.write('#PBS -M ljrivest@go.olemiss.edu\n\n')

        pbsfile.write('cd $PBS_O_WORKDIR\n\n')

        pbsfile.write('# Load any necessary modules\n')
        pbsfile.write('module load python\n')
        pbsfile.write(f'source activate {conda_env_name}\n\n')

        pbsfile.write('# Execute the Python script\n')
        for id_params in self.batches[batch_num-1]:
            string_Params=tuple_to_filename_string(id_params)
            pbsfile.write(f'python3 ./{python_script_path} '+string_Params+' '+str(Int_Nsteps)+'\n')
        
        fname=pbsfile.name
        pbsfile.close()
        print('fname:',fname)
        return fname


# Unstructured Database
# create csv of key info and h5 filenames
# Each run, or set of geodesics, has an h5 file containing a metadata file, whose contents are copied to a row of the csv file, and directories containing the phase trajectories and strain data arrays
# Later functionality will be added to query this database for post-processing. For now I will process some example data.



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




# if __name__ == '__main__':
    # Example of using the generate_pbs_script function:
# ID_index = 0

# PertG={}
# [p,e,x]= [10.0,0,1] # Resonance controller here?
My_Spacetimes=SpacetimeDataBase()

batch_num=0
for batch in My_Spacetimes.batches:
# while PertST_keys:
#     batches.append((batch_num,[]))
#     for num in range(num_per_batch):
#         batches[batch_num][1].append(PertST_keys.)
    batch_num+=1
    print(f'batch {batch_num}: {batch}')
    # for id in batch:

    # >>> from collections import namedtuple

    # >>> Person = namedtuple("Person", "name age height weight")
    # >>> jane = Person("Jane", 26, 1.75, 67)
    # >>> for field, value in zip(jane._fields, jane):
    # ...     print(field, "->", value)
    # ...
    # name -> Jane
    # age -> 26
    # height -> 1.75
    # weight -> 67

    # fparams=tuple_to_filename_string(key)

        # for keyname in ST_ref_dict: # or ST_IDs
        #     names[i]+=f' {ST_ref_dict[keyname][i]}'#{keyname}=
    
    # ST_f=h5py.File('ST'+fparams+'.h5','x')
    # PertG[key]={}
    # PertG[key][(a,p,e,x)]=PertST[key].GeodesicConstructor(p,e,x)
    # #chunks?
    # Gset=ST_f.create_group(f'Set_{a}_{p}_{e}_{x}',track_order=True)
    # for attr, value in PertST[key].__dict__.items(): # __dict__ or vars(), not dir() # import collections for dictionary unpacking?
    #     print(f"Attribute: {attr}, Value: {value}")
    #     Gset.attrs.create(f'{attr}',value)
    #     if isinstance(attr,(Perturber,Geodesic_Set)):
    #         for inner_attr, inner_value in PertST[key][attr].__dict__.items():
    #             print(f"\tInner Attribute: {inner_attr}, Value: {inner_value}")
    #             Gset.attrs.create(f'{inner_attr}',inner_value)
    # ST_f.close()

        # "".join([c if c.isalnum() else "_" for c in name])
    pbs_filename = My_Spacetimes.Spacetime_pbs_script(batch_num=batch_num) #
    
    
    #Submit the file in queue for running on cluster
    if submit_flag:os.system('qsub '+pbs_filename) #pbsfile{0}.pbs'.format(str(i))
    
    # ID_index+=1

    # requests.get from import requests ?

# #import multiprocessing as mproc
# #import numpy as np
# import ray
# from timebudget import timebudget
# from DataHandler import *
# # import SpacetimeController
# #import Integrator


# # Should each remote call have it's own main loop, or should the main loop contain all the remote calls?
# # I want the DataHandler to store or take dictionaries and create arrays and calculations for results and plots
# # nested for loops can make a dictionary with a long tuple as a key value, or it could make a nested directory
# # I could also create a function or super dictionary that returns desired key values and then creates an entry for that key value if it doesn't already exist
# # the simplest thing may be something like a pandas pivot table

# a1=np.linspace(0,1,3)
# b1=np.linspace(0,1,10)
# a2=np.linspace(1.1,2,10)
# b2=np.linspace(1.1,2,10)

# # DH=DataHandler
# def rev_ls_dec(ls_func):
#     def reverser(*args):
#         print("in decor")
#         ret_val=ls_func(*args)
#         ret_val.reverse()
#         print('back in decor',ret_val)
#         return ret_val
#     return reverser

# def lister(arr:np.ndarray)->list:
#     ls=[]
#     for el in arr:
#         ls.append(el)
#     return ls
# a1list=lister(a1)
# a1list.reverse()

# @rev_ls_dec
# def revlister(arr:np.ndarray)->list:
#     ls=[]
#     print('in func')
#     for el in arr:
#         ls.append(el)
#     return ls

# print(a1,a1list,revlister(a1))

# # cpu1
# for a in a1:
#     for b in b1:
#         pass#f(a,b) give data to DH

# #cpu2
# for a in a2:
#     for b in b2:
#         pass#f(a,b) give data to DH

# #cpu3
# # store/manage DH data, run DH functions