#Create_DataBase.py
#
import time
from collections import deque
import sys#,requests
from typing import List,Tuple
# import parser
import h5py
# from GeodesicController import Geodesic_Set
# from SpacetimeController import Perturber
from SpacetimeController import *
# import zDriver
start_time=time.perf_counter()
# maybe have a default Nsteps here
MY_USAGE = (f"Usage: {sys.argv[0]} "
         "[--help] | [-i <input_file.h5>] | [<spin:float>_<model:str>_<pert strength:float>_<polar angle position:float>_<azimuthal angle position:float=0>] <Inegration_Time_or_Nsteps:int>")
HELP_MESSAGE='''Help not Finished

Creates a Spacetime Database HDF5 file and fills in the trajectories and strains of a set of geodesics for the given list of parameters. 
For now, only one Group with p,e,x is generated, so either an input file or a '_'-separated list of params can be given.
Note: the list of parameters must be separated by '_'s.
If no file is given, but the file that would be created already exists, the existing on with be appended.
If given a .h5 file containing a database of spacetimes, it will check if a Group already exists for a Geodesics with parameters exists. If not and a new Group with be added to it.
'''
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

def parse(args: List[str]) -> Tuple[str, List[str]]:
    arguments = deque(args)
    filename=''
    #separator = " "
    #operands: List[int] = []
    operands= []
    while arguments:
        arg = arguments.popleft()
        if not len(operands):
            if arg == "--help":
                print(MY_USAGE)
                print(HELP_MESSAGE)
                sys.exit(0)
            if arg in ("-i", "--input_file"):
                filename = arguments.popleft() if arguments else ''#None
                break#continue
        # try:
        operands.append(arg) #validation will be worked out later
        # except ValueError:
        #     raise SystemExit(MY_USAGE)
        # if len(operands) > 3:
        #     raise SystemExit(MY_USAGE)
    
    try: 
        int(operands[-1]) #default Nsteps value?
    except TypeError:
        raise SystemExit('no time integer given\n'+MY_USAGE)
    
    return filename, operands


def toRealNum(s:str):
    '''If the string represents a real number, then return a float with the proper sign. Else, give back the string.'''

    try : 
        res=float(s) #positive float
        return res
    except : 
        if s.startswith('-'):
            try :
                new_s=s[1:] 
                res=float(new_s) #positive float
                return -res
            except : 
                return s
        return s
    
def tuple_to_filename_string(input_tuple)->str:
    '''
    Function to convert a tuple into a string suitable for a filename.
 
    Parameters:
    - input_tuple: tuple
        The input tuple that needs to be converted into a filename string.
 
    Returns:
    - str:
        A string suitable for a filename based on the input tuple. '-'s are replaced with 'n's if a word begines with '-'
 
    Raises:
    - TypeError:
        Raises an error if the input is not a tuple.
    '''
 
    # Check if the input is a tuple
    if not isinstance(input_tuple, tuple):
        raise TypeError("Input should be a tuple.")
 
    # Convert the tuple elements to strings
    tuple_str = [str(elem) for elem in input_tuple] # also handle negative values #.replace(' -',' n') 

    # Join the tuple elements with underscores to create the filename string
    filename_str = "_".join(tuple_str)
 
    return filename_str

def SpacetimeConstructor(a,ptype,eps,theta,phi,p,e,x,this_PertST,PertG):# *args,**kwargs):
    # if tuple([a,ptype,eps,theta,phi]) not in this_PertST:
    this_PertST=Spacetime(a)
    this_PertST.Add_Perturbation_Source(ptype,theta,phi,Epsilon=float(eps))
    if tuple([a,p,e,x]) not in PertG:
        PertG[(a,p,e,x)]={}
        PertG[(a,p,e,x)]=this_PertST.GeodesicConstructor(p,e,x)
    return this_PertST,PertG

def write_db_groups(ST_f,this_PertST,PertG,GsetsList:List,err_message=True):
    if not this_PertST: raise ValueError(f'Spacetime {this_PertST=} not created yet')
    if ST_f: #check that it's already open
        if not GsetsList:
            #chunks?
            set_index=0
            db_time=time.perf_counter()
            for key in PertG:
                a,p,e,x=key
                GsetsList.append(ST_f.create_group(f'Set_{a}_{p}_{e}_{x}',track_order=True))
                for attr, value in this_PertST.__dict__.items(): # __dict__ or vars(), not dir() # import collections for dictionary unpacking?
                    if value is None:value='None'
                    if isinstance(value,(Perturber,Geodesic_Set)):
                        for inner_attr, inner_value in this_PertST[attr].__dict__.items():
                            if inner_value is None:inner_value='None'
                            print(f"{attr}: Inner Attribute: {inner_attr}, Value: {inner_value}")
                            GsetsList[set_index].attrs.create(f'{attr}.{inner_attr}',inner_value)
                    else:
                        print(f"Attribute: {attr}, Value: {value}")
                        GsetsList[set_index].attrs.create(f'{attr}',value)
                set_index+=1
            print('group creation time:',time.perf_counter()-db_time)
            return GsetsList
        elif err_message:raise FileExistsError(f'write_db_groups Already Ran: {len(GsetsList)=}')
    print('write_db_groups Failed: ST_f is closed.')


def fill_database(ST_f,t_len,this_PertST,PertG,GsetsList,PertTraj,PertStrain):#,group_names:List[str]
    if not this_PertST: raise ValueError(f'Spacetime {this_PertST=} not created yet')
    if not ST_f: raise FileNotFoundError(f'{ST_f=} File is not open.')
    # try:
    Gsets:List=write_db_groups(ST_f,this_PertST,PertG,GsetsList)#,err_message=False
    # except FileExistsError:
    # try:
    set_index=0
    integ_time=time.perf_counter()
    for key in PertG:
        # a,p,e,x=key
        tau_array=np.linspace(0,t_len,t_len)
        Gsets[set_index].create_dataset('tau_array',data=tau_array)
        for psi_0_key in PertG[key].Trajectories: # groups and datasets
            Traj=this_PertST.run_Trajectory(PertG[key].Trajectories[psi_0_key],t_eval=tau_array)
            PertTraj[psi_0_key]=Gsets[set_index].create_dataset('Traj_'+tuple_to_filename_string(psi_0_key),data=Traj)
            Strain=this_PertST.calc_Strain(PertG[key].Trajectories[psi_0_key])
            PertStrain[psi_0_key]=Gsets[set_index].create_dataset('Strain_'+tuple_to_filename_string(psi_0_key),data=Strain)
        print('time for set',set_index,':',time.perf_counter()-integ_time)
        set_index+=1
    print('to fill one h5:',time.perf_counter()-integ_time)
    return this_PertST,PertG,Gsets,PertTraj,PertStrain
    # except: print(f'fill_database Failed: {PertG=}')


def main():#get_File_and_Group_names():
    this_PertST=None;PertG={};GsetsList=[];PertTraj={};PertStrain={}
    h5name, operands = parse(sys.argv[1:])
    [p,e,x]= [10.0,0,1] #These and the z's should be reported for reproducibility. Don't let x=0. [float(x) for x in args[1:4]] if len(args) else

    # opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
    # args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]# # fortunately, spin is not allowed to be negative

    # if ('-from_h5_file' in opts) or ('-i_h5' in opts)

    # params= operands[0].replace('_n','_-') #params are negative if they start with 'n' # used if I ever take p,e,x as args
    # ID_index=args.pop()
    # ST_id can be created either from the name (args) or the ST_IDs_df

    if h5name=='':
        if not operands:raise SystemExit('no operands\n'+MY_USAGE)
    
        # params= operands[0]#"".join([c if c.isalnum() else "_" for c in str(ST_id)])
    
        paramList=[toRealNum(arg) for arg in operands[0].split('_')]
        if len(paramList)<4 or len(paramList)>5: 
            print(f'need 4 or 5 parameters. {len(paramList)} were given: {paramList=}')
            raise SystemExit(MY_USAGE)
        elif len(paramList)==4:
            paramList.append(0)
        # assert len(paramList)==5
        # ST_id=tuple(paramList)
        h5name='ST'+operands[0]+'.hdf5'
        try:    
            ST_f=h5py.File(h5name,'x')
        except FileExistsError:
            print("A file with these parameters already exists.")
            ST_f=h5py.File(h5name,'r+') # need some strong validation
        
        try: # should validate all
            a:float=paramList[0] 
            ptype:str=paramList[1]
            eps:float=paramList[2];theta:float=paramList[3];phi:float=paramList[4]
        except:
            a,ptype,eps,theta,phi=paramList
            raise SystemExit(f'Tried:\n{[a,ptype,eps,theta,phi,p,e,x]=}, \n'+MY_USAGE)
        print('before constructor:',time.perf_counter()-start_time)
        this_PertST,PertG=SpacetimeConstructor(a,ptype,eps,theta,phi,p,e,x,this_PertST,PertG)
        print('before fill_database:',time.perf_counter()-start_time)
        fill_database(ST_f,int(operands[1]),this_PertST,PertG,GsetsList,PertTraj,PertStrain)
        print('Done:',time.perf_counter()-start_time)
        # try:
        #     Result:tuple=fill_database(ST_f,int(operands[1]),this_PertST,PertG,GsetsList,PertTraj,PertStrain)
        #     this_PertST,PertG,GsetsList,PertTraj,PertStrain=Result
        # except TypeError:
        #     fill_database(ST_f,int(operands[1]),this_PertST,PertG,GsetsList,PertTraj,PertStrain)
        #     print('still thinks it might be NoneType')
        # spin,model,strength,theta_p,phi_p=
        # ### Resonance Controller goes here?

        

    else:
        raise ValueError('input_file option not yet functional.')
        # if not h5name.endswith('.h5'):
        #     print('Currently only works with h5 files.')
        #     raise SystemExit(MY_USAGE)
        # # Not currently validating input
        
        # ST_f=h5py.File(h5name,'a') # need some strong validation
        # group_names= list(ST_f.keys())
        # try:
        #     a:float= spin
        #     SpacetimeConstructor(a,ptype,eps,theta,phi,p,e,x)
            
        #     # spin,model,strength,theta_p,phi_p=ST_id # validate eventually
        # except:raise SystemExit(MY_USAGE)
    
    
    # ST_f.close()

    # return h5name, group_names

    # Gset=ST_f[f'Set_{a}_{p}_{e}_{x}']

    
    # args[0] could be tau_end or Nsteps, or args[4] since the vector graphs dont require it
    # options: one source, n-sources
    #OR I could have one function to calculate z's then pipe that output to the input
    # maybe args should each be arg_name=value






if __name__ == "__main__":
    main()