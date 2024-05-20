#main.py
# import os
#import multiprocessing as mproc
#import numpy as np
import ray
from timebudget import timebudget
from DataHandler import *
# import SpacetimeController
#import Integrator


# Should each remote call have it's own main loop, or should the main loop contain all the remote calls?
# I want the DataHandler to store or take dictionaries and create arrays and calculations for results and plots
# nested for loops can make a dictionary with a long tuple as a key value, or it could make a nested directory
# I could also create a function or super dictionary that returns desired key values and then creates an entry for that key value if it doesn't already exist
# the simplest thing may be something like a pandas pivot table

a1=np.linspace(0,1,3)
b1=np.linspace(0,1,10)
a2=np.linspace(1.1,2,10)
b2=np.linspace(1.1,2,10)

# DH=DataHandler
def rev_ls_dec(ls_func):
    def reverser(*args):
        print("in decor")
        ret_val=ls_func(*args)
        ret_val.reverse()
        print('back in decor',ret_val)
        return ret_val
    return reverser

def lister(arr:np.ndarray)->list:
    ls=[]
    for el in arr:
        ls.append(el)
    return ls
a1list=lister(a1)
a1list.reverse()

@rev_ls_dec
def revlister(arr:np.ndarray)->list:
    ls=[]
    print('in func')
    for el in arr:
        ls.append(el)
    return ls

print(a1,a1list,revlister(a1))

# cpu1
for a in a1:
    for b in b1:
        pass#f(a,b) give data to DH

#cpu2
for a in a2:
    for b in b2:
        pass#f(a,b) give data to DH

#cpu3
# store/manage DH data, run DH functions