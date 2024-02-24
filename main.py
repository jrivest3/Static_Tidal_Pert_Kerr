#main.py
# import os
import numpy as np
import DataHandler
# import SpacetimeController
#import Integrator

a1=np.linspace(0,1,.1)
b1=np.linspace(0,1,.1)
a2=np.linspace(1.1,2,.1)
b2=np.linspace(1.1,2,.1)

# DH=DataHandler

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