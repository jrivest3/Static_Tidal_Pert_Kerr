# testvectorize.py
import numpy as np

X,Y,Z=np.meshgrid(np.arange(0,2,1),np.arange(0,2,1),np.arange(0,2,1))

def myinnerfunc(a,x,y,z,v):return a*v*z+y

def myfunc(x,y,z):
    
    return myinnerfunc(.01,x,y,z,x*np.array([[.1*y,.2+3j],[1,2]]))


vfunc=np.vectorize(myfunc,signature='(),(),()->(a,a)')
vals=vfunc(X,Y,Z)
print(X.shape,vals.shape)
print(vals)