#Kluge.py
import numpy as np
# from numpy import cos, sin
from scipy.spatial.transform import Rotation

# import pickle

# from time import time

# from matplotlib import pyplot as plt

def constructKlugeWaveformsFromKinematics(kinematics,observedLongitude = np.pi / 4, observedInclination = np.pi / 6,observerDistance=1,mode='spherical'):
    if mode=='spherical': #convert spherical kinematics to cartesian
        # project onto cartesian coordinates in flat space
        # kinematics is assumed to have the shape (4(vars),4(dim),Nsteps)     ####Old Method:(Nsteps,vars,dim) made by np.reshape(np.transpose(sol.y),(len(sol.t),4,4)) with sol.y.shape=(16,Nsteps) as returned by scipy.integrate.solve_ivp
        t, r, θ, ϕ = kinematics[0,0],kinematics[0,1],kinematics[0,2],kinematics[0,3]
        tdot, rdot, θdot, ϕdot = kinematics[1,0],kinematics[1,1],kinematics[1,2],kinematics[1,3]
        tddot, rddot, θddot, ϕddot = kinematics[2,0],kinematics[2,1],kinematics[2,2],kinematics[2,3]
        tdddot, rdddot, θdddot, ϕdddot = kinematics[3,0],kinematics[3,1],kinematics[3,2],kinematics[3,3]
        
        x = r * np.sin(θ) * np.cos(ϕ) 
        y = r * np.sin(θ) * np.sin(ϕ) 
        z = r * np.cos(θ)  

        # compute various derivatives of x_{p}^{μ} wrt τ
        dx = (np.cos(ϕ) * (np.sin(θ) * rdot + np.cos(θ) * r * θdot) - r * np.sin(θ) * np.sin(ϕ) * ϕdot) / tdot  # Eq. D.1
        d2x = (np.cos(ϕ) * np.sin(θ) * rddot + 2 * rdot * (np.cos(θ) * np.cos(ϕ) * θdot - np.sin(θ) * np.sin(ϕ) * ϕdot) + r * (np.cos(θ) * (-2 * np.sin(ϕ) * θdot * ϕdot + np.cos(ϕ) * θddot) - np.sin(θ) * (np.cos(ϕ) * (θdot**2 + ϕdot**2) + np.sin(ϕ) * ϕddot))) / (tdot**2)  # Eq. D.2
        d3x = (3 * (θddot * np.cos(θ) - np.sin(θ) * θdot**2) * (rdot * np.cos(ϕ) - r * ϕdot * np.sin(ϕ)) - 3 * θdot * np.cos(θ) * (np.cos(ϕ) * (r * ϕdot**2 - rddot) + np.sin(ϕ) * (2 * rdot * ϕdot + r * ϕddot)) + np.sin(θ) * (np.cos(ϕ) * (rdddot - 3 * ϕdot * (rdot * ϕdot + r * ϕddot)) + np.sin(ϕ) * (r * (ϕdot**3 - ϕdddot) - 3 * (rddot * ϕdot + rdot * ϕddot))) + r * np.cos(ϕ) * ((θdddot - θdot**3) * np.cos(θ) - 3 * θdot * θddot * np.sin(θ))) / (tdot**3)  # Eq. D.3

        dy = (np.sin(ϕ) * (np.sin(θ) * rdot + np.cos(θ) * r * θdot) + np.cos(ϕ) * r * np.sin(θ) * ϕdot) / tdot  # Eq. D.4
        d2y = (np.sin(ϕ) * (np.sin(θ) * (rddot - r * (θdot**2 + ϕdot**2)) + np.cos(θ) * r * θddot) + 2 * rdot * (np.cos(θ) * np.sin(ϕ) * θdot + np.cos(ϕ) * np.sin(θ) * ϕdot) + np.cos(ϕ) * r * (2 * np.cos(θ) * θdot * ϕdot + np.sin(θ) * ϕddot)) / (tdot**2)  # Eq. D.5
        d3y = (3 * (rdot * np.sin(θ) + r * θdot * np.cos(θ)) * (ϕddot * np.cos(ϕ) - np.sin(ϕ) * ϕdot**2) + 3 * ϕdot * np.cos(ϕ) * (rddot * np.sin(θ) + 2 * θdot * rdot * np.cos(θ) + r * (θddot * np.cos(θ) - np.sin(θ) * θdot**2)) + np.sin(ϕ) * (rdddot * np.sin(θ) + 3 * rddot * θdot * np.cos(θ) + 3 * rdot * (θddot * np.cos(θ) - np.sin(θ) * θdot**2) + r * ((θdddot - θdot**3) * np.cos(θ) - 3 * θdot * θddot * np.sin(θ))) + r * np.sin(θ) * ((ϕdddot - ϕdot**3) * np.cos(ϕ) - 3 * ϕdot * ϕddot * np.sin(ϕ))) / (tdot**3)  # Eq. D.6

        dz = (np.cos(θ) * rdot - r * np.sin(θ) * θdot) / tdot  # Eq. D.7
        d2z = (np.cos(θ) * (-r * θdot**2 + rddot) - np.sin(θ) * (2 * rdot * θdot + r * θddot)) / (tdot**2)  # Eq. D.8
        d3z = (rdddot * np.cos(θ) - 3 * rddot * θdot * np.sin(θ) - 3 * rdot * (np.cos(θ) * θdot**2 + θddot * np.sin(θ)) + r * ((θdot**3 - θdddot) * np.sin(θ) - 3 * θdot * θddot * np.cos(θ))) / (tdot**3)  # Eq. D.9

        # xμ = [x, y, z]
        # vμ = [dx, dy, dz]
        # aμ = [d2x, d2y, d2z]
        # jerkμ = [d3x, d3y, d3z]
        
        Cartesian_Kinematics= np.transpose(np.array([[x, y, z],[dx, dy, dz],[d2x, d2y, d2z],[d3x, d3y, d3z]]),axes=[2,0,1]) # change shape from (4,3,nsteps) to (nsteps,4,3) 

    elif mode=='3d': 
        # print(kinematics)
        Cartesian_Kinematics=kinematics # must have shape (Nsteps,4-[pos,vel,acc,jerk],3-[x,y,z])

    rhat = Rotation.from_euler('xz',[-observedInclination, -np.pi / 2 + observedLongitude])
    rhat = np.transpose(rhat.as_matrix())
    ex,ey,ez = rhat[0],rhat[1],rhat[2]

    z_projector = np.identity(3) - np.outer(ez, ez)

    #The following is taken from Alejandro's notes (sent via email).
    #Only the l=2 contribution is constructed
    mu = 1 #object mass
    R = observerDistance # Doesn't matter, because we are comparing phase-mismatch, rather than amplitude.

    h_plusList = []
    h_crossList = []

    for x, v, a, j in Cartesian_Kinematics:
        #Construct the l=2 I tensor
        I2 = mu * (np.outer(a, x) + 2 * np.outer(v, v) + np.outer(x, a))
        I2STF = I2 - np.trace(I2)/3 * np.identity(3)

        #construct the l=2 J tensor
        #cross product is used instead of levi-civita
        va = np.cross(v, a)
        xa = np.cross(x, a)
        xj = np.cross(x, j)
        xv = np.cross(x, v)

        J2 = np.outer(x, va) + 2 * np.outer(v, xa)\
            + np.outer(x, xj) + np.outer(a, xv)
        # J2 *= mu
        J2STF = (J2+np.transpose(J2))/2 - np.trace(J2)/3 * np.identity(3)

        #J tensor contribution to h_ij, forced symmetric here
        #(numpy docs say the cross product is performed on the last axis)
        J_contribution = 4/3/R * np.cross(J2STF, ez)
        J_contribution += np.transpose(J_contribution)

        #h_ij and the traceless transverse form
        h_ij = 2/R * I2STF + J_contribution
        hTT_ij = z_projector.dot(h_ij.dot(z_projector))\
            - 1/2 * z_projector * np.trace(z_projector.dot(h_ij))

        #The final polarizations
        h_plus = (hTT_ij.dot(ex).dot(ex) - hTT_ij.dot(ey).dot(ey))/2
        h_cross = (hTT_ij.dot(ex).dot(ey) + hTT_ij.dot(ey).dot(ex))/2

        h_plusList.append(h_plus)
        h_crossList.append(h_cross)

    return h_plusList, h_crossList

# def convertKinematicsDataToStrain(inputFile, outputFile,
# 	observedLongitude = np.pi / 4, observedInclination = np.pi / 6):

# 	timeList = []
# 	kinematics = []
# 	with open(inputFile, 'rb') as handle:
# 		timeList, minoTimeList, kinematics = pickle.load(handle)

# 	h_plusList, h_crossList =\
# 		constructKlugeWaveformsFromKinematics(kinematics,
# 			observedLongitude, observedInclination)

# 	with open(outputFile, 'wb') as handle:
# 		pickle.dump((timeList, h_plusList, h_crossList), handle,
# 			protocol=pickle.HIGHEST_PROTOCOL)
if __name__=='__main__':
    nsteps=10
    test_shpr_traj=np.zeros((nsteps,4,4))
    test_shpr_traj[:,:,0]=1 # t
    test_shpr_traj[:,:,1]=np.pi/2 # r
    test_shpr_traj[:,:,2]=np.pi/3 # theta, phi's are zero
    print(constructKlugeWaveformsFromKinematics(test_shpr_traj,mode='spherical'))#