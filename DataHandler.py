#DataHandler
from SpacetimeController import *
import pandas as pd
import matplotlib.pyplot as plt

import time
start=time.perf_counter()
restart=start
M=1
# class DataHandler:
#     def __init__(self) -> None:
#         self.Spacetimes={}

#     def SpacetimeConstructor(self,a):# *args,**kwargs):
#         if tuple([a,z1,z2,z3,z4,z5]) not in self.Spacetimes:
#             NewST=Spacetime(a)

Base_Convergence_graphs=False
Pert_Stability_graphs=True
Pert_Convergence_graphs=False

a=0.0001
p,e,x=10.0,0,1 #These and the z's should be reported for reproducibility. Don't let x=0.
KerrST=Spacetime(a)
KerrG=KerrST.GeodesicConstructor(p,e,x)

tau_end,Nsteps=3*10**3,5 # tau_end= 10000*ceil(T_orbit_for_r_min)
times=np.linspace(0,tau_end,tau_end)#*10)# tau_end/floor(T_orbit_for_r_min)*20 or 10000orbits*20 or more

exps=np.arange(8,14.) #exponents for the integration tolerances
eps_ar=[0.,10**(-9)] #,10**(-12),10**(-6),10**(-5)] 1e-5 diverges, 1e-6 diverges for off equator perturbers

if Base_Convergence_graphs:
    #Test Stability and Convergence when unperturbed 
    del_r_max=np.zeros(6) # tracker of max variance of base orbit from r0 for when e=0, contains the max vals of del_r for each exp of tolerance
    t_delr_max=np.zeros(6) # the times for del_r_max
    del_theta_max=np.zeros(6) # tracker of max variance of base orbit from th0 for when x=+/-1, contains the max vals of del_theta for each exp of tolerance
    t_delth_max=np.zeros(6) # the times for del_theta_max
    base_t_arrays=[]
    base_r_arrays=[]
    base_theta_arrays=[]
    base_phi_arrays=[]
    # del_r_arrays=[]
    # del_theta_arrays=[]
    for i in range(exps.size):
        restart=time.perf_counter()
        rtol, atol = 10.**(-exps[i]), 10.**(-exps[i])
        base_sol=ODE(KerrST.IntRHS(),[0,tau_end],KerrG.ICs.flatten(),t_eval=times,rtol=rtol,atol=atol)# Pre-set times makes it take significantly longer
        print(f"Unperturbed with tol={rtol}: time to integrate {tau_end*10} steps: {time.perf_counter()-restart}")
        # base_t_arrays.append(base_sol.t)
        base_r_arrays.append(base_sol.y[1])
        base_theta_arrays.append(base_sol.y[2])
        base_phi_arrays.append(base_sol.y[3])
        # t_array_len= base_sol.t.size
        # print("tau_end=",base_sol.t[-1],"base_sol runtime=",time.perf_counter()-start,"sol.t length=",base_sol.t.shape)
        # print("base size=",base_sol.y.shape,base_sol.nfev)
        # if e==0:
        #     r0=p
        #     del_r=[0.]
        #     for j in range(1,base_sol.t.size):
        #         r=base_sol.y[1][j]

        #         if abs(del_r[-1])<abs(r-r0): 
        #             del_r_max[i]=r-r0
        #             t_delr_max[i]=base_sol.t[j]
                    
        #         del_r.append(r-r0)
        #     del_r_arrays.append(del_r)

        # if x*x==1:
        #     th0=np.pi/2
        #     del_theta=[0.]
        #     for j in range(1,base_sol.t.size):
        #         theta=base_sol.y[2][j]
                
        #         if abs(del_theta[-1])<abs(theta-th0): 
        #             del_theta_max[i]=theta-th0
        #             t_delth_max[i]=base_sol.t[j]
                
        #         del_theta.append(theta-th0)
        #     del_theta_arrays.append(del_theta)
        # # print("Del_r=",del_r)

    #fig_base_sols: 8x6 base_sol.y[j][i] vs base_t_arrays[i]
    fig_base_sols, axs_base=plt.subplots(3,1)
    # for c in range(4): # columns (not t),r,th,ph
    for exp in range(exps.size): # tol exp
        axs_base[0].plot(times,base_r_arrays[exp],label="1e-%d" % exps[exp])
        axs_base[1].plot(times,base_theta_arrays[exp],label="1e-%d" % exps[exp])
        axs_base[2].plot(times,base_phi_arrays[exp],label="1e-%d" % exps[exp])
            # axs_base[1,c].plot(base_t_arrays[exp],del_r_arrays[exp],label="1e-%f" % exps[exp])# row for velocity
    axs_base[0].set_xlabel('proper time')
    axs_base[0].set_ylabel('radius')
    axs_base[0].legend()
    axs_base[1].set_xlabel('proper time')
    axs_base[1].set_ylabel('$\\theta$')
    axs_base[1].legend()
    axs_base[2].set_xlabel('proper time')
    axs_base[2].set_ylabel('$\phi$')
    axs_base[2].legend()
    fig_base_sols.suptitle(f'Unperturbed Convergence: a,p,e,x = {a,p,e,x}')
    fig_base_sols.set_figheight(15)
    fig_base_sols.set_figwidth(18)
    plt.savefig(f'Unpert_params_{a}_{p}_{e}_{x}.png')
    plt.show()
            
    #fig_del_r for exp in exps: 6 plots of del_r_arrays[i] vs base_t_arrays[i]

    #fig_del_th for exp in exps: 6 plots of del_theta_arrays[i] vs base_t_arrays[i]
    fig_del, (ax_del_r,ax_del_th,ax_del_ph)=plt.subplots(3,1)
    for i in range(exps.size-1):
        ax_del_r.plot(times,base_r_arrays[i]-base_r_arrays[-1],label="1e-%d" % exps[i])#del_r_arrays[i]
        ax_del_th.plot(times,base_theta_arrays[i]-base_theta_arrays[-1],label="1e-%d" % exps[i])#del_theta_arrays[i]
        ax_del_ph.plot(times,base_phi_arrays[i]-base_phi_arrays[-1],label="1e-%d" % exps[i])#del_theta_arrays[i]
    ax_del_r.set_xlabel('proper time')
    ax_del_r.set_ylabel('$\delta$r')
    ax_del_th.set_xlabel('proper time')
    ax_del_th.set_ylabel('$\delta\\theta$')
    ax_del_ph.set_xlabel('proper time')
    ax_del_ph.set_ylabel('$\delta\phi$')
    ax_del_r.legend()
    ax_del_th.legend()
    ax_del_ph.legend()
    fig_del.suptitle(f"Unperturbed Deviations from tol=1e-{exps[-1]}: a,p,e,x = {a,p,e,x}")
    fig_del.set_figheight(15)
    fig_del.set_figwidth(18)
    plt.savefig(f'UnpertDeviations_params_{a}_{p}_{e}_{x}.png')
    plt.show()

    # #fig_max_dels 2 scatter plots: del_r_max and t_delr_max vs tol and del_theta_max and t_delth_max vs tol
    # fig_max_dels, (ax_dr_max,ax_dth_max)=plt.subplots(1,2)
    # sizes=10*exps
    # cb1=ax_dr_max.scatter(t_delr_max,del_r_max,s=sizes,c=exps,cmap='RdBu_r')
    # ax_dr_max.set_title('Max $\delta$r')
    # fig_max_dels.colorbar(cb1,ax=ax_dr_max)
    # if x*x==1: 
    #     cb2=ax_dth_max.scatter(t_delth_max,del_theta_max,s=sizes,c=exps,cmap='RdBu_r')
    #     ax_dth_max.set_title('Max $\delta\\theta$')
    #     fig_max_dels.colorbar(cb2,ax=ax_dth_max)
    # fig_max_dels.suptitle("Tolerance Exp = marker size")
    # # fig_max_dels.set_figheight(15)
    # # fig_max_dels.set_figwidth(18)
    # plt.savefig('UnpertMaxDevs.png')
    # plt.show()

PertST={}
PertG={}
for (theta_p,phi_p) in [(np.pi-np.arccos(1/np.sqrt(3)),-np.pi/3)]:#,(np.pi/2,np.pi/4),(np.pi/4,0),(0,0)]:
    PertST[(theta_p,phi_p)]={}
    PertG[(theta_p,phi_p)]={}
    for i in range(len(eps_ar)):    
        PertST[(theta_p,phi_p)][eps_ar[i]]=Spacetime(a)
        PertST[(theta_p,phi_p)][eps_ar[i]].Add_Perturbation_Source('point',theta_p,phi_p,Epsilon=eps_ar[i])# -10*eps_ar[i]*np.array([100*(1+0.j),0,0,0,100*(1-0.j)]))
        PertG[(theta_p,phi_p)][eps_ar[i]]=PertST[(theta_p,phi_p)][eps_ar[i]].GeodesicConstructor(p,e,x)
        # print(PertST[(theta_p,phi_p)][eps_ar[i]].NetPerturbation.z_array,PertG[(theta_p,phi_p)][eps_ar[i]].zs)
    if Pert_Stability_graphs:
    # Test Stability for different Perturbation strengths
        pert_t_arrays=[]
        pert_r_arrays=[]
        pert_theta_arrays=[]
        pert_phi_arrays=[]
            # pert_ys=[]
        t_div=np.zeros(len(eps_ar)) # time what integration of perturbed orbit ended noting when the integration ends early do to divergence
        y_div=np.zeros((len(eps_ar),8)) # pert_sol.y at t_div for each eps
        for i in range(len(eps_ar)):
            rtol, atol = 10.**(-12), 10.**(-12)

            restart=time.perf_counter()
            pert_sol=ODE(PertST[(theta_p,phi_p)][eps_ar[i]].IntRHS(),[0,tau_end],PertG[(theta_p,phi_p)][eps_ar[i]].ICs.flatten(),rtol=rtol,atol=atol)#,t_eval=times,args=PertG[(theta_p,phi_p)][eps].zs
            print(f"Perturbed with eps={eps_ar[i]}: time to integrate to {tau_end}: {time.perf_counter()-restart}")
            PertG[(theta_p,phi_p)][eps_ar[i]].Trajectory=pert_sol
            pert_t_arrays.append(pert_sol.t)
            pert_r_arrays.append(pert_sol.y[1])
            pert_theta_arrays.append(pert_sol.y[2])
            pert_phi_arrays.append(pert_sol.y[3])

            # pert_ys.append(pert_sol.y)
            pt_array_len= pert_sol.t.size
            if pert_sol.t[-1]<tau_end:
                t_div[i]=pert_sol.t[-1]
                y_div[i]=[val[-1] for val in pert_sol.y]
            else:
                t_div[i]=np.nan
                y_div[i]=[np.nan for _ in pert_sol.y]
        print((start-time.perf_counter())/60,"min")
        #fig_pert_sols: 8x6 pert_sol.y[j][i] vs pert_t_arrays[i]
        fig_pert_sols, axs_pert=plt.subplots(3,1)
        # for c in range(4): # columns (not t),r,th,ph
        for i in range(len(eps_ar)):
            axs_pert[0].plot(pert_t_arrays[i],pert_r_arrays[i],label=f"$\epsilon$={eps_ar[i]}")
            axs_pert[1].plot(pert_t_arrays[i],pert_theta_arrays[i],label=f"$\epsilon$={eps_ar[i]}")
            axs_pert[2].plot(pert_t_arrays[i],pert_phi_arrays[i],label=f"$\epsilon$={eps_ar[i]}")
                # axs_pert[1,c].plot(pert_t_arrays[exp],del_r_arrays[exp],label="1e-%f" % exps[exp])# row for velocity
        axs_pert[0].set_xlabel('proper time')
        axs_pert[0].set_ylabel('radius')
        axs_pert[0].legend()
        axs_pert[1].set_xlabel('proper time')
        axs_pert[1].set_ylabel('$\\theta$')
        axs_pert[1].legend()
        axs_pert[2].set_xlabel('proper time')
        axs_pert[2].set_ylabel('$\phi$')
        axs_pert[2].legend()
        fig_pert_sols.set_figheight(15)
        fig_pert_sols.set_figwidth(18)
        # fig_pert_sols.suptitle(f"Point-Like Perturber ($\\theta_p,\phi_p$)={theta_p,phi_p}\n tol=1e-12 a,p,e,x={a,p,e,x}") 
        # plt.savefig(f'PLpert_zp{np.cos(theta_p)}_phip{phi_p}_aftercorrections.png')
        plt.show()
        # #fig_pert_sols: 8x6 pert_sol.y[j][i] vs pert_t_arrays[i]
        # fig_pert_sols, axs_pert=plt.subplots(1,3)
        # # for c in range(4): # columns (not t),r,th,ph
        # for i in range(len(eps_ar)):
        #     axs_pert[0].plot(times,pert_r_arrays[i]-base_r_arrays[4],label=f"$\epsilon$={eps_ar[i]}")
        #     axs_pert[1].plot(times,pert_theta_arrays[i]-base_theta_arrays[4],label=f"$\epsilon$={eps_ar[i]}")
        #     axs_pert[2].plot(times,pert_phi_arrays[i]-base_phi_arrays[4],label=f"$\epsilon$={eps_ar[i]}")
        #         # axs_pert[1,c].plot(pert_t_arrays[exp],del_r_arrays[exp],label="1e-%f" % exps[exp])# row for velocity
        # axs_pert[0].set_xlabel('proper time')
        # axs_pert[0].set_ylabel('radius')
        # axs_pert[0].legend()
        # axs_pert[1].set_xlabel('proper time')
        # axs_pert[1].set_ylabel('$\\theta$')
        # axs_pert[1].legend()
        # axs_pert[2].set_xlabel('proper time')
        # axs_pert[2].set_ylabel('$\phi$')
        # axs_pert[2].legend()
        # fig_pert_sols.suptitle(f"Perturbed minus Unperturbed")
        # plt.savefig(f'PLpert_zp{np.cos(theta_p)}minusBase.png')
        # plt.show()

        # #fig_t_div plot of t_div vs tol_exp
        # fig_t_div, ax_t_div=plt.subplots()
        # fig_t_div.suptitle('Proper Time to Divergence')
        # ax_t_div.stem(eps_ar,t_div,use_line_collection=True)
        # ax_t_div.set_xlabel('$\epsilon$')
        # ax_t_div.set_ylabel('proper time')
        # # fig_t_div.set_figheight(10)
        # # fig_t_div.set_figwidth(14)
        # plt.savefig(f'PLpert_zp{np.cos(theta_p)}DivergenceTimes.png')
        # plt.show()

        # #fig_pert_div: 8 subplots y_div vs t_div scatter plot labeled by exp
        # fig_y_div, axs_y_div=plt.subplots(1,2)
        # # sizes=10*
        # cbrdiv=axs_y_div[0].scatter(t_div,[val[1] for val in y_div],c=eps_ar,label='r')
        # # axs_y_div[0].scatter(t_div,[val[5] for val in y_div],marker='d',c=exps,label='dr/d$\\tau$')
        # axs_y_div[0].set_xlabel('time of divergence')
        # axs_y_div[0].legend()
        # fig_y_div.colorbar(cbrdiv,ax=axs_y_div[0],extend='both')
        # cbthdiv=axs_y_div[1].scatter(t_div,[val[2] for val in y_div],c=eps_ar,label='$\\theta$')
        # # axs_y_div[1].scatter(t_div,[val[6] for val in y_div],marker='d',s=sizes,c=exps,label='d$\\theta$/d$\\tau$')
        # axs_y_div[1].set_xlabel('time of divergence')
        # axs_y_div[1].legend()
        # fig_y_div.colorbar(cbthdiv,ax=axs_y_div[1],extend='both')
        # plt.savefig(f'PLpert_zp{np.cos(theta_p)}ValuesAtDivergence.png')
        # plt.show()

        # 
        # eps=10**(-5) #eps_ar[3]
        # # strong perturbations will cause accelerations too great for scipy's solver to integrate
        # # 1e-5 diverges around t=6500 for r=12.8  
        # theta_p,phi_p=np.pi/2,0
        # PertST[(theta_p,phi_p)][eps]=Spacetime(a)
        # PertST[(theta_p,phi_p)][eps].Add_Perturbation_Source('point',theta_p,phi_p,Epsilon=eps)
        # PertG[(theta_p,phi_p)][eps]=PertST[(theta_p,phi_p)][eps].GeodesicConstructor(p,e,x)
        # print(PertST[(theta_p,phi_p)][eps].NetPerturbation.z_array,PertG[(theta_p,phi_p)][eps].zs)  

    if Pert_Convergence_graphs:
        # Test Convergence of Perturbed Trajectories
        # exps=np.arange(8,14.) #exponents for the integration tolerances
        eps=10**(-9)
        pert_t_arrays=[]
        pert_r_arrays=[]
        pert_theta_arrays=[]
        pert_phi_arrays=[]
            # pert_ys=[]
        for i in range(exps.size):
            rtol, atol = 10.**(-exps[i]), 10.**(-exps[i])

            restart=time.perf_counter()
            pert_sol=ODE(PertST[(theta_p,phi_p)][eps].IntRHS(),[0,tau_end],PertG[(theta_p,phi_p)][eps].ICs.flatten(),t_eval=times,rtol=rtol,atol=atol)#,args=PertG[(theta_p,phi_p)][eps].zs
            print(f"Perturbed with tol={rtol}: time to integrate {tau_end*10} steps: {time.perf_counter()-restart}")
            pert_t_arrays.append(pert_sol.t)
            pert_r_arrays.append(pert_sol.y[1])
            pert_theta_arrays.append(pert_sol.y[2])
            pert_phi_arrays.append(pert_sol.y[3])
        print((start-time.perf_counter())/60,"min")
        #fig_pert_sols: 8x6 pert_sol.y[j][i] vs pert_t_arrays[i]
        fig_pert_sols, axs_pert=plt.subplots(3,1)
        # for c in range(4): # columns (not t),r,th,ph
        for exp in range(exps.size):# tol exp
            axs_pert[0].plot(pert_t_arrays[exp],pert_r_arrays[exp],label=f"tol={exps[exp]}")# f"$\epsilon$={eps_ar[exp]}")
            axs_pert[1].plot(pert_t_arrays[exp],pert_theta_arrays[exp],label=f"tol={exps[exp]}")# f"$\epsilon$={eps_ar[exp]}")
            axs_pert[2].plot(pert_t_arrays[exp],pert_phi_arrays[exp],label=f"tol={exps[exp]}")# f"$\epsilon$={eps_ar[exp]}")
                # axs_pert[1,c].plot(pert_t_arrays[exp],del_r_arrays[exp],label="1e-%f" % exps[exp])# row for velocity
        axs_pert[0].set_xlabel('proper time')
        axs_pert[0].set_ylabel('radius')
        axs_pert[0].legend()
        axs_pert[1].set_xlabel('proper time')
        axs_pert[1].set_ylabel('$\\theta$')
        axs_pert[1].legend()
        axs_pert[2].set_xlabel('proper time')
        axs_pert[2].set_ylabel('$\phi$')
        axs_pert[2].legend()
        fig_pert_sols.set_figheight(15)
        fig_pert_sols.set_figwidth(18)
        fig_pert_sols.suptitle(f"Point-Like Perturber ($\\theta_p,\phi_p$)={theta_p,phi_p}\n $\epsilon$={eps} a,p,e,x={a,p,e,x}") 
        plt.savefig(f'PLpert_zp{np.cos(theta_p)}_Convergence_for_eps{eps}.png')
        plt.show()
            # #fig_pert_sols: 8x6 pert_sol.y[j][i] vs pert_t_arrays[i]
            # fig_pert_sols, axs_pert=plt.subplots(1,3)
            # # for c in range(4): # columns (not t),r,th,ph
            # for exp in range(exps.size): # tol exp
            #     axs_pert[0].plot(times,pert_r_arrays[exp]-base_r_arrays[exp],label="1e-%d" % exps[exp])
            #     axs_pert[1].plot(times,pert_theta_arrays[exp]-base_theta_arrays[exp],label="1e-%d" % exps[exp])
            #     axs_pert[2].plot(times,pert_phi_arrays[exp]-base_phi_arrays[exp],label="1e-%d" % exps[exp])
            #         # axs_pert[1,c].plot(pert_t_arrays[exp],del_r_arrays[exp],label="1e-%f" % exps[exp])# row for velocity
            # axs_pert[0].set_xlabel('proper time')
            # axs_pert[0].set_ylabel('radius')
            # axs_pert[0].legend()
            # axs_pert[1].set_xlabel('proper time')
            # axs_pert[1].set_ylabel('$\\theta$')
            # axs_pert[1].legend()
            # axs_pert[2].set_xlabel('proper time')
            # axs_pert[2].set_ylabel('$\phi$')
            # axs_pert[2].legend()
            # fig_pert_sols.suptitle("Perturbed Deltas for eps=%f"%eps)
            # plt.show()

        # Deltas from lowest tolerance
        fig_pert_sols, axs_pert=plt.subplots(3,1)
        # for c in range(4): # columns (not t),r,th,ph
        for exp in range(exps.size-1): # tol exp
            axs_pert[0].plot(pert_t_arrays[exp],pert_r_arrays[-1]-pert_r_arrays[exp],label="tol=1e-%d" % exps[exp])
            axs_pert[1].plot(pert_t_arrays[exp],pert_theta_arrays[-1]-pert_theta_arrays[exp],label="tol=1e-%d" % exps[exp])
            axs_pert[2].plot(pert_t_arrays[exp],pert_phi_arrays[-1]-pert_phi_arrays[exp],label="tol=1e-%d" % exps[exp])
                # axs_pert[1,c].plot(pert_t_arrays[exp],del_r_arrays[exp],label="1e-%f" % exps[exp])# row for velocity
        axs_pert[0].set_xlabel('proper time')
        axs_pert[0].set_ylabel('radius')
        axs_pert[0].legend()
        axs_pert[1].set_xlabel('proper time')
        axs_pert[1].set_ylabel('$\\theta$')
        axs_pert[1].legend()
        axs_pert[2].set_xlabel('proper time')
        axs_pert[2].set_ylabel('$\phi$')
        axs_pert[2].legend()
        fig_pert_sols.set_figheight(15)
        fig_pert_sols.set_figwidth(18)
        fig_pert_sols.suptitle(f"Point-Like Perturber ($\\theta_p,\phi_p$)={theta_p,phi_p}\n  $\Delta$(tol=1e-13) $\epsilon$={eps} a,p,e,x={a,p,e,x}") 
        plt.savefig(f'PL_del13_zp{np.cos(theta_p)}_Convergence_for_eps{eps}.png')
        plt.show()

print((start-time.perf_counter())/60,"min")
# print("tau_end=",pert_sol.t[-1],"pert_sol runtime=",time.perf_counter()-restart,"sol.t length=",pert_sol.t.shape)
# print("pert size=",pert_sol.y.shape)
# print(pert_sol.nfev)
# if pert_sol.t[-1]<tau_end: print("Integration ended early because the acceleration exceeded the solver's capacity. This probably means the perturbation destabilized the orbit.")


# from math import floor
# print("base_sol.phi: {}".format([base_sol.y[3][floor(n)] for n in range(0,t_array_len,floor(t_array_len/Nsteps) )]))
# print("pert_sol.phi: {}".format([pert_sol.y[3][floor(n)] for n in range(0,pt_array_len,floor(pt_array_len/Nsteps) )]))
# print("base_sol.r: {}".format([base_sol.y[1][floor(n)] for n in range(0,t_array_len,floor(t_array_len/Nsteps) )]))
# print("pert_sol.r: {}".format([pert_sol.y[1][floor(n)] for n in range(0,pt_array_len,floor(pt_array_len/Nsteps) )]))
