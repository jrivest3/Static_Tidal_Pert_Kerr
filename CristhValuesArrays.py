import time


import numpy as np
import spherical
import quaternionic
from scipy.integrate import solve_ivp as ODE

# Relavent sYlm's
ell=2
ell_max=ell
ell_min=ell
def swsh_array(l, theta, phi):
    rotor = quaternionic.array.from_spherical_coordinates(theta, phi)
    sign_mat = [-1 if (-neg_s %2) else +1
                   for m in range(-l, l+1)
                   for neg_s in range(-l, l+1)
                ]
    sign_mat = np.array(sign_mat)
    W = spherical.Wigner(l, l)
    D_mat = W.D(rotor)
    return sign_mat*np.sqrt((2*l+1.)/(4.*np.pi)) * D_mat

sqrtof2=np.sqrt(2)
sqrtof3=np.sqrt(3)
sqrtof5=np.sqrt(5)
sqrtof6=sqrtof2*sqrtof3

start=time.perf_counter()
restart=start

M=1

def TotalGamma(a,radius,theta,phi,z_array=np.array([0,0,0,0,0])):
    '''
    Claculate the Christoffel symbols for Kerr spacetime (-+++) and an arbitrary, static (no time-depenence) tidal (l=2) perturbation.
    The tidal perturbation is defined by (5) real parameters z_m for m in range(-l,l). 
    Default z_array is zero, which returns only background Christoffel symbols.
    Returns a 3-d numpy array in 3+1 dimensions
    '''
    CosTH=np.cos(theta)
    SinTH=np.sin(theta)

    rho=(-1/(radius - 1.j*a*CosTH))
    rhobar=np.conj(rho) # (-1/(radius + 1.j*a*CosTH))

    PsiCDe2=M*rho*rho*rho
    PsiCDe2bar=M*rhobar*rhobar*rhobar

    Delta= a*a - 2*M*radius + radius*radius
    Deltap= 2*(radius-M)
    Deltapp= 2
    Sigma= 1/rho/rhobar

    # Non-zero Spin Coefficients (in terms of rho)
    tau= -rho*rhobar*(1.j*a*SinTH)/sqrtof2
    taubar= -tau
    pi= taubar*rho/rhobar # is np.conj faster?
    pibar= tau*rhobar/rho
    beta= (-rhobar/(2*sqrtof2)*CosTH/SinTH)
    # betabar= beta*rho/rhobar
    # eps=0
    # epsbar=0
    mu= rho*rho*rhobar*Delta/2
    mubar= mu*rhobar/rho
    gamma= mu + rho*rhobar*Deltap/4 # Computation time may be shortened by a simplification that removes gamma and gammabar
    gammabar= gamma-mu+mubar
    # alpha= pi-betabar
    # alphabar= pibar-beta


#def Gamma(radius,theta,phi): # Calculated by BHPToolkit, which uses(-+++). (+---) gets overall minus sign? maybe?
    Chr=np.zeros((4,4,4),np.complex128)
    
    Chr[0,0,1]=-1/32 * M /mu /rhobar * ( 2 * rho * \
        rhobar + ( ( rhobar*rhobar ) + ( rho*rho ) * ( 1 + 4 * ( a*a ) * ( \
        rhobar*rhobar ) ) ) ) * ( ( rho + ( rhobar + 2 * a * rho * rhobar ) ) \
        * ( -rhobar + rho * ( -1 + 2 * a * rhobar ) ) + 8 * ( taubar*taubar ) )
    Chr[0,0,2]=( a*a ) * M * rho * rhobar * ( rho + rhobar ) * CosTH * SinTH

    Chr[0,1,0]=Chr[0,0,1]

    Chr[0,1,3]=-1/16 /a * M /mu /( rho*rho ) \
        * ( 1/( rhobar*rhobar*rhobar ) ) * ( taubar*taubar ) * ( 8 * ( \
        a*a*a*a ) * ( rho*rho*rho*rho ) * ( rhobar*rhobar*rhobar*rhobar ) + ( \
        -6 * ( a*a ) * ( rho*rho ) * ( rhobar*rhobar ) * ( ( rho + rhobar )*( \
        rho + rhobar ) ) + ( -3 * ( ( rho + rhobar )*( rho + rhobar )*( rho + \
        rhobar )*( rho + rhobar ) ) + 4 * ( a + radius ) * rho * rhobar * ( rho + \
        ( rhobar + 2 * a * rho * rhobar ) ) * ( ( a*a ) * ( rho*rho ) * ( \
        rhobar*rhobar ) + 4 * ( taubar*taubar ) ) ) ) )
    Chr[0,2,0]=Chr[0,0,2]

    Chr[0,2,3]=-( a*a*a ) * M * rho * rhobar * ( rho + rhobar ) * CosTH * ( SinTH*SinTH*SinTH )

    Chr[0,3,1]=Chr[0,1,3]
    Chr[0,3,2]=Chr[0,2,3]

    Chr[1,0,0]=0.5 * M * mu /rho * ( 2 * rho * rhobar + ( ( \
        rhobar*rhobar ) + ( ( rho*rho ) * ( 1 -4 * ( a*a ) * ( rhobar*rhobar ) ) -8 * ( taubar*taubar ) ) ) )

    Chr[1,0,3]=2 /a * mu * ( 1/( rho*rho*rho ) ) * ( 1/( \
        rhobar*rhobar*rhobar ) ) * ( ( pi*pi ) + ( taubar*taubar ) ) * PsiCDe2bar

    Chr[1,1,1]=0.5 * ( -( 1/( mu ) ) * ( 2 * gammabar + mu ) * rho + rhobar )
    Chr[1,1,2]=1/sqrtof2 * ( ( 1/( rhobar ) ) - ( 1/( rho ) ) ) * tau

    Chr[1,2,1]=Chr[1,1,2]
    Chr[1,2,2]=0.5 * Delta * ( rho + rhobar )
 
    Chr[1,3,0]=Chr[1,0,3]

    Chr[1,3,3]=-2/( a*a ) * mu * ( Sigma**5/rho ) \
        * ( taubar*taubar ) * ( ( 2 * ( gammabar -mubar ) * ( ( rho \
        -rhobar )*( rho -rhobar ) ) -2 * ( ( a*a ) -M * radius \
        ) * ( rho*rho ) * ( rhobar*rhobar ) * ( rho + rhobar ) ) * ( \
        taubar*taubar ) + ( Delta + 2 * M * radius ) * ( rho*rho*rho*rho ) * ( \
        rhobar*rhobar*rhobar*rhobar ) * ( rho + rhobar ) * ( ( radius*radius ) + ( a*a ) * np.cos( 2 * theta ) ) )
    Chr[2,0,0]=( a*a ) * M * ( rho*rho ) * ( rhobar*rhobar ) * ( rho + rhobar ) * CosTH * SinTH

    Chr[2,0,3]=- a * M * ( ( a*a ) + ( radius*radius ) ) * ( rho*rho ) * ( rhobar*rhobar ) * ( rho + rhobar ) * CosTH * SinTH

    Chr[2,1,1]=0.5 /sqrtof2 / mu   * rho * ( rho -rhobar ) * taubar
    Chr[2,1,2]=-0.5 * ( rho + rhobar )

    Chr[2,2,1]=Chr[2,1,2]
    Chr[2,2,2]=1/sqrtof2 * ( ( 1/( rhobar ) ) - ( 1/( rho ) ) ) * tau

    Chr[2,3,0]=Chr[2,0,3]

    Chr[2,3,3]=- 1/sqrtof2 /( a*a ) /( rho*rho*rho )  * ( rho - rhobar ) /( rhobar*rhobar*rhobar )  * \
        ( -2* mu * rhobar + M * ( ( ( a*a ) + ( radius*radius ) )*( ( a*a ) + ( \
        radius*radius ) ) ) * ( rho*rho*rho ) * ( rhobar*rhobar*rhobar ) * ( rho + rhobar ) ) * tau

    Chr[3,0,1]=1/8 * a * M /mu * ( rho*rho ) * rhobar * ( 2 * rho * \
        rhobar + ( ( rhobar*rhobar ) + ( ( rho*rho ) * ( 1 -4 * ( a*a ) * ( rhobar*rhobar ) ) -8 * ( taubar*taubar ) ) ) )
    Chr[3,0,2]=-2 * sqrtof2 * a * M * beta * rho * ( rho + rhobar )

    Chr[3,1,0]=Chr[3,0,1]

    Chr[3,1,3]=1/32 /mu /rhobar * ( -( 2 * M -radius ) * \
        ( ( rho + rhobar )*( rho + rhobar )*( rho + rhobar )*( rho + rhobar ) \
        ) + ( -2 * ( a*a ) * rho * rhobar * ( ( rho + rhobar )*( rho + rhobar \
        ) ) * ( rho + ( rhobar + 3 * M * rho * rhobar ) ) + ( ( a*a*a*a ) * ( \
        rho*rho*rho ) * ( rhobar*rhobar*rhobar ) * ( 2 * M * rho * rhobar + \
        -3 * ( rho + rhobar ) ) + ( -4 * ( a*a ) * ( Delta + radius * ( M + radius ) ) \
        * ( rho*rho*rho ) * ( rhobar*rhobar*rhobar ) * ( rho + rhobar ) * \
        np.cos( 2 * theta ) + 4 * ( a*a*a*a ) * ( gammabar -mubar ) * ( \
        rho*rho*rho ) * ( rhobar*rhobar*rhobar ) * np.cos( 4 * theta ) ) ) ) )
    Chr[3,2,0]=Chr[3,3,1]=Chr[3,0,2]

    Chr[3,2,3]=-2 * sqrtof2 * beta /rhobar - a*a * M * rho * rhobar * ( rho + rhobar ) * CosTH * SinTH

    Chr[3,3,1]=Chr[3,1,3]
    Chr[3,3,2]=Chr[3,2,3]

#def dGamma(zs):#,radius,theta,phi):
    if np.linalg.norm(z_array)!=0:
        # Kinnersley Tetrad in Outgoing? Kerr-Newman Coordinates
        lKN=np.array([0,1,0,0])
        nKN=np.array([(radius*radius+a*a),-Delta/1,0,a])*rho*rhobar# decide if Delta is divided by 2 or not
        mKN=-rhobar/sqrtof2*np.array([1.j*a*SinTH,0,1,1.j/SinTH])
        mbKN=np.conj(mKN)#-rho/sqrtof2*np.array([-1.j*a*SinTH,0,1,-1.j/SinTH])

        lKNd=np.transpose(np.array([-1, 0, 0, a* SinTH*SinTH])) # negate each (-1*) for (+---), but doesn't change end result
        nKNd=np.transpose(np.array([(-Delta*rho*rhobar)/1, -1, 0, (1/1)*Delta*rho*rhobar * a*SinTH*SinTH]))
        mKNd=np.transpose(np.array([-1.j*rhobar* a*SinTH/sqrtof2, 0, 1/(sqrtof2*rho), 1.j*rhobar*(radius*radius + a*a)*SinTH/sqrtof2]))
        mbKNd=np.conj(mKNd)#np.transpose(np.array([ 1.j*rho* a*SinTH/sqrtof2, 0, 1/(sqrtof2*rhobar), -1.j*rho*(radius*radius + a*a)*SinTH/sqrtof2]))

        Tet=np.array([lKN,nKN,mKN,mbKN])
        InvTet=np.array([lKNd,nKNd,mKNd,mbKNd])

        sY2m=swsh_array(2, theta, phi).reshape((2*ell+1,2*ell+1))
        # Assume z's are all real, so the conjugate relation is maintained.
        #Yn22m= z_array*sY2m[:,0] # s=-2 # not currently used in the calculation
        Yn12m= z_array*sY2m[:,1] # s=-1
        Y02m = z_array*sY2m[:,2] # s= 0
        Y12m = z_array*sY2m[:,3] # s= 1
        Y22m = z_array*sY2m[:,4] # s= 2

        # Not Dependent on q
        # q=0
        # ell=2

        # for ell>2, these need to be either arrays or functions of m
        zGqlm= Delta*Delta/12
        # my h-components are in a form that only contain zG, ThornzG, and ThornThornzG.
        # dChr will contain first derivatives of the h-comp's which are also in terms of Thorn's only
        ThornzGqlm= Deltap*Delta/6
        ThornThornzGqlm= (Deltapp*Delta+Deltap*Deltap)/6
        ThornThornThornzGqlm= Deltap#(*Deltapp/2)

        # Pieces of dGamma

        # hnn_2m= -Delta*Delta*(sqrtof6*rhobar*rhobar*Y02m+2*sqrtof2*rhobar*tau*Yn12m)
        # hnmb_2m= sqrtof2*rhobar*(rho+rhobar)*Delta*Delta*Yn12m + 2*Deltap*Delta*(sqrtof2*rhobar*Yn12m + (tau-pibar)*Yn22m)
        # hmbmb_2m= -(2*Deltap^2+Delta*(4+4*rho*Deltap))*Yn22m



        # At some point, we must sum over m
        # if ell==2
        hnndag=np.sum(2 * zGqlm * ( 2 * Y12m * pi * rhobar - sqrtof6 * Y02m * ( rho*rho ) )
            )
        hnm=np.sum(( 2 * Y12m * rho * ( ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ThornzGqlm * Y22m * ( pi - taubar ) )
            )
        hmm=np.sum(0 - Y22m * ( ThornThornzGqlm + 2 * ThornzGqlm * rhobar )
            )
        

        # Each Thornp of h contains explicit m-dependence, so I initialize them to be summed by a loop.
        Thornphnndag=0
        Thornphnm=0
        Thornphmm=0

        Thornhnndag=np.sum(2 * rho * ( - sqrtof6 * Y02m * rho * ( ThornzGqlm + 2 \
        * zGqlm * rho ) + 2 * Y12m * ( ThornzGqlm + zGqlm * ( rhobar + 2 * \
        rho ) ) * taubar ))
        
        Ethhnndag=np.sum(2 * zGqlm * ( 2 * rho * ( -3 * Y12m * rhobar * rho + ( 2 * \
        Y22m * rhobar * taubar + sqrtof6 * Y02m * ( rhobar + rho ) * taubar ) \
        ) + Y12m * ( 3 * PsiCDe2 * rhobar + ( PsiCDe2bar * rho + ( 2 * mu * \
        rhobar * ( rho - rhobar ) -2 * ( 2 * rhobar + rho ) * ( \
        taubar*taubar ) ) ) ) ))
        Ethphnndag=np.sum(4 * zGqlm * rho * ( 3 * Yn12m * ( rho*rho ) + taubar * ( \
        - sqrtof6 * Y02m * rho + Y12m * taubar ) ))
        Thornhnm=np.sum(( 2 * Y12m * rho * ( ThornThornzGqlm + ( ThornzGqlm * ( \
        rhobar + 2 * rho ) + zGqlm * ( ( rhobar*rhobar ) + rho * ( rhobar + 2 \
        * rho ) ) ) ) + Y22m * ( ThornThornzGqlm * ( pi - taubar ) + \
        ThornzGqlm * ( 2 * pi * rho - ( rhobar + rho ) * taubar ) ) ))
        
        Ethhnm=np.sum(( 4 * Y22m * rhobar * rho * ( ThornzGqlm + zGqlm * ( rhobar + \
        rho ) ) + ( -2 * Y12m * ( ThornzGqlm * ( 2 * rhobar + rho ) + zGqlm * \
        ( rhobar + rho ) * ( rhobar + 2 * rho ) ) * taubar + 0.5 * ThornzGqlm \
        * Y22m /rhobar /rho * ( 2 * mu * rhobar * ( ( \
        rho - rhobar )*( rho - rhobar ) ) + ( PsiCDe2bar * rho * ( \
        3 * rhobar + rho ) + ( PsiCDe2 * rhobar * ( rhobar + 3 * rho ) -2 * \
        ( -2 * ( rhobar*rhobar ) + ( 3 * rhobar * rho + ( rho*rho ) ) ) * ( \
        taubar*taubar ) ) ) ) ) ))
        Ethphnm=np.sum(( 1/( rhobar ) ) * ( -2 * sqrtof6 * Y02m * rhobar * ( rho*rho \
        ) * ( ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ( -2 * Y12m * rho * ( \
        - rhobar + 2 * rho ) * ( ThornzGqlm + zGqlm * ( rhobar + rho ) ) * \
        taubar + ThornzGqlm * Y22m * taubar * ( -2 * pi * rho + ( - rhobar \
        * taubar + rho * taubar ) ) ) ))
        Thornhmm=np.sum(- Y22m * ( ThornThornThornzGqlm + 2 * rhobar * ( \
        ThornThornzGqlm + ThornzGqlm * rhobar ) ))
        
        Ethhmm=np.sum(2 * Y22m * rhobar * ( ThornThornzGqlm + ThornzGqlm * rhobar ) \
        /rho * taubar)
        Ethphmm=np.sum(( 2 * ThornThornzGqlm /rhobar * rho * ( Y12m * \
        rhobar + Y22m * taubar ) + 2 * ThornzGqlm * ( 2 * Y12m * rhobar * rho \
        + Y22m * ( - rhobar + 2 * rho ) * taubar ) ))

        # only Thornp has explicit m-dependence
        #if q==0: 
        for m in [-2,-1,0,1,2]:
            Thornphnndag+= ( complex( 0,-2 ) * zGqlm * a * m * rhobar * ( rho*rho ) \
                    * ( sqrtof6 * Y02m[m+2] * rho -2 * Y12m[m+2] * taubar ) + ( 2 * sqrtof6 * \
                    Y02m[m+2] /rhobar * rho * ( ThornzGqlm * mu * rhobar + zGqlm * \
                    ( PsiCDe2bar * rho + ( rhobar * ( PsiCDe2 + ( 4 * mu * rhobar -2 * \
                    ( 4 * gammabar + mu ) * rho ) ) -2 * ( rhobar + rho ) * ( \
                    taubar*taubar ) ) ) ) + 2 * Y12m[m+2] /rho * ( -2 * ThornzGqlm \
                    * mu * rho * taubar + zGqlm * ( PsiCDe2 * ( pibar * rhobar - rho \
                    * taubar ) + 2 * rho * taubar * ( -4 * mu * rhobar + ( 6 * gammabar * \
                    rho + ( mu * rho + taubar * ( pi + taubar ) ) ) ) ) ) ) )

            Thornphnm+=( complex( 0,2 ) * Y12m[m+2] * a * m * rhobar * ( rho*rho ) * ( \
                    ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ( ThornThornzGqlm * Y22m[m+2] * \
                    ( - mubar + mu ) * ( 1/( ( rhobar - rho ) ) ) * ( pi - \
                    taubar ) + ( -2 * Y12m[m+2] /rhobar * ( 1/( ( rho - rhobar \
                    ) ) ) * ( zGqlm * ( rhobar * ( -3 * mu * ( rhobar*rhobar*rhobar ) + ( \
                    2 * ( 3 * gammabar + mu ) * ( rhobar*rhobar ) * rho + ( ( rho*rho ) * \
                    ( PsiCDe2 -6 * gammabar * rho ) + rhobar * rho * ( - PsiCDe2bar \
                    + mu * rho ) ) ) ) -2 * ( rho - rhobar ) * ( ( rhobar + rho \
                    )*( rhobar + rho ) ) * ( taubar*taubar ) ) + ( rho - rhobar ) * \
                    ( ThornThornzGqlm * mu * rhobar + ThornzGqlm * ( PsiCDe2bar * rho + ( \
                    rhobar * ( PsiCDe2 + ( 3 * mu * rhobar -6 * gammabar * rho ) ) -2 \
                    * ( rhobar + rho ) * ( taubar*taubar ) ) ) ) ) + 0.5 * ThornzGqlm * ( \
                    1/( rhobar ) ) /rho * ( complex( 0,2 ) * Y22m[m+2] * a * m * \
                    rhobar * ( rho*rho ) * ( rho - rhobar ) * taubar + Y22m[m+2] * ( \
                    taubar * ( PsiCDe2bar * rho + ( 2 * ( rho - rhobar ) * ( 4 * \
                    gammabar * rho - mu * ( 2 * rhobar + rho ) ) + ( 2 * pi * rho * \
                    taubar -2 * rhobar * ( taubar*taubar ) ) ) ) + PsiCDe2 * ( - rho \
                    * taubar + rhobar * ( pibar + taubar ) ) ) ) ) ) )

            Thornphmm+= ( complex( 0,-1 ) * Y22m[m+2] * a * m * rhobar * ( \
                    ThornThornzGqlm + 2 * ThornzGqlm * rhobar ) * rho + Y22m[m+2] * ( 1/( \
                    rhobar ) ) /rho * ( ThornThornThornzGqlm * mu * rhobar + ( \
                    2 * ThornzGqlm * rhobar * ( PsiCDe2bar * rho + ( rhobar * ( PsiCDe2 + \
                    ( mu * rhobar -4 * gammabar * rho ) ) -2 * ( rhobar + rho ) * ( \
                    taubar*taubar ) ) ) + ThornThornzGqlm * ( PsiCDe2bar * rho + ( rhobar \
                    * ( PsiCDe2 + ( 2 * mu * rhobar -4 * gammabar * rho ) ) -2 * ( \
                    rhobar + rho ) * ( taubar*taubar ) ) ) ) ) )

        # else:
        #     #q=0 
        #     for m in [-2,-1,0,1,2]:
        #         Thornphnndag += -2 * ( complex( 0,1 ) * zGqlm * a * m * rhobar * ( \
        #             rho*rho ) * ( sqrtof6 * Y02m[m+2] * rho -2 * Y12m[m+2] * taubar ) + ( ( 3/2 \
        #             )**( 0.5 ) * Y02m[m+2] /rhobar * rho * ( -2 * ThornzGqlm * mu * \
        #             rhobar + zGqlm * ( PsiCDe2bar * ( -2 + q ) * rho + ( rhobar * ( \
        #             PsiCDe2 * ( -2 + q ) + ( 2 * ( -4 + q ) * mu * rhobar + ( -4 * ( -4 + \
        #             q ) * gammabar * rho + ( 4 -2 * q ) * mu * rho ) ) ) -2 * ( -2 + \
        #             q ) * ( rhobar + rho ) * ( taubar*taubar ) ) ) ) + Y12m[m+2] * ( 1/( rho ) \
        #             ) * ( 2 * ThornzGqlm * mu * rho * taubar + zGqlm * ( PsiCDe2 * ( -1 + \
        #             q ) * ( pibar * rhobar - rho * taubar ) + rho * taubar * ( ( 8 + \
        #             -2 * q ) * mu * rhobar + ( 4 * ( -3 + q ) * gammabar * rho + ( 2 * ( \
        #             -1 + q ) * mu * rho + 2 * ( -1 + q ) * taubar * ( pi + taubar ) ) ) ) \
        #             ) ) ) )

        #         Thornphnm+= ( -2 * rho * ( Y12m[m+2] * ( ( -4 + q ) * gammabar + ( -2 + q ) \
        #             * gamma ) + complex( 0,-1 ) * Y12m[m+2] * a * m * rhobar * rho ) * ( \
        #             ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ( ThornzGqlm * ( - Y22m[m+2] \
        #             * ( -4 * gammabar + q * ( gammabar + gamma ) ) + complex( 0,1 ) * \
        #             Y22m[m+2] * a * m * rhobar * rho ) * ( pi - taubar ) + ( Y12m[m+2] * ( 1/( \
        #             rhobar ) ) * ( 1/( ( rho - rhobar ) ) ) * ( ( rho - rhobar \
        #             ) * ( -2 * ThornThornzGqlm * mu * rhobar + ThornzGqlm * ( PsiCDe2 * ( \
        #             -2 + q ) * rhobar + ( PsiCDe2bar * ( -2 + q ) * rho + ( -2 * mu * \
        #             rhobar * ( rhobar + 2 * rho ) -2 * ( -2 + q ) * ( rhobar + rho ) * \
        #             ( taubar*taubar ) ) ) ) ) + zGqlm * ( PsiCDe2 * ( -2 + q ) * rhobar * \
        #             ( rho*rho ) + ( 2 * mu * rhobar * ( ( rhobar*rhobar*rhobar ) + ( \
        #             rhobar * ( rho*rho ) -2 * ( rho*rho*rho ) ) ) + ( -2 + q ) * ( - \
        #             PsiCDe2bar * ( rhobar*rhobar ) * rho -2 * ( rho - rhobar ) * ( \
        #             ( rhobar + rho )*( rhobar + rho ) ) * ( taubar*taubar ) ) ) ) ) + \
        #             Y22m[m+2] * ( ThornzGqlm * ( 2 * mubar * taubar - mu * ( pi + taubar \
        #             ) ) + ( pi - taubar ) * ( ThornThornzGqlm * ( mubar - mu ) \
        #             * ( 1/( ( rho - rhobar ) ) ) + 0.5 * ThornzGqlm * ( -1 + q ) * ( \
        #             1/( rhobar ) ) /rho * ( PsiCDe2 * rhobar + ( PsiCDe2bar * \
        #             rho + 2 * ( rhobar + rho ) * taubar * tau ) ) ) ) ) ) )

        #         Thornphmm+= ( - ( ThornThornzGqlm + 2 * ThornzGqlm * rhobar ) * ( -1 \
        #             * Y22m[m+2] * ( -4 * gammabar + q * ( gammabar + gamma ) ) + complex( 0,1 \
        #             ) * Y22m[m+2] * a * m * rhobar * rho ) -0.5 * Y22m[m+2] /rhobar * \
        #             ( 1/( rho ) ) * ( -2 * ThornThornThornzGqlm * mu * rhobar + ( \
        #             ThornThornzGqlm * ( rhobar * ( PsiCDe2 * ( -2 + q ) -4 * mu * \
        #             rhobar ) + ( PsiCDe2bar * ( -2 + q ) * rho -2 * ( -2 + q ) * ( \
        #             rhobar + rho ) * ( taubar*taubar ) ) ) + 2 * ThornzGqlm * rhobar * ( \
        #             rhobar * ( PsiCDe2 * ( -2 + q ) -2 * mu * rhobar ) + ( PsiCDe2bar * \
        #             ( -2 + q ) * rho -2 * ( -2 + q ) * ( rhobar + rho ) * ( \
        #             taubar*taubar ) ) ) ) ) )

        hnn=np.conj(hnndag)
        hnmb=np.conj(hnm)
        hmbmb=np.conj(hmm)

        Thornhnn=np.conj(Thornhnndag)
        # Thornhnmb=np.conj(Thornhnm)
        Thornhmbmb=np.conj(Thornhmm)

        Thornphnn=np.conj(Thornphnndag)
        Thornphnmb=np.conj(Thornphnm)
        Thornphmbmb=np.conj(Thornphmm)

        Ethhnn=np.conj(Ethphnndag)
        Ethhnmb=np.conj(Ethphnm)
        Ethhmbmb=np.conj(Ethphmm)

        Ethphnn=np.conj(Ethhnndag)
        # Ethphnmb=np.conj(Ethhnm)
        Ethphmbmb=np.conj(Ethhmm)

        # PDhnn_m= 
        # PDhnmb_m= 
        # PDhmbmb_m= 



        dChr=np.zeros((4,4,4),np.complex128) # complex numbers are being cast as real
        # Some terms have explicit q-dependence, but all boost-weight dependence vanishes in the final summation.
        dChr[0,0,1]=0.5 * ( Thornhnn + Thornhnndag) - ( hnmb * pibar + hnm * pi )
        dChr[0,0,2]=0.5 * ( Thornhnm +  hnm * rhobar - hmm * ( pi + taubar ) )
        dChr[0,0,3]=np.conj(dChr[0,0,2]) # 0.5 * ( Thornhnmb + ( hnmb * rho - hmbmb * ( pibar + tau ) ) )
        dChr[0,1,0]=dChr[0,0,1]
        dChr[0,1,1]=0.5 * ( Thornphnn + Thornphnndag )
        dChr[0,1,2]=0.5 * ( Ethhnn + Ethhnndag ) - hnm * mu
        dChr[0,1,3]=np.conj(dChr[0,1,2])
        dChr[0,2,0]=dChr[0,0,2]
        dChr[0,2,1]=dChr[0,1,2]
        dChr[0,2,2]=( Ethhnm -0.5 * Thornphmm - ( hmm * mu + hnm * tau ) )
        dChr[0,2,3]=0.5 * ( Ethhnmb + Ethphnm + ( hnn + hnndag ) * ( rhobar + rho ) - ( hnm * taubar + hnmb * tau ) )
        dChr[0,3,0]=dChr[0,0,3]
        dChr[0,3,1]=dChr[0,1,3]
        dChr[0,3,2]=dChr[0,2,3]
        dChr[0,3,3]=np.conj(dChr[0,2,2]) # ( Ethphnmb + ( -0.5 * Thornphmbmb + ( 0 - hmbmb * mubar - hnmb * taubar ) ) )

        dChr[1,1,1]=( 0.5 * ( 0 -  Thornhnn - Thornhnndag ) + ( hnm * ( pi + taubar ) + hnmb * ( pibar + tau ) ) )
        dChr[1,1,2]=0.5 * ( 0 - Thornhnm + ( hnm * rhobar + hmm * ( pi + taubar ) ) )
        dChr[1,1,3]=np.conj(dChr[1,1,2])

        dChr[1,2,1]=dChr[1,1,2]
        dChr[1,2,2]=( -0.5 * Thornhmm + hmm * rhobar )

        dChr[1,3,1]=dChr[1,1,3]

        dChr[1,3,3]=np.conj(dChr[1,2,2])

        dChr[2,0,1]=dChr[1,1,3] - hmbmb*tau # 0.5 * ( 0 - Thornhnmb + ( hnmb * rho + hmbmb * ( pibar - tau ) ) )

        dChr[2,0,3]=-0.5 * Thornhmbmb
        dChr[2,1,0]=dChr[2,0,1]
        dChr[2,1,1]=0.5 * ( Ethphnn + Ethphnndag ) - ( Thornphnmb + hnmb * mubar + ( hnn + hnndag ) * taubar )
        dChr[2,1,2]=0.5 * ( Ethphnm - Ethhnmb + ( hnn + hnndag ) * ( rho - rhobar ) - ( hnm * taubar + hnmb * tau ) )
        dChr[2,1,3]=( -0.5 * Thornphmbmb - hnmb * taubar )

        dChr[2,2,1]=dChr[2,1,2] # 0.5 * ( 0 - Ethhnmb + ( Ethphnm + ( 0 - hnn * rhobar + ( 0 - hnndag * rhobar + ( hnn * rho + ( hnndag * rho + ( 0 - hnm * taubar - hnmb * tau ) ) ) ) ) ) )
        dChr[2,2,2]=( 0.5 * Ethphmm + hnm * ( rho - rhobar ) )
        dChr[2,2,3]=( -0.5 * Ethhmbmb - hnmb * rhobar )
        dChr[2,3,0]=dChr[2,0,3]
        dChr[2,3,1]=dChr[2,1,3]
        dChr[2,3,2]=dChr[2,2,3]
        dChr[2,3,3]=-0.5 * Ethphmbmb

        dChr[3,0,1]=np.conj(dChr[2,0,1]) #dChr[1,1,2] -hmm*taubar # 0.5 * ( 0 - Thornhnm + ( hnm * rhobar + hmm * ( pi - taubar ) ) )
        dChr[3,0,2]=-0.5 * Thornhmm # np.conj(dChr[2,0,3])

        dChr[3,1,0]=dChr[3,0,1]
        dChr[3,1,1]=np.conj(dChr[2,1,1])
        dChr[3,1,2]=np.conj(dChr[2,1,3]) # ( -0.5 * Thornphmm - hnm * tau )
        dChr[3,1,3]=np.conj(dChr[2,1,2])
        dChr[3,2,0]=-0.5 * Thornhmm
        dChr[3,2,1]=dChr[3,1,2]
        dChr[3,2,2]=-0.5 * Ethhmm
        dChr[3,2,3]=np.conj(dChr[2,3,2]) # ( -0.5 * Ethphmm - hnm * rho )

        dChr[3,3,1]=dChr[3,1,3]
        dChr[3,3,2]=dChr[3,2,3]
        dChr[3,3,3]=np.conj(dChr[2,2,2]) # ( 0.5 * Ethhmbmb + hnmb * ( rhobar - rho ) )

        # def etaouter(){

        # }

        # transform to coordinate basis
        dChrKN=np.zeros((4,4,4),np.complex128)
        for ii in range(4): 
            for jj in range(4): 
                for kk in range(4): 
                    lnis=[1,0,3,2] # apply eta matrix for (-+++), change ">1"->"<2" for (+---) or give overall minus
                    signi = 1 if ii>1 else -1 
                    signj = 1 if jj>1 else -1 
                    signk = 1 if kk>1 else -1 
                    dChrKN += signi*signj*signk*np.einsum('i,j,k',Tet[lnis[ii]],InvTet[lnis[jj]],InvTet[lnis[kk]]) * dChr[ii, jj, kk]
        
        return Chr+dChrKN
        #return np.real(dChrKN)

## Average computation time for this function is > ~2ms
#print((time.perf_counter()-start)/10000)

    return Chr

    # totalG=Chr+dChr

# Average time for totalG is about the same
# print((time.perf_counter()-start)/10000)


a=.9
eps_ar=[10**(-5),10**(-6),10**(-9),10**(-12)]
# 1e-5 diverges around t=6500 for r=12.8
eps=eps_ar[3]
zs = eps*np.array([2+1.1j,-3.222,1.03+2.331j,-.1-3.2j,-1.54]) # eps X z_array with Norm ~ O(1)
# strong perturbations will cause accelerations too great for scipy's solver to integrate

# imtolfactor=10**6 # imaginary parts less than imtolfactor*machineprecision will be suppressed.
# #As far as I have currently found, the largest variances of dChr[a,b,c] from dChr[a,c,b], or any Im[] from zero, or q=0 from q!=0 are ~e-11.

# for it in range(0,10000):

#     t=0.001*it
#     radius=4+0.0012*it
#     theta=0.0001+3.14*it/10000
#     phi=4.4*it/10000    
#     ChrTot=TotalGamma(a,radius,theta,phi,zs)
#     # Check for complex values
#     it1, it2, it3 = 0,0,0
#     for x in ChrTot: 
#         for y in x: 
#             for z in y:
#                 if np.iscomplex(np.real_if_close(z,imtolfactor)):
#                     print((it1,it2,it3),z,ChrTot[it1,it2,it3]-ChrTot[it1,it3,it2])
#                 it3+=1
#                 it3=it3%4
#             it2+=1
#             it2=it2%4
#         it1+=1
#         it1=it1%4

            

# print(np.real_if_close(dChrKN[0,0,0],100),np.real_if_close(dChrKN[0,0,1],100),np.real_if_close(dChrKN[0,0,2],100),np.real_if_close(dChrKN[0,0,3],100))
# print(np.real_if_close(dChrKN[0,1,0],100),np.real_if_close(dChrKN[0,1,1],100),np.real_if_close(dChrKN[0,1,2],100),np.real_if_close(dChrKN[0,1,3],100))
# print(np.real_if_close(dChrKN[0,2,0],100),np.real_if_close(dChrKN[0,2,1],100),np.real_if_close(dChrKN[0,2,2],100),np.real_if_close(dChrKN[0,2,3],100))
# print(np.real_if_close(dChrKN[0,3,0],100),np.real_if_close(dChrKN[0,3,1],100),np.real_if_close(dChrKN[0,3,2],100),np.real_if_close(dChrKN[0,3,3],100))
# dChrKN = np.real_if_close(dChrKN,1)
# print("dChrKN calculated ",time.perf_counter()-start,time.perf_counter()-restart)

p,e,x=5.7,0,.5
r1,r2=p/(1+e),p/(1-e)
zm = np.sqrt(1 - x*x)
#Equatorial and Circular
if e==0 and x*x==1:
    En=((-2 + p)* np.sqrt(p) + a/x)/np.sqrt(2* a/x *p**(3/2) + (-3 + p)* p*p)
    Lz=(a*a - 2 *a/x *np.sqrt(p) + p*p)/np.sqrt(2 *a/x + (-3 + p)* np.sqrt(p))/p**(3/4)
    Q=0
elif e==0: 
    En=np.sqrt(
        (
         (-3 + p)* (-2 + p)*(-2 + p) *p**5 - 2 *a**5 *x *(-1 + x*x)*
                np.sqrt(p*p*p + a*a *p *(-1 + x*x))
         + a*a*a*a *p*p *(-1 + x*x) *
                (4 - 5* p* (-1 + x*x) + 3* p*p *(-1 + x*x))
         - a**6 *(-1 + x*x)*(-1 + x*x) *
                (x*x + p*p *(-1 + x*x) - p *(1 + 2 *x*x))
         + a*a *p*p*p *
                (4 - 4 *x*x + p *(12 - 7 *x*x) - 3 *p*p*p *(-1 + x*x) + p*p *(-13 + 10 *x*x))
         + a*
                (-2* p**(9/2) *x *
                  np.sqrt(p*p + a*a *(-1 + x*x))
                 + 4* p*p*p *x *
                  np.sqrt(p*p*p + a*a *p *(-1 + x*x))
                )
         + 2 *a*a*a *
                (2 *p *x *(-1 + x*x) *
                  np.sqrt(p*p*p + a*a *p *(-1 + x*x))
                 - x*x*x *
                  np.sqrt(p**7 + a*a *p**5 * (-1 + x*x))
                )
        )/(
            (p*p - a*a *(-1 + x*x)) *
                ((-3 + p)*(-3 + p) *p*p*p*p - 2 *a*a *p*p *
                   (3 + 2 *p - 3 *x*x + p*p *(-1 + x*x))
                 + a*a*a*a *(-1 + x*x) *
                   (-1 + x*x + p*p *(-1 + x*x) - 2 *p *(1 + x*x))
                )
          )
    )
    
    g = 2 *a *p
    d = (a*a + (-2 + p) *p) *(p*p - a*a *(-1 + x*x))
    h = ((-2 + p) *p - a*a *(-1 + x*x))/(x*x)
    f = p*p*p*p + a*a *(p *(2 + p) - (a*a + (-2 + p) *p) *(-1 + x*x))
    Lz= (-En *g + x *np.sqrt((-d *h + En*En *(g*g + f *h))/(x*x)))/h
    Q = zm*zm *(a*a * (1 - En*En) + Lz*Lz/(1 - zm*zm))

Phisq = Lz *Lz / (1 - En *En)
q = Q / (1 - En *En) # is mu^2 1 or -1?
U0sq = 1-x*x # zm^2
asqU1sq = ((a *a + q + Phisq) + np.sqrt((a *a + q + Phisq) ** 2 - 4 * a *a * q)) / (2) #* a**2); #u0 and u1 just need to be swapped. now u1^2>1

t0,r0,th0,ph0=0,r1,np.arccos(zm),0

usq=np.cos(th0)*np.cos(th0)
Delta=r0 * r0 - 2 * r0 + a * a
Sigma=r0*r0+a*a*usq

def base_RHS(t, y):
    time,radius,theta,phi=y[0:4]
    u=y[4:8]
    #xdot=u
    udot=-1*np.einsum('ijk,j,k->i',TotalGamma(a,radius,theta,phi),u,u)
    return np.concatenate((u,udot))
def pert_RHS(t, y, z1,z2,z3,z4,z5):
    time,radius,theta,phi=y[0:4]
    u=y[4:8]
    #xdot=u
    udot=-1*np.einsum('ijk,j,k->i',TotalGamma(a,radius,theta,phi,np.array([z1,z2,z3,z4,z5])),u,u)
    return np.concatenate((u,udot))


td0,rd0,thd0,phd0=(a * (Lz - a * En * (1 - usq)) + (r0 *r0 + a *a) * ((r0 *r0 + a *a) * En - Lz * a) / Delta) / Sigma,0,0,(2 * M * r0 * a * En + (Sigma - 2 * M * r0) * Lz / (1 - usq)) / (Delta * Sigma)


#print(time.perf_counter()-start)

ICs=[t0,r0,th0,ph0,td0,rd0,thd0,phd0]

tau_end,Nsteps=10000,5


import matplotlib.pyplot as plt

#make plots testing errors and divergences for correlation with integrator tolerance level
exps=np.arange(7.,13.) #exponents for the integration tolerances
t_div=np.zeros(exps.size) # time what integration of perturbed orbit ended noting when the integration ends early do to divergence
y_div=np.zeros((exps.size,8)) # pert_sol.y at t_div fpr each tolerance level
# del_r_max=np.zeros(6) # tracker of max variance of base orbit from r0 for when e=0, contains the max vals of del_r for each exp of tolerance
# t_delr_max=np.zeros(6) # the times for del_r_max
# del_theta_max=np.zeros(6) # tracker of max variance of base orbit from th0 for when x=+/-1, contains the max vals of del_theta for each exp of tolerance
# t_delth_max=np.zeros(6) # the times for del_theta_max
del_r_arrays=[]
del_theta_arrays=[]
base_t_arrays=[]
base_r_arrays=[]
base_theta_arrays=[]
base_phi_arrays=[]
pert_t_arrays=[]
pert_r_arrays=[]
pert_theta_arrays=[]
pert_phi_arrays=[]
# pert_ys=[]
times=np.linspace(0,tau_end,tau_end)
for i in range(exps.size):
    rtol, atol = 10.**(-exps[i]), 10.**(-exps[i])
    # start=time.perf_counter()
    # base_sol=ODE(base_RHS,[0,tau_end],ICs,t_eval=times,rtol=rtol,atol=atol)
    # # base_t_arrays.append(base_sol.t)
    # base_r_arrays.append(base_sol.y[1])
    # base_theta_arrays.append(base_sol.y[2])
    # base_phi_arrays.append(base_sol.y[3])
    # # t_array_len= base_sol.t.size
    # # print("tau_end=",base_sol.t[-1],"base_sol runtime=",time.perf_counter()-start,"sol.t length=",base_sol.t.shape)
    # # print("base size=",base_sol.y.shape,base_sol.nfev)
    # if e==0:
    #     del_r=[0.]
    #     for j in range(1,base_sol.t.size):
    #         r=base_sol.y[1][j]

    #         # if abs(del_r[-1])<abs(r-r0): 
    #         #     del_r_max[i]=r-r0
    #         #     t_delr_max[i]=base_sol.t[j]
            
    #         del_r.append(r-r0)
    #     del_r_arrays.append(del_r)

    # if x*x==1:
    #     del_theta=[0.]
    #     for j in range(1,base_sol.t.size):
    #         theta=base_sol.y[2][j]
            
    #         # if abs(del_theta[-1])<abs(theta-th0): 
    #         #     del_theta_max[i]=theta-th0
    #         #     t_delth_max[i]=base_sol.t[j]
            
    #         del_theta.append(theta-th0)
    #     del_theta_arrays.append(del_theta)
    # # print("Del_r=",del_r)

    # restart=time.perf_counter()

    pert_sol=ODE(pert_RHS,[0,tau_end],ICs,t_eval=times,args=zs,rtol=rtol,atol=atol)#) Doing this makes it take significantly longer
    pert_t_arrays.append(pert_sol.t)
    pert_r_arrays.append(pert_sol.y[1])
    pert_theta_arrays.append(pert_sol.y[2])
    pert_phi_arrays.append(pert_sol.y[3])
    
    # pert_ys.append(pert_sol.y)
    # pt_array_len= pert_sol.t.size
    # if pert_sol.t[-1]<tau_end:
    #     t_div[i]=pert_sol.t[-1]
    #     y_div[i]=[val[-1] for val in pert_sol.y]
    # else:
    #     t_div[i]=np.nan
    #     y_div[i]=np.nan

    # if e==0:
    #     del_r=[0.]
    #     for j in range(1,pert_sol.t.size):
    #         r=pert_sol.y[1][j]

    #         if abs(del_r[-1])<abs(r-r0): 
    #             del_r_max[i]=r-r0
    #             t_delr_max[i]=base_sol.t[j]
            
    #         del_r.append(r-r0)
    #     del_r_arrays.append(del_r)

    # if x*x==1:
    #     del_theta=[0.]
    #     for j in range(1,base_sol.t.size):
    #         theta=base_sol.y[2][j]
            
    #         if abs(del_theta[-1])<abs(theta-th0): 
    #             del_theta_max[i]=theta-th0
    #             t_delth_max[i]=base_sol.t[j]
            
    #         del_theta.append(theta-th0)
    #     del_theta_arrays.append(del_theta)

# #fig_del_r for exp in exps: 6 plots of del_r_arrays[i] vs base_t_arrays[i]

# #fig_del_th for exp in exps: 6 plots of del_theta_arrays[i] vs base_t_arrays[i]
# fig_del, (ax_del_r,ax_del_th)=plt.subplots(1,2)
# for exp in range(exps.size):
#     ax_del_r.plot(times,del_r_arrays[exp],label="1e-%d" % exps[exp])
#     if x*x==1: ax_del_th.plot(times,del_theta_arrays[exp],label="1e-%d" % exps[exp])
# ax_del_r.set_xlabel('proper time')
# ax_del_r.set_ylabel('$\delta$r')
# ax_del_th.set_xlabel('proper time')
# ax_del_th.set_ylabel('$\delta\\theta$')
# ax_del_r.legend()
# ax_del_th.legend()
# fig_del.suptitle("Unperturbed Deviations")
# plt.show()

# #fig_max_dels 2 scatter plots: del_r_max and t_delr_max vs tol and del_theta_max and t_delth_max vs tol
# fig_max_dels, (ax_dr_max,ax_dth_max)=plt.subplots(1,2)
sizes=10*exps
# cb1=ax_dr_max.scatter(t_delr_max,del_r_max,s=sizes,c=exps,cmap='RdBu_r')
# ax_dr_max.set_title('Max $\delta$r')
# fig_max_dels.colorbar(cb1,ax=ax_dr_max)
# if x*x==1: 
#     cb2=ax_dth_max.scatter(t_delth_max,del_theta_max,s=sizes,c=exps,cmap='RdBu_r')
#     ax_dth_max.set_title('Max $\delta\\theta$')
#     fig_max_dels.colorbar(cb2,ax=ax_dth_max)
# fig_max_dels.suptitle("Tolerance Exp = marker size")
# plt.show()

#fig_pert_sols: 8x6 pert_sol.y[j][i] vs pert_t_arrays[i]
fig_pert_sols, axs_pert=plt.subplots(3,1)
# for c in range(4): # columns (not t),r,th,ph
for e in range(exps.size): # tol exp
    axs_pert[0].plot(pert_t_arrays[e],pert_r_arrays[e],label="1e-%d" % exps[e])
    axs_pert[1].plot(pert_t_arrays[e],pert_theta_arrays[e],label="1e-%d" % exps[e])
    axs_pert[2].plot(pert_t_arrays[e],pert_phi_arrays[e],label="1e-%d" % exps[e])
        # axs_pert[1,c].plot(pert_t_arrays[e],del_r_arrays[e],label="1e-%f" % exps[e])# row for velocity
axs_pert[0].set_xlabel('proper time')
axs_pert[0].set_ylabel('radius')
axs_pert[0].legend()
axs_pert[1].set_xlabel('proper time')
axs_pert[1].set_ylabel('$\\theta$')
axs_pert[1].legend()
axs_pert[2].set_xlabel('proper time')
axs_pert[2].set_ylabel('$\phi$')
axs_pert[2].legend()
fig_pert_sols.suptitle("Perturbed Solution for eps=%f"%eps)
plt.show()

# #fig_pert_sols: 8x6 pert_sol.y[j][i] vs pert_t_arrays[i]
# fig_pert_sols, axs_pert=plt.subplots(1,3)
# # for c in range(4): # columns (not t),r,th,ph
# for e in range(exps.size): # tol exp
#     axs_pert[0].plot(times,pert_r_arrays[e]-base_r_arrays[e],label="1e-%d" % exps[e])
#     axs_pert[1].plot(times,pert_theta_arrays[e]-base_theta_arrays[e],label="1e-%d" % exps[e])
#     axs_pert[2].plot(times,pert_phi_arrays[e]-base_phi_arrays[e],label="1e-%d" % exps[e])
#         # axs_pert[1,c].plot(pert_t_arrays[e],del_r_arrays[e],label="1e-%f" % exps[e])# row for velocity
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

fig_pert_sols, axs_pert=plt.subplots(3,1)
# for c in range(4): # columns (not t),r,th,ph
for e in range(exps.size-1): # tol exp
    axs_pert[0].plot(pert_t_arrays[e],pert_r_arrays[-1]-pert_r_arrays[e],label="1e-%d" % exps[e])
    axs_pert[1].plot(pert_t_arrays[e],pert_theta_arrays[-1]-pert_theta_arrays[e],label="1e-%d" % exps[e])
    axs_pert[2].plot(pert_t_arrays[e],pert_phi_arrays[-1]-pert_phi_arrays[e],label="1e-%d" % exps[e])
        # axs_pert[1,c].plot(pert_t_arrays[e],del_r_arrays[e],label="1e-%f" % exps[e])# row for velocity
axs_pert[0].set_xlabel('proper time')
axs_pert[0].set_ylabel('radius')
axs_pert[0].legend()
axs_pert[1].set_xlabel('proper time')
axs_pert[1].set_ylabel('$\\theta$')
axs_pert[1].legend()
axs_pert[2].set_xlabel('proper time')
axs_pert[2].set_ylabel('$\phi$')
axs_pert[2].legend()
fig_pert_sols.suptitle("Perturbed by eps=%f Delta from tol e-12"%eps)
plt.show()

# #fig_t_div plot of t_div vs tol_exp
# fig_t_div, ax_t_div=plt.subplots()
# fig_t_div.suptitle('Proper Time to Divergence')
# ax_t_div.stem(exps,t_div,use_line_collection=True)
# ax_t_div.set_xlabel('tol=1e-X')
# ax_t_div.set_ylabel('proper time')
# plt.show()

# #fig_pert_div: 8 subplots y_div vs t_div scatter plot labeled by exp
# fig_y_div, axs_y_div=plt.subplots(1,2)
# cbrdiv=axs_y_div[0].scatter(t_div,[val[1] for val in y_div],s=sizes,c=exps,label='r')
# # axs_y_div[0].scatter(t_div,[val[5] for val in y_div],marker='d',s=sizes,c=exps,label='dr/d$\\tau$')
# axs_y_div[0].set_xlabel('time of divergence')
# axs_y_div[0].legend()
# fig_y_div.colorbar(cbrdiv,ax=axs_y_div[0],extend='both')
# cbthdiv=axs_y_div[1].scatter(t_div,[val[2] for val in y_div],s=sizes,c=exps,label='$\\theta$')
# # axs_y_div[1].scatter(t_div,[val[6] for val in y_div],marker='d',s=sizes,c=exps,label='d$\\theta$/d$\\tau$')
# axs_y_div[1].set_xlabel('time of divergence')
# axs_y_div[1].legend()
# fig_y_div.colorbar(cbthdiv,ax=axs_y_div[1],extend='both')
# plt.show()

# print("tau_end=",pert_sol.t[-1],"pert_sol runtime=",time.perf_counter()-restart,"sol.t length=",pert_sol.t.shape)
# print("pert size=",pert_sol.y.shape)
# print(pert_sol.nfev)
# if pert_sol.t[-1]<tau_end: print("Integration ended early because the acceleration exceeded the solver's capacity. This probably means the perturbation destabilized the orbit.")


# from math import floor
# print("base_sol.phi: {}".format([base_sol.y[3][floor(n)] for n in range(0,t_array_len,floor(t_array_len/Nsteps) )]))
# print("pert_sol.phi: {}".format([pert_sol.y[3][floor(n)] for n in range(0,pt_array_len,floor(pt_array_len/Nsteps) )]))
# print("base_sol.r: {}".format([base_sol.y[1][floor(n)] for n in range(0,t_array_len,floor(t_array_len/Nsteps) )]))
# print("pert_sol.r: {}".format([pert_sol.y[1][floor(n)] for n in range(0,pt_array_len,floor(pt_array_len/Nsteps) )]))
