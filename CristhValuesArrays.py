import time


import numpy as np
import spherical
import quaternionic



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
# restart=start

M=1
a=.9

zs=[223,-47,100,-559,-15]

imtolfactor=10*6 # imaginary parts less than imtolfactor*machineprecision will be suppressed.
# As far as I have currently found, the largest variances of dChr[a,b,c] from dChr[a,c,b], or any Im[] from zero, or q=0 from q!=0 are ~e-11.

for it in range(5999,6000):

    t=0.001*it
    radius=4+0.005*it
    theta=1.0001+3.14*it/10000
    phi=4.4*it/10000

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

    # Kinnersley Tetrad in Outgoing? Kerr-Newman Coordinates
    lKN=np.array([0,1,0,0])
    nKN=np.array([(radius*radius+a*a),-Delta/1,0,a])*rho*rhobar# decide if Delta is divided by 2 or not
    mKN=-rhobar/sqrtof2*np.array([1.j*a*SinTH,0,1,1.j/SinTH])
    mbKN=np.conj(mKN)#-rho/sqrtof2*np.array([-1.j*a*SinTH,0,1,-1.j/SinTH])

    lKNd=np.transpose(np.array([1, 0, 0, -a* SinTH*SinTH]))
    nKNd=np.transpose(np.array([(Delta*rho*rhobar)/1, 1, 0, -(1/1)*Delta*rho*rhobar * a*SinTH*SinTH]))
    mKNd=np.transpose(np.array([1.j*rhobar* a*SinTH/sqrtof2, 0, -1/(sqrtof2*rho), -1.j*rhobar*(radius*radius + a*a)*SinTH/sqrtof2]))
    mbKNd=np.conj(mKNd)#np.transpose(np.array([ -1.j*rho* a*SinTH/sqrtof2, 0, -1/(sqrtof2*rhobar), 1.j*rho*(radius*radius + a*a)*SinTH/sqrtof2]))

    Tet=np.array([lKN,nKN,mKN,mbKN])
    InvTet=np.array([lKNd,nKNd,mKNd,mbKNd])

    # Non-zero Spin Coefficients (in terms of rho)
    tau= -rho*rhobar*(1.j*a*SinTH)/sqrtof2
    taubar= -tau
    pi= taubar*rho/rhobar # is np.conj faster?
    pibar= tau*rhobar/rho
    beta= (-rhobar/(2*sqrtof2)*CosTH/SinTH)
    betabar= beta*rho/rhobar
    # eps=0
    # epsbar=0
    mu= rho*rho*rhobar*Delta/2
    mubar= mu*rhobar/rho
    gamma= mu + rho*rhobar*Deltap/4
    gammabar= gamma-mu+mubar
    alpha= pi-betabar
    alphabar= pibar-beta




    # print("Ylms func ",time.perf_counter()-start,time.perf_counter()-restart)
    # Ylmtimer=time.perf_counter()


    sY2m=swsh_array(2, theta, phi).reshape((2*ell+1,2*ell+1))

    # Assume z's are all real.
    Yn22m=zs*sY2m[:,0]
    Yn12m=zs*sY2m[:,1]
    Y02m=zs*sY2m[:,2]
    Y12m=zs*sY2m[:,3]
    Y22m=zs*sY2m[:,4]

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

    Thornhnndag=np.sum(2 * rho * ( -1 * sqrtof6 * Y02m * rho * ( ThornzGqlm + 2 \
    * zGqlm * rho ) + 2 * Y12m * ( ThornzGqlm + zGqlm * ( rhobar + 2 * \
    rho ) ) * taubar ))
    
    Ethhnndag=np.sum(2 * zGqlm * ( 2 * rho * ( -3 * Y12m * rhobar * rho + ( 2 * \
    Y22m * rhobar * taubar + sqrtof6 * Y02m * ( rhobar + rho ) * taubar ) \
    ) + Y12m * ( 3 * PsiCDe2 * rhobar + ( PsiCDe2bar * rho + ( 2 * mu * \
    rhobar * ( -1 * rhobar + rho ) + -2 * ( 2 * rhobar + rho ) * ( \
    taubar*taubar ) ) ) ) ))
    Ethphnndag=np.sum(4 * zGqlm * rho * ( 3 * Yn12m * ( rho*rho ) + taubar * ( \
    -1 * sqrtof6 * Y02m * rho + Y12m * taubar ) ))
    Thornhnm=np.sum(( 2 * Y12m * rho * ( ThornThornzGqlm + ( ThornzGqlm * ( \
    rhobar + 2 * rho ) + zGqlm * ( ( rhobar*rhobar ) + rho * ( rhobar + 2 \
    * rho ) ) ) ) + Y22m * ( ThornThornzGqlm * ( pi + -1 * taubar ) + \
    ThornzGqlm * ( 2 * pi * rho + -1 * ( rhobar + rho ) * taubar ) ) ))
    
    Ethhnm=np.sum(( 4 * Y22m * rhobar * rho * ( ThornzGqlm + zGqlm * ( rhobar + \
    rho ) ) + ( -2 * Y12m * ( ThornzGqlm * ( 2 * rhobar + rho ) + zGqlm * \
    ( rhobar + rho ) * ( rhobar + 2 * rho ) ) * taubar + 0.5 * ThornzGqlm \
    * Y22m * ( 1/( rhobar ) ) * ( 1/( rho ) ) * ( 2 * mu * rhobar * ( ( \
    -1 * rhobar + rho )*( -1 * rhobar + rho ) ) + ( PsiCDe2bar * rho * ( \
    3 * rhobar + rho ) + ( PsiCDe2 * rhobar * ( rhobar + 3 * rho ) + -2 * \
    ( -2 * ( rhobar*rhobar ) + ( 3 * rhobar * rho + ( rho*rho ) ) ) * ( \
    taubar*taubar ) ) ) ) ) ))
    Ethphnm=np.sum(( 1/( rhobar ) ) * ( -2 * sqrtof6 * Y02m * rhobar * ( rho*rho \
    ) * ( ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ( -2 * Y12m * rho * ( \
    -1 * rhobar + 2 * rho ) * ( ThornzGqlm + zGqlm * ( rhobar + rho ) ) * \
    taubar + ThornzGqlm * Y22m * taubar * ( -2 * pi * rho + ( -1 * rhobar \
    * taubar + rho * taubar ) ) ) ))
    Thornhmm=np.sum(-1 * Y22m * ( ThornThornThornzGqlm + 2 * rhobar * ( \
    ThornThornzGqlm + ThornzGqlm * rhobar ) ))
    
    Ethhmm=np.sum(2 * Y22m * rhobar * ( ThornThornzGqlm + ThornzGqlm * rhobar ) \
    * ( 1/( rho ) ) * taubar)
    Ethphmm=np.sum(( 2 * ThornThornzGqlm * ( 1/( rhobar ) ) * rho * ( Y12m * \
    rhobar + Y22m * taubar ) + 2 * ThornzGqlm * ( 2 * Y12m * rhobar * rho \
    + Y22m * ( -1 * rhobar + 2 * rho ) * taubar ) ))

    # only Thornp has explicit m-dependence
    #if q==0: 
    for m in [-2,-1,0,1,2]:
        Thornphnndag+= ( complex( 0,-2 ) * zGqlm * a * m * rhobar * ( rho*rho ) \
                * ( sqrtof6 * Y02m[m+2] * rho + -2 * Y12m[m+2] * taubar ) + ( 2 * sqrtof6 * \
                Y02m[m+2] * ( 1/( rhobar ) ) * rho * ( ThornzGqlm * mu * rhobar + zGqlm * \
                ( PsiCDe2bar * rho + ( rhobar * ( PsiCDe2 + ( 4 * mu * rhobar + -2 * \
                ( 4 * gammabar + mu ) * rho ) ) + -2 * ( rhobar + rho ) * ( \
                taubar*taubar ) ) ) ) + 2 * Y12m[m+2] * ( 1/( rho ) ) * ( -2 * ThornzGqlm \
                * mu * rho * taubar + zGqlm * ( PsiCDe2 * ( pibar * rhobar + -1 * rho \
                * taubar ) + 2 * rho * taubar * ( -4 * mu * rhobar + ( 6 * gammabar * \
                rho + ( mu * rho + taubar * ( pi + taubar ) ) ) ) ) ) ) )

        Thornphnm+=( complex( 0,2 ) * Y12m[m+2] * a * m * rhobar * ( rho*rho ) * ( \
                ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ( ThornThornzGqlm * Y22m[m+2] * \
                ( -1 * mubar + mu ) * ( 1/( ( rhobar + -1 * rho ) ) ) * ( pi + -1 * \
                taubar ) + ( -2 * Y12m[m+2] * ( 1/( rhobar ) ) * ( 1/( ( -1 * rhobar + rho \
                ) ) ) * ( zGqlm * ( rhobar * ( -3 * mu * ( rhobar*rhobar*rhobar ) + ( \
                2 * ( 3 * gammabar + mu ) * ( rhobar*rhobar ) * rho + ( ( rho*rho ) * \
                ( PsiCDe2 + -6 * gammabar * rho ) + rhobar * rho * ( -1 * PsiCDe2bar \
                + mu * rho ) ) ) ) + -2 * ( -1 * rhobar + rho ) * ( ( rhobar + rho \
                )*( rhobar + rho ) ) * ( taubar*taubar ) ) + ( -1 * rhobar + rho ) * \
                ( ThornThornzGqlm * mu * rhobar + ThornzGqlm * ( PsiCDe2bar * rho + ( \
                rhobar * ( PsiCDe2 + ( 3 * mu * rhobar + -6 * gammabar * rho ) ) + -2 \
                * ( rhobar + rho ) * ( taubar*taubar ) ) ) ) ) + 0.5 * ThornzGqlm * ( \
                1/( rhobar ) ) * ( 1/( rho ) ) * ( complex( 0,2 ) * Y22m[m+2] * a * m * \
                rhobar * ( rho*rho ) * ( -1 * rhobar + rho ) * taubar + Y22m[m+2] * ( \
                taubar * ( PsiCDe2bar * rho + ( 2 * ( -1 * rhobar + rho ) * ( 4 * \
                gammabar * rho + -1 * mu * ( 2 * rhobar + rho ) ) + ( 2 * pi * rho * \
                taubar + -2 * rhobar * ( taubar*taubar ) ) ) ) + PsiCDe2 * ( -1 * rho \
                * taubar + rhobar * ( pibar + taubar ) ) ) ) ) ) )

        Thornphmm+= ( complex( 0,-1 ) * Y22m[m+2] * a * m * rhobar * ( \
                ThornThornzGqlm + 2 * ThornzGqlm * rhobar ) * rho + Y22m[m+2] * ( 1/( \
                rhobar ) ) * ( 1/( rho ) ) * ( ThornThornThornzGqlm * mu * rhobar + ( \
                2 * ThornzGqlm * rhobar * ( PsiCDe2bar * rho + ( rhobar * ( PsiCDe2 + \
                ( mu * rhobar + -4 * gammabar * rho ) ) + -2 * ( rhobar + rho ) * ( \
                taubar*taubar ) ) ) + ThornThornzGqlm * ( PsiCDe2bar * rho + ( rhobar \
                * ( PsiCDe2 + ( 2 * mu * rhobar + -4 * gammabar * rho ) ) + -2 * ( \
                rhobar + rho ) * ( taubar*taubar ) ) ) ) ) )

    # else:
    #     #q=0 
    #     for m in [-2,-1,0,1,2]:
    #         Thornphnndag += -2 * ( complex( 0,1 ) * zGqlm * a * m * rhobar * ( \
    #             rho*rho ) * ( sqrtof6 * Y02m[m+2] * rho + -2 * Y12m[m+2] * taubar ) + ( ( 3/2 \
    #             )**( 0.5 ) * Y02m[m+2] * ( 1/( rhobar ) ) * rho * ( -2 * ThornzGqlm * mu * \
    #             rhobar + zGqlm * ( PsiCDe2bar * ( -2 + q ) * rho + ( rhobar * ( \
    #             PsiCDe2 * ( -2 + q ) + ( 2 * ( -4 + q ) * mu * rhobar + ( -4 * ( -4 + \
    #             q ) * gammabar * rho + ( 4 + -2 * q ) * mu * rho ) ) ) + -2 * ( -2 + \
    #             q ) * ( rhobar + rho ) * ( taubar*taubar ) ) ) ) + Y12m[m+2] * ( 1/( rho ) \
    #             ) * ( 2 * ThornzGqlm * mu * rho * taubar + zGqlm * ( PsiCDe2 * ( -1 + \
    #             q ) * ( pibar * rhobar + -1 * rho * taubar ) + rho * taubar * ( ( 8 + \
    #             -2 * q ) * mu * rhobar + ( 4 * ( -3 + q ) * gammabar * rho + ( 2 * ( \
    #             -1 + q ) * mu * rho + 2 * ( -1 + q ) * taubar * ( pi + taubar ) ) ) ) \
    #             ) ) ) )

    #         Thornphnm+= ( -2 * rho * ( Y12m[m+2] * ( ( -4 + q ) * gammabar + ( -2 + q ) \
    #             * gamma ) + complex( 0,-1 ) * Y12m[m+2] * a * m * rhobar * rho ) * ( \
    #             ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ( ThornzGqlm * ( -1 * Y22m[m+2] \
    #             * ( -4 * gammabar + q * ( gammabar + gamma ) ) + complex( 0,1 ) * \
    #             Y22m[m+2] * a * m * rhobar * rho ) * ( pi + -1 * taubar ) + ( Y12m[m+2] * ( 1/( \
    #             rhobar ) ) * ( 1/( ( -1 * rhobar + rho ) ) ) * ( ( -1 * rhobar + rho \
    #             ) * ( -2 * ThornThornzGqlm * mu * rhobar + ThornzGqlm * ( PsiCDe2 * ( \
    #             -2 + q ) * rhobar + ( PsiCDe2bar * ( -2 + q ) * rho + ( -2 * mu * \
    #             rhobar * ( rhobar + 2 * rho ) + -2 * ( -2 + q ) * ( rhobar + rho ) * \
    #             ( taubar*taubar ) ) ) ) ) + zGqlm * ( PsiCDe2 * ( -2 + q ) * rhobar * \
    #             ( rho*rho ) + ( 2 * mu * rhobar * ( ( rhobar*rhobar*rhobar ) + ( \
    #             rhobar * ( rho*rho ) + -2 * ( rho*rho*rho ) ) ) + ( -2 + q ) * ( -1 * \
    #             PsiCDe2bar * ( rhobar*rhobar ) * rho + -2 * ( -1 * rhobar + rho ) * ( \
    #             ( rhobar + rho )*( rhobar + rho ) ) * ( taubar*taubar ) ) ) ) ) + \
    #             Y22m[m+2] * ( ThornzGqlm * ( 2 * mubar * taubar + -1 * mu * ( pi + taubar \
    #             ) ) + ( pi + -1 * taubar ) * ( ThornThornzGqlm * ( mubar + -1 * mu ) \
    #             * ( 1/( ( -1 * rhobar + rho ) ) ) + 0.5 * ThornzGqlm * ( -1 + q ) * ( \
    #             1/( rhobar ) ) * ( 1/( rho ) ) * ( PsiCDe2 * rhobar + ( PsiCDe2bar * \
    #             rho + 2 * ( rhobar + rho ) * taubar * tau ) ) ) ) ) ) )

    #         Thornphmm+= ( -1 * ( ThornThornzGqlm + 2 * ThornzGqlm * rhobar ) * ( -1 \
    #             * Y22m[m+2] * ( -4 * gammabar + q * ( gammabar + gamma ) ) + complex( 0,1 \
    #             ) * Y22m[m+2] * a * m * rhobar * rho ) + -0.5 * Y22m[m+2] * ( 1/( rhobar ) ) * \
    #             ( 1/( rho ) ) * ( -2 * ThornThornThornzGqlm * mu * rhobar + ( \
    #             ThornThornzGqlm * ( rhobar * ( PsiCDe2 * ( -2 + q ) + -4 * mu * \
    #             rhobar ) + ( PsiCDe2bar * ( -2 + q ) * rho + -2 * ( -2 + q ) * ( \
    #             rhobar + rho ) * ( taubar*taubar ) ) ) + 2 * ThornzGqlm * rhobar * ( \
    #             rhobar * ( PsiCDe2 * ( -2 + q ) + -2 * mu * rhobar ) + ( PsiCDe2bar * \
    #             ( -2 + q ) * rho + -2 * ( -2 + q ) * ( rhobar + rho ) * ( \
    #             taubar*taubar ) ) ) ) ) )

    hnn=np.conj(hnndag)
    hnmb=np.conj(hnm)
    hmbmb=np.conj(hmm)

    Thornhnn=np.conj(Thornhnndag)
    Thornhnmb=np.conj(Thornhnm)
    Thornhmbmb=np.conj(Thornhmm)

    Thornphnn=np.conj(Thornphnndag)
    Thornphnmb=np.conj(Thornphnm)
    Thornphmbmb=np.conj(Thornphmm)

    Ethhnn=np.conj(Ethphnndag)
    Ethhnmb=np.conj(Ethphnm)
    Ethhmbmb=np.conj(Ethphmm)

    Ethphnn=np.conj(Ethhnndag)
    Ethphnmb=np.conj(Ethhnm)
    Ethphmbmb=np.conj(Ethhmm)

    # PDhnn_m= 
    # PDhnmb_m= 
    # PDhmbmb_m= 



    dChr=np.zeros((4,4,4),np.complex128) # complex numbers are being cast as real
    # Some terms have explicit q-dependence
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
                lnis=[1,0,3,2] # apply eta matrix
                if ii>1: signi= -1 
                else: signi = 1
                if jj>1: signj= -1 
                else: signj = 1
                if kk>1: signk= -1 
                else: signk = 1
                dChrKN += signi*signj*signk*np.einsum('i,j,k',Tet[lnis[ii]],InvTet[lnis[jj]],InvTet[lnis[kk]]) * dChr[ii, jj, kk]

print((time.perf_counter()-start)/10000)

it1, it2, it3 = 0,0,0
for x in dChrKN: 
    for y in x: 
        for z in y:
            if np.iscomplex(np.real_if_close(z,imtolfactor)):
                print((it1,it2,it3),z,dChrKN[it1,it2,it3]-dChrKN[it1,it3,it2])
            it3+=1
            it3=it3%4
        it2+=1
        it2=it2%4
    it1+=1
    it1=it1%4
            

# print(np.real_if_close(dChrKN[0,0,0],100),np.real_if_close(dChrKN[0,0,1],100),np.real_if_close(dChrKN[0,0,2],100),np.real_if_close(dChrKN[0,0,3],100))
# print(np.real_if_close(dChrKN[0,1,0],100),np.real_if_close(dChrKN[0,1,1],100),np.real_if_close(dChrKN[0,1,2],100),np.real_if_close(dChrKN[0,1,3],100))
# print(np.real_if_close(dChrKN[0,2,0],100),np.real_if_close(dChrKN[0,2,1],100),np.real_if_close(dChrKN[0,2,2],100),np.real_if_close(dChrKN[0,2,3],100))
# print(np.real_if_close(dChrKN[0,3,0],100),np.real_if_close(dChrKN[0,3,1],100),np.real_if_close(dChrKN[0,3,2],100),np.real_if_close(dChrKN[0,3,3],100))
# dChrKN = np.real_if_close(dChrKN,1)
# print("dChrKN calculated ",time.perf_counter()-start,time.perf_counter()-restart)

