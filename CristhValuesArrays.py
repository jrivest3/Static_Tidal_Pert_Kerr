# import time

import numpy as np
import spherical
import quaternionic
# from scipy.integrate import solve_ivp as ODE



sqrtof2=np.sqrt(2)
sqrtof3=np.sqrt(3)
sqrtof5=np.sqrt(5)
sqrtof6=sqrtof2*sqrtof3

# start=time.perf_counter()
# restart=start

def Gamma(M,a,radius,theta): # Calculated by BHPToolkit
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

    #outgoing Kerr coords
    ChrKN=np.zeros((4,4,4),np.complex128)
    ChrKN[0,0,0]=-1/8 * M * ( 1/( rho ) ) * ( 1/( rhobar ) ) * ( ( \
    rho*rho ) + ( rhobar*rhobar ) ) * ( 2 * rho * rhobar + ( ( \
    rhobar*rhobar ) + ( rho*rho ) * ( 1 + 4 * ( a*a ) * ( rhobar*rhobar ) \
    ) ) )
    
    ChrKN[0,0,2]=( a*a ) * M * rho * rhobar * ( rho + rhobar ) * CosTH * SinTH
    ChrKN[0,0,3]=-1/4 * ( 1/( a ) ) * M * ( 1/( rho*rho*rho ) ) * ( 1/( \
    rhobar*rhobar*rhobar ) ) * ( ( rho*rho ) + ( rhobar*rhobar ) ) * ( 2 \
    * rho * rhobar + ( ( rhobar*rhobar ) + ( rho*rho ) * ( 1 + 4 * ( a*a \
    ) * ( rhobar*rhobar ) ) ) ) * ( taubar*taubar )
    
    ChrKN[0,1,2]= 1/sqrtof2 * ( 1/( rhobar ) - 1/( rho ) ) * taubar# 1/sqrtof2 * ( 1/( rho ) ) * ( rho - rhobar ) * ( 1/( rhobar ) ) * taubar
    ChrKN[0,1,3]=-0.5 * a * ( rho + rhobar ) * ( SinTH*SinTH )
    ChrKN[0,2,0]=ChrKN[0,0,2]
    ChrKN[0,2,1]=ChrKN[0,1,2]
    ChrKN[0,2,2]=-0.5 * ( Delta + 2 * M * radius ) * ( rho + rhobar )
    ChrKN[0,2,3]=- ( a*a*a ) * M * rho * rhobar * ( rho + rhobar ) * \
    CosTH * ( SinTH*SinTH*SinTH )
    ChrKN[0,3,0]=ChrKN[0,0,3]
    ChrKN[0,3,1]=ChrKN[0,1,3]
    ChrKN[0,3,2]=ChrKN[0,2,3]
    ChrKN[0,3,3]=1/16 * ( 1/( a*a ) ) * ( 1/( rho*rho*rho*rho ) ) * ( 1/( \
    rhobar*rhobar*rhobar*rhobar ) ) * ( ( a*a ) + 1/4 * ( 1/( rho*rho ) ) \
    * ( 1/( rhobar*rhobar ) ) * ( ( rho + rhobar )*( rho + rhobar ) ) ) * \
    ( -4 * ( gammabar - mubar ) * ( ( rho - rhobar )*( rho -1 \
    * rhobar )*( rho - rhobar )*( rho - rhobar ) ) + ( 8 * ( \
    a*a ) * M * ( rho*rho*rho ) * ( rhobar*rhobar*rhobar ) * ( ( rho + \
    rhobar )*( rho + rhobar ) ) + ( ( ( rho + rhobar ) )**( 5 ) + 2 * ( ( \
    rho - rhobar )*( rho - rhobar ) ) * ( 4 * ( a*a ) * M * ( \
    rho*rho*rho ) * ( rhobar*rhobar*rhobar ) - ( ( rho + rhobar )*( \
    rho + rhobar ) ) * ( rho + ( rhobar - M * rho * rhobar ) ) ) ) ) \
    ) * ( taubar*taubar )
    ChrKN[1,0,0]=1/4 * M * ( ( rho*rho ) + ( rhobar*rhobar ) ) * ( 2 * ( \
    a*a ) * rho * rhobar + ( 2 * M * ( rho + rhobar ) - radius * ( rho + \
    rhobar ) ) )
    ChrKN[1,0,1]=1/4 * M * ( 2 * rho * rhobar + ( ( rhobar*rhobar ) + ( ( \
    rho*rho ) * ( 1 -4 * ( a*a ) * ( rhobar*rhobar ) ) -8 * ( \
    taubar*taubar ) ) ) )
    
    ChrKN[1,0,3]=2 * ( 1/( a ) ) * mu * ( 1/( rho*rho*rho ) ) * ( 1/( \
    rhobar*rhobar*rhobar ) ) * ( ( pi*pi ) + ( taubar*taubar ) ) * \
    PsiCDe2bar
    ChrKN[1,1,0]=ChrKN[1,0,1]
    
    ChrKN[1,1,2]= -ChrKN[0,1,2]# 1/sqrtof2 * ( - ( 1/( rho ) ) + ( 1/( rhobar ) ) ) * tau # tau=-taubar
    ChrKN[1,1,3]=( 1/( a ) ) * ( 1/( rho*rho*rho ) ) * ( 1/( \
    rhobar*rhobar*rhobar ) ) * ( 4 * rhobar * ( gammabar * rho - mu \
    * rhobar ) * ( taubar*taubar ) + rho * ( ( pi + taubar )*( pi + \
    taubar ) ) * PsiCDe2bar )
    
    ChrKN[1,2,1]=ChrKN[1,1,2]
    ChrKN[1,2,2]=0.5 * Delta * ( rho + rhobar )
    
    ChrKN[1,3,0]=ChrKN[1,0,3]
    ChrKN[1,3,1]=ChrKN[1,1,3]
    
    ChrKN[1,3,3]=-1/8 * ( 1/( a*a ) ) * mu * ( 1/( rho )**( 6 ) ) * ( 1/( \
    rhobar )**( 5 ) ) * ( -4 * ( gammabar - mubar ) * ( ( rho - \
    rhobar )*( rho - rhobar )*( rho - rhobar )*( rho - \
    rhobar ) ) + ( 8 * ( a*a ) * M * ( rho*rho*rho ) * ( \
    rhobar*rhobar*rhobar ) * ( ( rho + rhobar )*( rho + rhobar ) ) + ( ( \
    ( rho + rhobar ) )**( 5 ) + 2 * ( ( rho - rhobar )*( rho - \
    rhobar ) ) * ( 4 * ( a*a ) * M * ( rho*rho*rho ) * ( \
    rhobar*rhobar*rhobar ) - ( ( rho + rhobar )*( rho + rhobar ) ) * \
    ( rho + ( rhobar - M * rho * rhobar ) ) ) ) ) ) * ( \
    taubar*taubar )
    ChrKN[2,0,0]=( a*a ) * M * ( rho*rho ) * ( rhobar*rhobar ) * ( rho + \
    rhobar ) * CosTH * SinTH
    
    ChrKN[2,0,3]=- a * M * ( Delta + 2 * M * radius ) * ( rho*rho ) * ( \
    rhobar*rhobar ) * ( rho + rhobar ) * CosTH * SinTH
    
    ChrKN[2,1,2]=-0.5 * ( rho + rhobar )
    ChrKN[2,1,3]=ChrKN[1,1,2]/a # 1/sqrtof2 * ( 1/( a ) ) * ( - ( 1/( rho ) ) + ( 1/( rhobar ) ) ) * tau
    
    ChrKN[2,2,1]=ChrKN[2,1,2]
    ChrKN[2,2,2]=ChrKN[1,1,2] # 1/sqrtof2 * ( - ( 1/( rho ) ) + ( 1/( rhobar ) ) ) * \tau
    
    ChrKN[2,3,0]=ChrKN[2,0,3]
    ChrKN[2,3,1]=ChrKN[2,1,3]
    
    ChrKN[2,3,3]=1/64 * 1/sqrtof2 * ( 1/( a*a ) ) * ( 1/( rho*rho*rho*rho \
    ) ) * ( - ( 1/( rho ) ) + ( 1/( rhobar ) ) ) * ( 1/( \
    rhobar*rhobar*rhobar*rhobar ) ) * tau * ( 8 * mu * ( ( rho - \
    rhobar )*( rho - rhobar )*( rho - rhobar )*( rho - \
    rhobar ) ) * rhobar + ( -8 * ( a*a ) * ( 2 * M + radius ) * ( rho*rho*rho \
    ) * ( rhobar*rhobar*rhobar ) * ( ( rho + rhobar )*( rho + rhobar )*( \
    rho + rhobar ) ) + ( ( ( rho + rhobar ) )**( 6 ) + ( -64 * ( a*a ) * \
    ( 2 * ( M*M ) * radius + Delta * ( M + radius ) ) * ( rho )**( 5 ) * ( rhobar \
    )**( 5 ) * ( rho + rhobar ) * ( CosTH*CosTH ) -64 * ( a*a ) * M * ( \
    Delta + 2 * M * radius ) * ( rho )**( 5 ) * ( rhobar )**( 5 ) * ( rho + \
    rhobar ) * ( SinTH*SinTH ) ) ) ) )
    ChrKN[3,0,0]=-0.5 * a * ( rhobar * PsiCDe2 + rho * PsiCDe2bar )
    
    ChrKN[3,0,2]=-2 * sqrtof2 * a * M * beta * rho * ( rho + rhobar )
    ChrKN[3,0,3]=- Sigma * ( 1/( rhobar ) ) * ( ( pi*pi ) + ( \
    taubar*taubar ) ) * PsiCDe2bar
    
    ChrKN[3,1,2]=-2 * sqrtof2 * a * beta * rho
    ChrKN[3,1,3]=ChrKN[2,1,2] # 0.5 * ( - rho - rhobar )
    ChrKN[3,2,0]=ChrKN[3,0,2]
    ChrKN[3,2,1]=ChrKN[3,1,2]
    ChrKN[3,2,2]=-0.5 * a * ( rho + rhobar )
    ChrKN[3,2,3]=-0.5 * 1/sqrtof2 * beta * ( rho*rho ) * rhobar * ( 3 * ( \
    a*a*a*a ) + ( -4 * ( a*a ) * M * ( 1/( rho ) ) * ( 1/( rhobar ) ) * ( \
    rho + rhobar ) + ( 2 * ( a*a ) * ( 1/( rho*rho ) ) * ( 1/( \
    rhobar*rhobar ) ) * ( ( rho + rhobar )*( rho + rhobar ) ) + ( 0.5 * ( \
    1/( rho*rho*rho*rho ) ) * ( 1/( rhobar*rhobar*rhobar*rhobar ) ) * ( ( \
    rho + rhobar )*( rho + rhobar )*( rho + rhobar )*( rho + rhobar ) ) + \
    ( 2 * ( 1/( rho*rho*rho ) ) * ( 1/( rhobar*rhobar*rhobar ) ) * ( 2 * \
    ( a*a ) * rho * rhobar - Deltap * ( rho + rhobar ) ) * ( ( a*a ) \
    * ( rho*rho ) * ( rhobar*rhobar ) + 4 * ( taubar*taubar ) ) + ( \
    a*a*a*a ) * np.cos( 4 * theta ) ) ) ) ) )
    ChrKN[3,3,0]=- Sigma * ( 1/( rhobar ) ) * ( ( pi*pi ) + ( \
    taubar*taubar ) ) * PsiCDe2bar
    ChrKN[3,3,1]=ChrKN[3,1,3]
    ChrKN[3,3,2]=ChrKN[3,2,3]
    ChrKN[3,3,3]=1/16 * ( 1/( a ) ) * ( 1/( rho*rho*rho*rho ) ) * ( 1/( \
    rhobar*rhobar*rhobar*rhobar ) ) * ( -4 * ( gammabar - mubar ) * \
    ( ( rho - rhobar )*( rho - rhobar )*( rho - rhobar )*( \
    rho - rhobar ) ) + ( 8 * ( a*a ) * M * ( rho*rho*rho ) * ( \
    rhobar*rhobar*rhobar ) * ( ( rho + rhobar )*( rho + rhobar ) ) + ( ( \
    ( rho + rhobar ) )**( 5 ) + 2 * ( ( rho - rhobar )*( rho - \
    rhobar ) ) * ( 4 * ( a*a ) * M * ( rho*rho*rho ) * ( \
    rhobar*rhobar*rhobar ) - ( ( rho + rhobar )*( rho + rhobar ) ) * \
    ( rho + ( rhobar - M * rho * rhobar ) ) ) ) ) ) * ( \
    taubar*taubar )
    #BL coords
    # Chr=np.zeros((4,4,4),np.complex128)
    
    # Chr[0,0,1]=-1/32 * M /mu /rhobar * ( 2 * rho * \
    #     rhobar + ( ( rhobar*rhobar ) + ( rho*rho ) * ( 1 + 4 * ( a*a ) * ( \
    #     rhobar*rhobar ) ) ) ) * ( ( rho + ( rhobar + 2 * a * rho * rhobar ) ) \
    #     * ( -rhobar + rho * ( -1 + 2 * a * rhobar ) ) + 8 * ( taubar*taubar ) )
    # Chr[0,0,2]=( a*a ) * M * rho * rhobar * ( rho + rhobar ) * CosTH * SinTH

    # Chr[0,1,0]=Chr[0,0,1]

    # Chr[0,1,3]=-1/16 /a * M /mu /( rho*rho ) \
    #     * ( 1/( rhobar*rhobar*rhobar ) ) * ( taubar*taubar ) * ( 8 * ( \
    #     a*a*a*a ) * ( rho*rho*rho*rho ) * ( rhobar*rhobar*rhobar*rhobar ) + ( \
    #     -6 * ( a*a ) * ( rho*rho ) * ( rhobar*rhobar ) * ( ( rho + rhobar )*( \
    #     rho + rhobar ) ) + ( -3 * ( ( rho + rhobar )*( rho + rhobar )*( rho + \
    #     rhobar )*( rho + rhobar ) ) + 4 * ( a + radius ) * rho * rhobar * ( rho + \
    #     ( rhobar + 2 * a * rho * rhobar ) ) * ( ( a*a ) * ( rho*rho ) * ( \
    #     rhobar*rhobar ) + 4 * ( taubar*taubar ) ) ) ) )
    # Chr[0,2,0]=Chr[0,0,2]

    # Chr[0,2,3]=-( a*a*a ) * M * rho * rhobar * ( rho + rhobar ) * CosTH * ( SinTH*SinTH*SinTH )

    # Chr[0,3,1]=Chr[0,1,3]
    # Chr[0,3,2]=Chr[0,2,3]

    # Chr[1,0,0]=0.5 * M * mu /rho * ( 2 * rho * rhobar + ( ( \
    #     rhobar*rhobar ) + ( ( rho*rho ) * ( 1 -4 * ( a*a ) * ( rhobar*rhobar ) ) -8 * ( taubar*taubar ) ) ) )

    # Chr[1,0,3]=2 /a * mu * ( 1/( rho*rho*rho ) ) * ( 1/( \
    #     rhobar*rhobar*rhobar ) ) * ( ( pi*pi ) + ( taubar*taubar ) ) * PsiCDe2bar

    # Chr[1,1,1]=0.5 * ( -( 1/( mu ) ) * ( 2 * gammabar + mu ) * rho + rhobar )
    # Chr[1,1,2]=1/sqrtof2 * ( ( 1/( rhobar ) ) - ( 1/( rho ) ) ) * tau

    # Chr[1,2,1]=Chr[1,1,2]
    # Chr[1,2,2]=0.5 * Delta * ( rho + rhobar )
 
    # Chr[1,3,0]=Chr[1,0,3]

    # Chr[1,3,3]=-2/( a*a ) * mu * ( Sigma**5/rho ) \
    #     * ( taubar*taubar ) * ( ( 2 * ( gammabar -mubar ) * ( ( rho \
    #     -rhobar )*( rho -rhobar ) ) -2 * ( ( a*a ) -M * radius \
    #     ) * ( rho*rho ) * ( rhobar*rhobar ) * ( rho + rhobar ) ) * ( \
    #     taubar*taubar ) + ( Delta + 2 * M * radius ) * ( rho*rho*rho*rho ) * ( \
    #     rhobar*rhobar*rhobar*rhobar ) * ( rho + rhobar ) * ( ( radius*radius ) + ( a*a ) * np.cos( 2 * theta ) ) )
    # Chr[2,0,0]=( a*a ) * M * ( rho*rho ) * ( rhobar*rhobar ) * ( rho + rhobar ) * CosTH * SinTH

    # Chr[2,0,3]=- a * M * ( ( a*a ) + ( radius*radius ) ) * ( rho*rho ) * ( rhobar*rhobar ) * ( rho + rhobar ) * CosTH * SinTH

    # Chr[2,1,1]=0.5 /sqrtof2 / mu   * rho * ( rho -rhobar ) * taubar
    # Chr[2,1,2]=-0.5 * ( rho + rhobar )

    # Chr[2,2,1]=Chr[2,1,2]
    # Chr[2,2,2]=1/sqrtof2 * ( ( 1/( rhobar ) ) - ( 1/( rho ) ) ) * tau

    # Chr[2,3,0]=Chr[2,0,3]

    # Chr[2,3,3]=- 1/sqrtof2 /( a*a ) /( rho*rho*rho )  * ( rho - rhobar ) /( rhobar*rhobar*rhobar )  * \
    #     ( -2* mu * rhobar + M * ( ( ( a*a ) + ( radius*radius ) )*( ( a*a ) + ( \
    #     radius*radius ) ) ) * ( rho*rho*rho ) * ( rhobar*rhobar*rhobar ) * ( rho + rhobar ) ) * tau

    # Chr[3,0,1]=1/8 * a * M /mu * ( rho*rho ) * rhobar * ( 2 * rho * \
    #     rhobar + ( ( rhobar*rhobar ) + ( ( rho*rho ) * ( 1 -4 * ( a*a ) * ( rhobar*rhobar ) ) -8 * ( taubar*taubar ) ) ) )
    # Chr[3,0,2]=-2 * sqrtof2 * a * M * beta * rho * ( rho + rhobar )

    # Chr[3,1,0]=Chr[3,0,1]

    # Chr[3,1,3]=1/32 /mu /rhobar * ( -( 2 * M -radius ) * \
    #     ( ( rho + rhobar )*( rho + rhobar )*( rho + rhobar )*( rho + rhobar ) \
    #     ) + ( -2 * ( a*a ) * rho * rhobar * ( ( rho + rhobar )*( rho + rhobar \
    #     ) ) * ( rho + ( rhobar + 3 * M * rho * rhobar ) ) + ( ( a*a*a*a ) * ( \
    #     rho*rho*rho ) * ( rhobar*rhobar*rhobar ) * ( 2 * M * rho * rhobar + \
    #     -3 * ( rho + rhobar ) ) + ( -4 * ( a*a ) * ( Delta + radius * ( M + radius ) ) \
    #     * ( rho*rho*rho ) * ( rhobar*rhobar*rhobar ) * ( rho + rhobar ) * \
    #     np.cos( 2 * theta ) + 4 * ( a*a*a*a ) * ( gammabar -mubar ) * ( \
    #     rho*rho*rho ) * ( rhobar*rhobar*rhobar ) * np.cos( 4 * theta ) ) ) ) )
    # Chr[3,2,0]=Chr[3,3,1]=Chr[3,0,2]

    # Chr[3,2,3]=-2 * sqrtof2 * beta /rhobar - a*a * M * rho * rhobar * ( rho + rhobar ) * CosTH * SinTH

    # Chr[3,3,1]=Chr[3,1,3]
    # Chr[3,3,2]=Chr[3,2,3]
    return ChrKN.real


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

def dGamma(M,a,radius,theta,phi_star,z_array):
    # if radius<=(M+np.sqrt(M*M-a*a)): return np.zeros((4,4,4)) # the test particle shouldn't cross this, this is for making 3D plots

    CosTH=np.cos(theta)
    SinTH=np.sin(theta)
    
    rho=(-1/(radius - 1.j*a*CosTH))
    rhobar=np.conj(rho) # (-1/(radius + 1.j*a*CosTH))

    PsiCDe2=M*rho*rho*rho
    PsiCDe2bar=M*rhobar*rhobar*rhobar

    Delta= a*a - 2*M*radius + radius*radius
    Deltap= 2*(radius-M)
    Deltapp= 2
    Sigma= 1/(rho*rhobar)

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
    gamma= mu + rho*rhobar*Deltap/4 # Computation time may be shortened by a simplification that removes gamma and gammabar
    gammabar= gamma-mu+mubar
    # alpha= pi-betabar
    # alphabar= pibar-beta

    # Kinnersley Tetrad in BL coordinates # Outgoing Kerr-Newman Coordinates
    lBL=np.array([(radius*radius+a*a)/Delta,1,0,a/Delta]) # lKN=np.array([0,1,0,0])
    nBL=np.array([(radius*radius+a*a),-Delta,0,a])/Sigma/2 # nKN=np.array([(radius*radius+a*a),-Delta/2,0,a])*rho*rhobar# decide if Delta is divided by 2 or not
    mBL=-rhobar/sqrtof2*np.array([1.j*a*SinTH,0,1,1.j/SinTH])
    mbBL=np.conj(mBL)#-rho/sqrtof2*np.array([-1.j*a*SinTH,0,1,-1.j/SinTH])

    #lower with g
    lBLd=np.array([-1, Sigma/Delta, 0, a* SinTH*SinTH]) # lKNd=np.array([-1, 0, 0, a* SinTH*SinTH]) 
    nBLd=np.array([-Delta/Sigma, -1, 0, Delta/Sigma * a*SinTH*SinTH])/2 # nKNd=np.array([(-Delta*rho*rhobar)/2, -1, 0, (1/2)*Delta*rho*rhobar * a*SinTH*SinTH])
    mBLd=np.array([1.j*rhobar* a*SinTH, 0, -1/(rho), -1.j*rhobar*(radius*radius + a*a)*SinTH])/sqrtof2
    mbBLd=np.conj(mBLd) # np.transpose(np.array([ 1.j*rho* a*SinTH/sqrtof2, 0, 1/(sqrtof2*rhobar), -1.j*rho*(radius*radius + a*a)*SinTH/sqrtof2]))

    Tet=np.array([lBL,nBL,mBL,mbBL]).transpose() # np.array([lKN,nKN,mKN,mbKN])
    # raise with eta.
    InvTet=np.array([-nBLd,-lBLd,mbBLd,mBLd]) # np.array([-nKNd,-lKNd,mbKNd,mKNd]) # negate each (-1*) for (+---), but doesn't change end result?
    # print(np.einsum('ij,jk',InvTet,Tet)) # these are correct
    
    sY2m=swsh_array(ell, theta, phi_star).reshape((2*ell+1,2*ell+1)) # Are these Goldberg of Penrose and Rindler?
    
    # The z's may be complex, but have a reality condition
    # z[-m]=(-1)**m * np.conj(z[m]),
    # while Y[-s]l[-m] = (-1)**(m+s) * np.conj(Y[s]l[m]) so the conjugate relation is maintained.
    # So, np.conj(z[m]*Y[s]l[m])=(-1)**(s) * z[-m] * Y[-s]l[-m]
    # The z's are taken to be constants, and every instance of |s|=<2 is the result of a derivative.
    # s=0 can be reached by derivatives of both s=+/-1, but we only take derivatives s=2 and then conjugate the results.
    #Yn22m= np.conj(z_array)*sY2m[:,0] # s=-2 # not currently used in the calculation
    Yn12m= z_array*sY2m[:,1] # s=-1 # This only appears in Eth'hnnbar where three spin-lowering's occur
    Y02m = z_array*sY2m[:,2] # s= 0
    Y12m = z_array*sY2m[:,3] # s= 1
    Y22m = z_array*sY2m[:,4] # s= 2

    # Not Dependent on q
    # q=0
    # ell=2

    # for ell>2, these need to be either arrays or functions of m
    Gqlm= Delta*Delta/12
    # my h-components are in a form that only contain zG, ThornzG, and ThornThornzG.
    # dChr will contain first derivatives of the h-comp's which are also in terms of Thorn's only
    ThornGqlm= Deltap*Delta/6
    ThornThornGqlm= (Deltapp*Delta+Deltap*Deltap)/6
    ThornThornThornGqlm= Deltap#(*Deltapp/2)

    # Pieces of dGamma

    # hnn_2m= -Delta*Delta*(sqrtof6*rhobar*rhobar*Y02m+2*sqrtof2*rhobar*tau*Yn12m)
    # hnmb_2m= sqrtof2*rhobar*(rho+rhobar)*Delta*Delta*Yn12m + 2*Deltap*Delta*(sqrtof2*rhobar*Yn12m + (tau-pibar)*Yn22m)
    # hmbmb_2m= -(2*Deltap^2+Delta*(4+4*rho*Deltap))*Yn22m



    # At some point, we must sum over m
    # if ell==2
    hnndag=np.sum(Gqlm * rho * (2 * sqrtof2 * Y12m * taubar - sqrtof6 * Y02m * rho)
                    )
    hnm=np.sum(( sqrtof2 * Y12m * rho * ( ThornGqlm + Gqlm * ( rhobar + rho ) ) + ThornGqlm * Y22m * ( pi - taubar ) )
                )
    hmm=np.sum(- Y22m * ( ThornThornGqlm + 2 * ThornGqlm * rhobar )
                )
    # hnndag=np.sum(2 * Gqlm * ( 2 * Y12m * pi * rhobar - sqrtof6 * Y02m * ( rho*rho ) )
    #     )
    # hnm=np.sum(( 2 * Y12m * rho * ( ThornGqlm + Gqlm * ( rhobar + rho ) ) + ThornGqlm * Y22m * ( pi - taubar ) )
    #     )
    # hmm=np.sum(0 - Y22m * ( ThornThornGqlm + 2 * ThornGqlm * rhobar )
    #     )
    hnn=np.conj(hnndag)
    hnmb=np.conj(hnm)
    hmbmb=np.conj(hmm)
    # hnn=np.sum(- Gqlm * ( sqrtof6 * ( rhobar*rhobar ) * Y02m + 2*sqrtof2 * tau * rhobar * Yn12m)
    #     )
    # hnmb=np.sum(-0.5 * ( 2 * sqrtof2 * Gqlm * Yn12m * rhobar * ( rhobar + \
        # rho ) + 2 * ThornGqlm * ( sqrtof2 * Yn12m * rhobar + Yn22m * ( - pibar + tau ) ) )
    #     )
    # hmbmb=np.sum(0 - Y22m * ( ThornThornGqlm + 2 * ThornGqlm * rhobar )
    #     )
    
    # hnndag=np.conj(hnn)
    # hnm=np.conj(hnmb)
    # hmm=np.conj(hmbmb)

    # Each Thornp of h contains explicit m-dependence, so I initialize them to be summed by a loop.
    Thornphnndag=0
    Thornphnm=0
    Thornphmm=0

    # Thornhnndag=0
    # Thornhnm=0
    # Thornhmm=0
    # Ethhnndag=0
    # Ethhnm=0
    # Ethhmm=0
    # Ethphnndag=0
    # Ethphnm=0
    # Ethphmm=0
    Thornhnndag=np.sum(rho * ( - sqrtof6 * Y02m * rho * ( ThornGqlm + 2 * Gqlm \
    * rho ) + 2 * sqrtof2 * Y12m * ( ThornGqlm + Gqlm * ( rhobar + 2 * \
    rho ) ) * taubar )
    )

    Ethhnndag=np.sum(Gqlm * ( 2 * rho * ( 2 * Y22m * rhobar + sqrtof6 * Y02m * ( \
    rhobar + rho ) ) * taubar + sqrtof2 * Y12m * ( PsiCDe2bar * rho + ( \
    rhobar * ( 3 * PsiCDe2 + ( -2 * mu * rhobar + ( 2 * mu * rho -3 * ( \
    rho*rho ) ) ) ) -2 * ( 2 * rhobar + rho ) * ( taubar*taubar ) ) ) )
    )

    Ethphnndag=np.sum(sqrtof2 * Gqlm * rho * ( 3 * Yn12m * ( rho*rho ) + 2 * \
    taubar * ( - sqrtof3 * Y02m * rho + Y12m * taubar ) )
    )
    
    Thornhnm=np.sum( sqrtof2 * Y12m * rho * ( ThornThornGqlm + ( ThornGqlm * \
    ( rhobar + 2 * rho ) + Gqlm * ( ( rhobar*rhobar ) + ( rhobar * rho + \
    2 * ( rho*rho ) ) ) ) ) + Y22m * ( ThornThornGqlm * ( pi - taubar ) \
    + ThornGqlm * ( 2 * pi * rho - ( rhobar + rho ) * taubar ) ) 
    )
    
    Ethhnm=np.sum( 2 * Y22m * rhobar * rho * ( ThornGqlm + Gqlm * ( rhobar + \
    rho ) ) + ( - sqrtof2 * Y12m * ( ThornGqlm * ( 2 * rhobar + rho ) + \
    Gqlm * ( rhobar + rho ) * ( rhobar + 2 * rho ) ) * taubar + 0.5 * \
    ThornGqlm * Y22m * ( 1/( rhobar ) ) * ( 1/( rho ) ) * ( 2 * mu * \
    rhobar * ( ( - rhobar + rho )*( - rhobar + rho ) ) + ( PsiCDe2bar * \
    rho * ( 3 * rhobar + rho ) + ( PsiCDe2 * rhobar * ( rhobar + 3 * rho \
    ) -2 * ( -2 * ( rhobar*rhobar ) + ( 3 * rhobar * rho + ( rho*rho ) ) \
    ) * ( taubar*taubar ) ) ) ) ) 
    )
    
    Ethphnm= -1/( rhobar ) * np.sum( sqrtof6 * Y02m * rhobar * ( rho*rho ) \
    * ( ThornGqlm + Gqlm * ( rhobar + rho ) ) + ( sqrtof2 * Y12m * rho \
    * ( - rhobar + 2 * rho ) * ( ThornGqlm + Gqlm * ( rhobar + rho ) ) \
    * taubar + ThornGqlm * Y22m * taubar * ( 2 * pi * rho + ( rhobar * \
    taubar - rho * taubar ) ) ) 
    )
    
    Thornhmm=- np.sum(Y22m * ( ThornThornThornGqlm + 2 * rhobar * ( \
    ThornThornGqlm + ThornGqlm * rhobar ) )
    )
    
    Ethhmm=np.sum(2 * Y22m * rhobar * ( ThornThornGqlm + ThornGqlm * rhobar ) * taubar / rho 
    )
    
    # Ethphmm=( sqrtof2 * Y12m * ( ThornThornGqlm + 2 * ThornGqlm * \
    # rhobar ) * rho + 2 * Y22m * ( 1/( rhobar ) ) * ( ThornThornGqlm * \
    # rho + ThornGqlm * rhobar * ( - rhobar + 2 * rho ) ) * taubar 
    # )
    Ethphmm=np.sum( sqrtof2 * Y12m * ( ThornThornGqlm + 2 * ThornGqlm * rhobar \
    ) * rho + 2 * Y22m * ( ThornThornGqlm * pi + ThornGqlm * ( - rhobar + \
    2 * rho ) * taubar ) )

    for m in [-2,-1,0,1,2]:
        ### Attempt in BL
        # Thornhnndag+=- sqrtof2 * rho * ( ThornGqlm * ( sqrtof3 * Y02m[m+2]* rho \
        #     -2 * Y12m[m+2]* taubar ) + Gqlm * ( sqrtof3 * Y02m[m+2]* rho * ( complex( 0,1 \
        #     ) * ( 1/( Delta ) ) * a * m + 2 * rho ) -2 * Y12m[m+2]* ( complex( 0,1 ) \
        #     * ( 1/( Delta ) ) * a * m + ( rhobar + 2 * rho ) ) * taubar ) )
        # Thornphnndag+=( sqrtof2 * ThornGqlm * mu * ( sqrtof3 * Y02m[m+2]* rho -2 * \
        #     Y12m[m+2]* taubar ) + Gqlm * ( sqrtof6 * Y02m[m+2]* ( 1/( rhobar ) ) * rho * \
        #     ( PsiCDe2 * rhobar + ( PsiCDe2bar * rho + rhobar * ( 4 * mu * rhobar \
        #     + ( -2 * ( 4 * gammabar + mu ) * rho + ( complex( 0,-0.5 ) * a * m * \
        #     rhobar * ( rho*rho ) -2 * taubar * ( pi + taubar ) ) ) ) ) ) + \
        #     sqrtof2 * Y12m[m+2]* ( 1/( rho ) ) * ( PsiCDe2 * ( pibar * rhobar - rho * \
        #     taubar ) + 2 * rho * taubar * ( -4 * mu * rhobar + ( rho * ( 6 * \
        #     gammabar + ( mu + complex( 0,0.5 ) * a * m * rhobar * rho ) ) + \
        #     taubar * ( pi + taubar ) ) ) ) ) )
        # Ethhnndag+=Gqlm * ( 2 * rho * ( 2 * Y22m[m+2]* rhobar + sqrtof6 * Y02m[m+2]* ( \
        #     rhobar + rho ) ) * taubar + sqrtof2 * Y12m[m+2]* ( 3 * PsiCDe2 * rhobar + \
        #     ( -2 * mu * ( rhobar*rhobar ) + ( PsiCDe2bar * rho + ( 2 * mu * \
        #     rhobar * rho + ( -3 * rhobar * ( rho*rho ) + ( -4 * rhobar * ( \
        #     taubar*taubar ) -2 * rho * ( taubar*taubar ) ) ) ) ) ) ) )
        # Ethphnndag+=sqrtof2 * Gqlm * rho * ( 3 * ( Yn12m[m+2] ) * ( rho*rho ) + \
        #     2 * taubar * ( - sqrtof3 * Y02m[m+2]* rho + Y12m[m+2]* taubar ) )
        # Thornhnm+=( 1/( Delta*Delta ) ) * ( Gqlm * ( sqrtof2 * Y12m[m+2]* rho * ( \
        #     complex( 0,-1 ) * Deltap * a * m + ( - ( a*a ) * ( m*m ) + ( complex( \
        #     0,1 ) * Delta * a * m * ( rhobar + 2 * rho ) + ( Delta*Delta ) * ( ( \
        #     rhobar*rhobar ) + ( rhobar * rho + 2 * ( rho*rho ) ) ) ) ) ) + \
        #     complex( 0,-1 ) * Y22m[m+2]* a * m * ( -2 * Delta * pi * rho + ( Deltap * \
        #     ( pi - taubar ) + ( complex( 0,-1 ) * a * m * ( pi - taubar ) + Delta \
        #     * ( rhobar + rho ) * taubar ) ) ) ) + Delta * ( sqrtof2 * Y12m[m+2]* rho \
        #     * ( Delta * ThornThornGqlm + ThornGqlm * ( complex( 0,2 ) * a * m + \
        #     Delta * ( rhobar + 2 * rho ) ) ) + Y22m[m+2]* ( Delta * ThornThornGqlm * \
        #     ( pi - taubar ) + ThornGqlm * ( 2 * Delta * pi * rho + ( complex( 0,2 \
        #     ) * a * m * ( pi - taubar ) - Delta * ( rhobar + rho ) * taubar ) ) ) \
        #     ) )
        # Thornphnm+=0.5 * ( -2 * ( 1/( Delta ) ) * ThornGqlm * ( - mubar + mu ) \
        #     * ( 1/( ( - rhobar + rho ) ) ) * ( sqrtof2 * Y12m[m+2]* rho * ( complex( \
        #     0,1 ) * a * m + Delta * ( rhobar + rho ) ) + complex( 0,1 ) * Y22m[m+2]* \
        #     a * m * ( pi - taubar ) ) + ( - ( 1/( Delta ) ) * Gqlm * ( 1/( rhobar \
        #     ) ) * ( complex( 0,-0.5 ) * Y22m[m+2]* a * m * rhobar * ( 4 * ( \
        #     betabar*betabar ) * ( 2 * gammabar + ( - mubar + mu ) ) + ( -6 * \
        #     mubar * pi + ( -4 * mu * pi + ( -4 * betabar * ( 2 * gammabar + ( - \
        #     mubar + mu ) ) * pi + ( -2 * mubar * ( pi*pi ) + ( 2 * mu * ( pi*pi ) \
        #     + ( 4 * gammabar * pi * ( 2 + pi ) + ( complex( 0,1 ) * a * m * pi * \
        #     rhobar * rho + ( 4 * mubar * taubar + 2 * mu * taubar ) ) ) ) ) ) ) ) \
        #     ) + ( sqrtof2 * Y12m[m+2]* ( 2 * Delta * PsiCDe2 * rhobar * rho + ( \
        #     complex( 0,2 ) * a * m * ( mubar + mu ) * rhobar * rho + ( 2 * Delta \
        #     * PsiCDe2bar * rho * ( rhobar + rho ) + ( -12 * Delta * gammabar * \
        #     rhobar * rho * ( rhobar + rho ) + ( complex( 0,-1 ) * Delta * a * m * \
        #     ( rhobar*rhobar ) * ( rho*rho ) * ( rhobar + rho ) + ( 2 * Delta * ( \
        #     rhobar*rhobar ) * ( PsiCDe2 + mu * ( 3 * rhobar + rho ) ) -4 * Delta \
        #     * ( ( rhobar + rho )*( rhobar + rho ) ) * ( taubar*taubar ) ) ) ) ) ) \
        #     ) + a * m * ( Y22m[m+2]* ( 0.5 * a * m * rhobar * rho * ( -2 * rhobar + \
        #     rho ) + complex( 0,1 ) * ( -4 * mubar * rhobar + ( -4 * gammabar * ( \
        #     -2 * rhobar + rho ) + mu * ( rhobar + rho ) ) ) ) * taubar + complex( \
        #     0,1 ) * sqrtof2 * Y12m[m+2]* ( PsiCDe2bar * rho + ( complex( 0,-1 ) * a * \
        #     m * ( rhobar*rhobar ) * ( rho*rho ) + ( rhobar * ( PsiCDe2 + ( 4 * mu \
        #     * rhobar -2 * ( 6 * gammabar + mu ) * rho ) ) -2 * ( rhobar + rho ) * \
        #     ( taubar*taubar ) ) ) ) ) ) ) + ( - ( 1/( rhobar ) ) * ( 1/( rho ) ) \
        #     * ( sqrtof2 * Y12m[m+2]* rho + Y22m[m+2]* ( pi - taubar ) ) * ( 2 * \
        #     ThornThornGqlm * mu * rhobar + ThornGqlm * ( PsiCDe2 * rhobar + ( \
        #     PsiCDe2bar * rho + 2 * ( rhobar + rho ) * taubar * tau ) ) ) + \
        #     ThornGqlm * ( 2 * Y22m[m+2]* ( ( 4 * gammabar + complex( 0,0.5 ) * a * m \
        #     * rhobar * rho ) * ( pi - taubar ) + ( 2 * mubar * taubar - mu * ( pi \
        #     + taubar ) ) ) - sqrtof2 * Y12m[m+2]* ( 1/( rhobar ) ) * ( PsiCDe2 * \
        #     rhobar + ( PsiCDe2bar * rho + 2 * rhobar * ( rho * ( -4 * gammabar + \
        #     ( -2 * gamma + ( mu + complex( 0,-0.5 ) * a * m * rhobar * rho ) ) ) \
        #     + ( pi + taubar ) * tau ) ) ) ) ) ) )
        # Ethhnm+=1/4 * ( 1/( Delta ) ) * ( 1/( rhobar ) ) * ( 1/( rho ) ) * ( \
        #     -2 * Delta * ThornGqlm * pibar * rhobar * rho * ( 2 * sqrtof2 * Y12m[m+2]\
        #     * rho + 2 * Y22m[m+2]* ( pi - taubar ) ) + ( Gqlm * ( 8 * Delta * Y22m[m+2]* \
        #     ( rhobar*rhobar ) * ( rho*rho ) * ( rhobar + rho ) + ( complex( 0,1 ) \
        #     * Y22m[m+2]* a * m * ( rho * ( PsiCDe2bar * rho + rhobar * ( 3 * PsiCDe2 \
        #     + ( 4 * ( betabar*betabar ) * pibar + ( -2 * mu * rhobar + ( 2 * mu * \
        #     rho + 2 * pi * ( 2 * beta + ( pibar * ( -2 + ( -2 * betabar + pi ) ) \
        #     - taubar ) ) ) ) ) ) ) + ( -2 * pibar * rhobar * rho * ( pi -2 * \
        #     taubar ) -8 * pi * ( rhobar*rhobar ) * taubar ) ) + ( -4 * sqrtof2 * \
        #     Y12m[m+2]* rhobar * rho * ( complex( 0,1 ) * a * m * pibar * rho + Delta \
        #     * ( rhobar + rho ) * ( rhobar + 2 * rho ) * taubar ) + complex( 0,1 ) \
        #     * a * m * ( 4 * rhobar * rho * ( 2 * Y22m[m+2]* rhobar * rho - sqrtof2 * \
        #     Y12m[m+2]* ( 3 * rhobar + rho ) * taubar ) + Y22m[m+2]* ( 2 * mu * rhobar * ( \
        #     -2 * rhobar + rho ) * ( - rhobar + rho ) + ( PsiCDe2bar * rho * ( 6 * \
        #     rhobar + rho ) + ( PsiCDe2 * rhobar * ( 2 * rhobar + 3 * rho ) -2 * ( \
        #     -6 * ( rhobar*rhobar ) + ( 4 * rhobar * rho + ( rho*rho ) ) ) * ( \
        #     taubar*taubar ) ) ) ) ) ) ) ) + 2 * Delta * ThornGqlm * ( Y22m[m+2]* ( \
        #     PsiCDe2bar * rho * ( 3 * rhobar + rho ) + rhobar * ( 2 * mu * ( ( - \
        #     rhobar + rho )*( - rhobar + rho ) ) + ( PsiCDe2 * ( rhobar + 3 * rho \
        #     ) -2 * taubar * ( pi * rho + ( -3 * rhobar * taubar + 4 * rho * \
        #     taubar ) ) ) ) ) + 2 * rhobar * rho * ( 2 * Y22m[m+2]* rhobar * rho + \
        #     sqrtof2 * Y12m[m+2]* ( 3 * rhobar + rho ) * tau ) ) ) )
        # Ethphnm+=-0.5 * ( 1/( Delta ) ) * ( 1/( rhobar ) ) * ( 2 * Delta * \
        #     ThornGqlm * ( sqrtof6 * Y02m[m+2]* rhobar * ( rho*rho ) + ( sqrtof2 * \
        #     Y12m[m+2]* rho * ( - rhobar + 2 * rho ) * taubar + Y22m[m+2]* taubar * ( 2 * \
        #     pi * rho + ( rhobar - rho ) * taubar ) ) ) + Gqlm * ( complex( 0,-1 ) \
        #     * Y22m[m+2]* a * m * ( 2 * ( betabar*betabar ) * pi * rhobar + ( - ( \
        #     pi*pi ) * rhobar + ( ( pi*pi*pi ) * rhobar + ( -2 * betabar * pi * ( \
        #     1 + pi ) * rhobar + ( 2 * pi * ( rhobar - rho ) * taubar -2 * rhobar \
        #     * ( taubar*taubar ) ) ) ) ) ) + sqrtof2 * rho * ( 2 * sqrtof3 * Y02m[m+2]\
        #     * rhobar * rho * ( complex( 0,1 ) * a * m + Delta * ( rhobar + rho ) \
        #     ) + Y12m[m+2]* ( 2 * Delta * ( - ( rhobar*rhobar ) + ( rhobar * rho + 2 * \
        #     ( rho*rho ) ) ) * taubar + complex( 0,1 ) * a * m * ( 3 * pi * rhobar \
        #     + ( -2 * rhobar * taubar + rho * taubar ) ) ) ) ) )
        # Thornhmm+=- ( 1/( Delta*Delta*Delta ) ) * Y22m[m+2]* ( complex( 0,-1 ) * \
        #     Gqlm * a * m * ( -2 * ( Deltap*Deltap ) + ( complex( 0,3 ) * Deltap * \
        #     a * m + ( ( a*a ) * ( m*m ) + ( -2 * ( Delta*Delta ) * ( \
        #     rhobar*rhobar ) + 2 * Delta * ( 1 + ( Deltap * rhobar + complex( 0,-1 \
        #     ) * a * m * rhobar ) ) ) ) ) ) + Delta * ( ThornGqlm * ( complex( \
        #     0,-3 ) * Deltap * a * m + ( -3 * ( a*a ) * ( m*m ) + ( complex( 0,4 ) \
        #     * Delta * a * m * rhobar + 2 * ( Delta*Delta ) * ( rhobar*rhobar ) ) \
        #     ) ) + Delta * ( Delta * ThornThornThornGqlm + ThornThornGqlm * ( \
        #     complex( 0,3 ) * a * m + 2 * Delta * rhobar ) ) ) )
        # Thornphmm+=Y22m[m+2]* ( ThornThornThornGqlm * ( - mubar + mu ) * ( 1/( ( - \
        #     rhobar + rho ) ) ) + ( complex( 0,1 ) * ( 1/( Delta*Delta ) ) * \
        #     ThornGqlm * a * m * ( - mubar + mu ) * ( - Deltap + ( complex( 0,1 ) \
        #     * a * m + 2 * Delta * rhobar ) ) * ( 1/( ( - rhobar + rho ) ) ) + ( \
        #     0.5 * ( 1/( Delta*Delta ) ) * ( Gqlm * a * m * ( complex( 0,1 ) * \
        #     Deltap + ( a * m + complex( 0,-2 ) * Delta * rhobar ) ) - Delta * ( \
        #     Delta * ThornThornGqlm + 2 * ThornGqlm * ( complex( 0,1 ) * a * m + \
        #     Delta * rhobar ) ) ) * ( 8 * gammabar + complex( 0,1 ) * a * m * \
        #     rhobar * rho ) + ( ( 1/( rho ) ) * ( 2 * ThornThornGqlm * mu * rhobar \
        #     + 2 * ThornGqlm * ( rhobar * ( PsiCDe2 + mu * rhobar ) + ( PsiCDe2bar \
        #     * rho -2 * ( rhobar + rho ) * ( taubar*taubar ) ) ) ) + ( \
        #     ThornThornGqlm * ( 1/( rhobar ) ) * ( 1/( rho ) ) * ( PsiCDe2bar * \
        #     rho + rhobar * ( PsiCDe2 + 2 * ( pi + taubar ) * tau ) ) + ( complex( \
        #     0,-1 ) * ( 1/( Delta*Delta ) ) * Gqlm * a * m * ( 1/( rho ) ) * ( 2 * \
        #     ( Deltap + complex( 0,-1 ) * a * m ) * ( mubar + mu ) * rho - Delta * \
        #     ( PsiCDe2 * rhobar + ( PsiCDe2bar * rho + rho * ( 4 * mubar * rhobar \
        #     + ( 2 * mu * rhobar + ( - rhobar * rho + ( 2 * pibar * taubar + 2 * \
        #     taubar * tau ) ) ) ) ) ) ) + complex( 0,1 ) * ( 1/( Delta ) ) * a * m \
        #     * ( 1/( rhobar ) ) * ( 1/( rho ) ) * ( 2 * ThornThornGqlm * mu * \
        #     rhobar + ThornGqlm * ( PsiCDe2 * rhobar + ( PsiCDe2bar * rho + 2 * ( \
        #     mubar * rhobar * rho + ( mu * rhobar * rho + ( rhobar + rho ) * \
        #     taubar * tau ) ) ) ) ) ) ) ) ) ) )
        # Ethhmm+=2 * ( 1/( Delta*Delta ) ) * Y22m[m+2]* ( 1/( rho ) ) * ( - Gqlm * \
        #     a * m * ( complex( 0,1 ) * Deltap + ( a * m + complex( 0,-1 ) * Delta \
        #     * rhobar ) ) * ( pibar * rho + 2 * rhobar * taubar ) + Delta * ( \
        #     Delta * ThornThornGqlm * rhobar * taubar + ThornGqlm * ( Delta * ( \
        #     rhobar*rhobar ) * taubar + complex( 0,1 ) * a * m * ( pibar * rho + 3 \
        #     * rhobar * taubar ) ) ) )
        # Ethphmm+=( 1/( Delta*Delta ) ) * ( complex( 0,-1 ) * Gqlm * a * m * ( \
        #     sqrtof2 * Y12m[m+2]* ( Deltap + ( complex( 0,-1 ) * a * m -2 * Delta * \
        #     rhobar ) ) * rho + 2 * Y22m[m+2]* ( Deltap * pi + ( complex( 0,-1 ) * a * \
        #     m * pi + ( - Delta * pi * rhobar + ( Delta * rhobar * taubar - Delta \
        #     * rho * taubar ) ) ) ) ) + Delta * ( sqrtof2 * Y12m[m+2]* ( Delta * \
        #     ThornThornGqlm + 2 * ThornGqlm * ( complex( 0,1 ) * a * m + Delta * \
        #     rhobar ) ) * rho + 2 * Y22m[m+2]* ( Delta * ThornThornGqlm * pi + \
        #     ThornGqlm * ( complex( 0,2 ) * a * m * pi + ( - Delta * rhobar * \
        #     taubar + 2 * Delta * rho * taubar ) ) ) ) )

        Thornphnndag+=( complex( 0,-1 ) * Gqlm * a * m * rhobar * ( rho*rho ) \
    * ( sqrtof6 * Y02m[m+2] * rho -2 * sqrtof2 * Y12m[m+2] * taubar ) + ( sqrtof2 * \
    Y12m[m+2] * ( 1/( rho ) ) * ( -2 * ThornGqlm * mu * rho * taubar + Gqlm \
    * ( PsiCDe2 * ( pibar * rhobar - rho * taubar ) + 2 * rho * taubar * \
    ( -4 * mu * rhobar + ( 6 * gammabar * rho + ( mu * rho + taubar * ( \
    pi + taubar ) ) ) ) ) ) + sqrtof6 * Y02m[m+2] * rho * ( ThornGqlm * mu + \
    Gqlm * ( PsiCDe2 + ( 2 * ( -2 * ( gammabar + gamma ) + mu ) * rho + \
    ( PsiCDe2bar * ( 1/( rhobar ) ) * rho + 2 * ( pi + taubar ) * tau ) ) \
    ) ) ) )
        
        dummythornphnm=( ThornGqlm * ( 4  * gammabar + complex( 0,1 ) * \
    a * m * rhobar * rho ) * ( pi - taubar ) * Y22m[m+2] + ( sqrtof2 * ( \
    complex( 0,1 ) * Y12m[m+2] * a * m * rhobar * ( rho*rho ) * ( ThornGqlm + \
    Gqlm * ( rhobar + rho ) ) - Y12m[m+2] * ( 1/( rhobar ) ) * ( 1/( ( - \
    rhobar + rho ) ) ) * ( Gqlm * ( rhobar * ( -3 * mu * ( \
    rhobar*rhobar*rhobar ) + ( 2 * ( 3 * gammabar + mu ) * ( \
    rhobar*rhobar ) * rho + ( ( rho*rho ) * ( PsiCDe2 -6 * gammabar * rho \
    ) + rhobar * rho * ( - PsiCDe2bar + mu * rho ) ) ) ) -2 * ( - rhobar \
    + rho ) * ( ( rhobar + rho )*( rhobar + rho ) ) * ( taubar*taubar ) ) \
    + ( - rhobar + rho ) * ( ThornThornGqlm * mu * rhobar + ThornGqlm * \
    ( PsiCDe2bar * rho + ( rhobar * ( PsiCDe2 + ( 3 * mu * rhobar -6 * \
    gammabar * rho ) ) -2 * ( rhobar + rho ) * ( taubar*taubar ) ) ) ) ) \
    ) + Y22m[m+2] * ( ThornGqlm * ( 2 * mubar * taubar - mu * ( pi + taubar ) \
    ) + ( pi - taubar ) * ( ThornThornGqlm * ( mubar - mu ) * ( 1/( ( - \
    rhobar + rho ) ) ) -0.5 * ThornGqlm * ( 1/( rhobar ) ) * ( 1/( rho ) \
    ) * ( PsiCDe2 * rhobar + ( PsiCDe2bar * rho + 2 * ( rhobar + rho ) * \
    taubar * tau ) ) ) ) ) ) # m= +/- 2 or 1 has a canellation problem near the equatorial plane due to round-off error
        
        if m==2 or m==-2:
            middummy=np.real(dummythornphnm)*np.cos(theta)
            if theta<=np.pi/2+.01 and theta>=np.pi/2-.01:
                Thornphnm+= np.imag(dummythornphnm)*1j
            else: Thornphnm+= middummy/np.cos(theta) + np.imag(dummythornphnm)*1j
        elif m==1 or m==-1:
            middummy=np.imag(dummythornphnm)*np.cos(theta)
            if theta<=np.pi/2+.01 and theta>=np.pi/2-.01:
                Thornphnm+= np.real(dummythornphnm)
            else: Thornphnm+= middummy*1j/np.cos(theta) + np.real(dummythornphnm)
        else: Thornphnm+= dummythornphnm

    Thornphmm+=( complex( 0,-1 ) * Y22m[m+2] * a * m * rhobar * ( \
    ThornThornGqlm + 2 * ThornGqlm * rhobar ) * rho + Y22m[m+2] * ( 1/( \
    rhobar ) ) * ( 1/( rho ) ) * ( ThornThornThornGqlm * mu * rhobar + ( \
    2 * ThornGqlm * rhobar * ( PsiCDe2bar * rho + ( rhobar * ( PsiCDe2 + \
    ( mu * rhobar -4 * gammabar * rho ) ) -2 * ( rhobar + rho ) * ( \
    taubar*taubar ) ) ) + ThornThornGqlm * ( PsiCDe2bar * rho + ( rhobar \
    * ( PsiCDe2 + ( 2 * mu * rhobar -4 * gammabar * rho ) ) -2 * ( rhobar \
    + rho ) * ( taubar*taubar ) ) ) ) ) )
        
    # print('Ethphnndag=',Ethphnndag,'Thornphnm=',Thornphnm,'Ethhnndag=',Ethhnndag,'Thornhnm=',Thornhnm)

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
    # print(dChr)
    # transform to coordinate basis
    # dChrKN=np.zeros((4,4,4),np.complex128)
    # for ii in range(4): 
    #     for jj in range(4): 
    #         for kk in range(4): 
    #             lnis=[1,0,3,2] # apply eta matrix for (-+++), change ">1"->"<2" for (+---) or give overall minus, but doesn't change end result
    #             signi = 1 #if ii>1 else -1 
    #             signj = 1 if jj>1 else -1 
    #             signk = 1 if kk>1 else -1 
    #             dChrKN += signi*signj*signk*np.einsum('i,j,k',Tet[ii],InvTet[lnis[jj]],InvTet[lnis[kk]]) * dChr[ii, jj, kk]
    dChrBL = np.einsum('ai,jb,kc,ijk->abc',Tet,InvTet,InvTet, dChr)
    return dChrBL.real
    # gives some 1e-12 differences and 1e-12 imaginary parts for very large perturbations


def TotalGamma(a,radius,theta,phi_star,z_array=None,M=1):#np.array([0,0,0,0,0])):
    '''
    Claculate the Christoffel symbols for Kerr spacetime (-+++) and an arbitrary, static (no time-depenence) tidal (l=2) perturbation.
    The tidal perturbation is defined by (5) real parameters z_m for m in range(-l,l). Uses Outgoing Kerr coordinates.
    Default z_array is zero, which returns only background Christoffel symbols.
    Returns a 3-d numpy array in 3+1 dimensions
    '''
    # CosTH=np.cos(theta)
    # SinTH=np.sin(theta)

    # rho=(-1/(radius - 1.j*a*CosTH))
    # rhobar=np.conj(rho) # (-1/(radius + 1.j*a*CosTH))

    # PsiCDe2=M*rho*rho*rho
    # PsiCDe2bar=M*rhobar*rhobar*rhobar

    # Delta= a*a - 2*M*radius + radius*radius
    # Deltap= 2*(radius-M)
    # Deltapp= 2
    # Sigma= 1/rho/rhobar

    # # Non-zero Spin Coefficients (in terms of rho)
    # tau= -rho*rhobar*(1.j*a*SinTH)/sqrtof2
    # taubar= -tau
    # pi= taubar*rho/rhobar # is np.conj faster?
    # pibar= tau*rhobar/rho
    # beta= (-rhobar/(2*sqrtof2)*CosTH/SinTH)
    # # betabar= beta*rho/rhobar
    # # eps=0
    # # epsbar=0
    # mu= rho*rho*rhobar*Delta/2
    # mubar= mu*rhobar/rho
    # gamma= mu + rho*rhobar*Deltap/4 # Computation time may be shortened by a simplification that removes gamma and gammabar
    # gammabar= gamma-mu+mubar
    # # alpha= pi-betabar
    # # alphabar= pibar-beta

    Chr=Gamma(M,a,radius,theta)

    if z_array is not None:#np.linalg.norm(z_array)!=0:
        dChr=dGamma(M,a,radius,theta,phi_star,z_array)
        return Chr+dChr

## Average computation time for this function is > ~2ms
#print((time.perf_counter()-start)/10000)

    return Chr

    # totalG=Chr+dChr
if __name__ == '__main__':
    # restart=time.perf_counter()
    # imtolfactor=10**3 # imaginary parts less than imtolfactor*machineprecision will be suppressed.
    # #As far as I have currently found, the largest variances of dChr[a,b,c] from dChr[a,c,b], or any Im[] from zero, or q=0 from q!=0 are ~e-11.
    # theta_p,phi_p=np.pi/2,np.pi/2 # phi_KN ~ phi_BL at large distances (such as the perturber distances), but
    # # expresions in BL vs Cartesion vs Harmonic vs outgoing KN
    # sqrtpiover5=np.sqrt(np.pi/5)
    # CosTHp=np.cos(theta_p)
    # SinTHp=np.sin(theta_p)
    # expPHp=np.exp(1.j*phi_p)
    # z0,z1,z2=2*sqrtof6*sqrtpiover5* ( 1 - 3 * CosTHp*CosTHp ), 12*sqrtpiover5 *expPHp *CosTHp*SinTHp, -6*sqrtpiover5 /expPHp/expPHp *SinTHp*SinTHp#Cartesion (BL?)
    # z_companion= np.array([np.conj(z2),-np.conj(z1),z0,z1,z2]) # Norm = 4*sqrt(6*pi/5)~7.8 for all theta and phi.
    # eps=10**(-5)
    zs= .5/10000*np.sqrt(np.pi/5)*np.array([0,-6,2,6,0])#eps*z_companion
    # for it in range(0,10000):

    #     t=0.001*it
    #     radius=6+0.0008*it
    #     theta=0.0001+3.14*3/2*it/10000
    #     phi=4.4*it/10000    
    #     ChrTot=TotalGamma(.75,radius,theta,phi,zs)
    #     # Check for complex values
    #     it1, it2, it3 = 0,0,0
    #     for x in ChrTot: 
    #         for y in x: 
    #             for z in y:
    #                 if np.iscomplex(np.real_if_close(z,imtolfactor)) or abs(ChrTot[it1,it2,it3]-ChrTot[it1,it3,it2])>10**(-12):
    #                     print((it1,it2,it3),z,ChrTot[it1,it2,it3]-ChrTot[it1,it3,it2])
    #                 it3+=1
    #                 it3=it3%4
    #             it2+=1
    #             it2=it2%4
    #         it1+=1
    #         it1=it1%4

                

    # print("dChrKN calculated ",time.perf_counter()-start,time.perf_counter()-restart)


    # # # Average time for totalG is about the same
    # print((time.perf_counter()-start)/10000)



    # def dGamma2(M,a,radius,theta,phi_star,z_array):
    #     # if radius<=(M+np.sqrt(M*M-a*a)): return np.zeros((4,4,4)) # the test particle shouldn't cross this, this is for making 3D plots

    #     CosTH=np.cos(theta)
    #     SinTH=np.sin(theta)

    #     rho=(-1/(radius - 1.j*a*CosTH))
    #     rhobar=np.conj(rho) # (-1/(radius + 1.j*a*CosTH))

    #     PsiCDe2=M*rho*rho*rho
    #     PsiCDe2bar=M*rhobar*rhobar*rhobar

    #     Delta= a*a - 2*M*radius + radius*radius
    #     Deltap= 2*(radius-M)
    #     Deltapp= 2
    #     Sigma= 1/rho/rhobar

    #     # Non-zero Spin Coefficients (in terms of rho)
    #     tau= -rho*rhobar*(1.j*a*SinTH)/sqrtof2
    #     taubar= -tau
    #     pi= taubar*rho/rhobar # is np.conj faster?
    #     pibar= tau*rhobar/rho
    #     beta= (-rhobar/(2*sqrtof2)*CosTH/SinTH)
    #     # betabar= beta*rho/rhobar
    #     # eps=0
    #     # epsbar=0
    #     mu= rho*rho*rhobar*Delta/2
    #     mubar= mu*rhobar/rho
    #     gamma= mu + rho*rhobar*Deltap/4 # Computation time may be shortened by a simplification that removes gamma and gammabar
    #     gammabar= gamma-mu+mubar
    #     # alpha= pi-betabar
    #     # alphabar= pibar-beta

    #     # Kinnersley Tetrad in Outgoing Kerr-Newman Coordinates
    #     lKN=np.array([0,1,0,0])
    #     nKN=np.array([(radius*radius+a*a),-Delta/2,0,a])*rho*rhobar# decide if Delta is divided by 2 or not
    #     mKN=-rhobar/sqrtof2*np.array([1.j*a*SinTH,0,1,1.j/SinTH])
    #     mbKN=np.conj(mKN)#-rho/sqrtof2*np.array([-1.j*a*SinTH,0,1,-1.j/SinTH])

    #     #lower with g
    #     lKNd=np.array([-1, 0, 0, a* SinTH*SinTH]) 
    #     nKNd=np.array([(-Delta*rho*rhobar)/2, -1, 0, (1/2)*Delta*rho*rhobar * a*SinTH*SinTH])
    #     mKNd=np.array([1.j*rhobar* a*SinTH/sqrtof2, 0, -1/(sqrtof2*rho), -1.j*rhobar*(radius*radius + a*a)*SinTH/sqrtof2])
    #     mbKNd=np.conj(mKNd)#np.transpose(np.array([ 1.j*rho* a*SinTH/sqrtof2, 0, 1/(sqrtof2*rhobar), -1.j*rho*(radius*radius + a*a)*SinTH/sqrtof2]))

    #     Tet=np.array([lKN,nKN,mKN,mbKN])
    #     # raise with eta.
    #     InvTet=-np.array([-nKNd,-lKNd,mbKNd,mKNd]) # negate each (-1*) for (+---), but doesn't change end result

    #     sY2m=swsh_array(ell, theta, phi_star).reshape((2*ell+1,2*ell+1))
        
    #     # The z's may be complex, but have a reality condition
    #     # z[-m]=(-1)**m * np.conj(z[m]),
    #     # while Y[-s]l[-m] = (-1)**(m+s) * np.conj(Y[s]l[m]) so the conjugate relation is maintained.
    #     # So, np.conj(z[m]*Y[s]l[m])=(-1)**(s) * z[-m] * Y[-s]l[-m]
    #     # The z's are taken to be constants, and every instance of |s|=<2 is the result of a derivative.
    #     # s=0 can be reached by derivatives of both s=+/-1, but we only take derivatives s=2 and then conjugate the results.
    #     #Yn22m= np.conj(z_array)*sY2m[:,0] # s=-2 # not currently used in the calculation
    #     Yn12m= z_array*sY2m[:,1] # s=-1 # This only appears in Eth'hnnbar where three spin-lowering's occur
    #     Y02m = z_array*sY2m[:,2] # s= 0
    #     Y12m = z_array*sY2m[:,3] # s= 1
    #     Y22m = z_array*sY2m[:,4] # s= 2

    #     # Not Dependent on q
    #     # q=0
    #     # ell=2

    #     # for ell>2, these need to be either arrays or functions of m
    #     Gqlm= Delta*Delta/12
    #     # my h-components are in a form that only contain zG, ThornzG, and ThornThornzG.
    #     # dChr will contain first derivatives of the h-comp's which are also in terms of Thorn's only
    #     ThornGqlm= Deltap*Delta/6
    #     ThornThornGqlm= (Deltapp*Delta+Deltap*Deltap)/6
    #     ThornThornThornGqlm= Deltap#(*Deltapp/2)

    #     # Pieces of dGamma

    #     # hnn_2m= -Delta*Delta*(sqrtof6*rhobar*rhobar*Y02m+2*sqrtof2*rhobar*tau*Yn12m)
    #     # hnmb_2m= sqrtof2*rhobar*(rho+rhobar)*Delta*Delta*Yn12m + 2*Deltap*Delta*(sqrtof2*rhobar*Yn12m + (tau-pibar)*Yn22m)
    #     # hmbmb_2m= -(2*Deltap^2+Delta*(4+4*rho*Deltap))*Yn22m



    #     # At some point, we must sum over m
    #     # if ell==2
    #     hnndag=np.sum(Gqlm * rho * (2 * sqrtof2 * Y12m * taubar - sqrtof6 * Y02m * rho)
    #                     )
    #     hnm=np.sum(( sqrtof2 * Y12m * rho * ( ThornGqlm + Gqlm * ( rhobar + rho ) ) + ThornGqlm * Y22m * ( pi - taubar ) )
    #                 )
    #     hmm=np.sum(- Y22m * ( ThornThornGqlm + 2 * ThornGqlm * rhobar )
    #                 )
    #     # hnndag=np.sum(2 * Gqlm * ( 2 * Y12m * pi * rhobar - sqrtof6 * Y02m * ( rho*rho ) )
    #     #     )
    #     # hnm=np.sum(( 2 * Y12m * rho * ( ThornGqlm + Gqlm * ( rhobar + rho ) ) + ThornGqlm * Y22m * ( pi - taubar ) )
    #     #     )
    #     # hmm=np.sum(0 - Y22m * ( ThornThornGqlm + 2 * ThornGqlm * rhobar )
    #     #     )
    #     hnn=np.conj(hnndag)
    #     hnmb=np.conj(hnm)
    #     hmbmb=np.conj(hmm)
    #     # hnn=np.sum(- Gqlm * ( sqrtof6 * ( rhobar*rhobar ) * Y02m + 2*sqrtof2 * tau * rhobar * Yn12m)
    #     #     )
    #     # hnmb=np.sum(-0.5 * ( 2 * sqrtof2 * Gqlm * Yn12m * rhobar * ( rhobar + \
    #         # rho ) + 2 * ThornGqlm * ( sqrtof2 * Yn12m * rhobar + Yn22m * ( - pibar + tau ) ) )
    #     #     )
    #     # hmbmb=np.sum(0 - Y22m * ( ThornThornGqlm + 2 * ThornGqlm * rhobar )
    #     #     )
        
    #     # hnndag=np.conj(hnn)
    #     # hnm=np.conj(hnmb)
    #     # hmm=np.conj(hmbmb)

    #     # Each Thornp of h contains explicit m-dependence, so I initialize them to be summed by a loop.
    #     Thornphnndag=0
    #     Thornphnm=0
    #     Thornphmm=0

    #     Thornhnndag=np.sum(rho * ( - sqrtof6 * Y02m * rho * ( ThornGqlm + 2 * Gqlm \
    #     * rho ) + 2 * sqrtof2 * Y12m * ( ThornGqlm + Gqlm * ( rhobar + 2 * \
    #     rho ) ) * taubar )
    #     )

    #     Ethhnndag=np.sum(0.5 * Gqlm * ( rho * ( -6 * sqrtof2 * Y12m * rhobar * rho \
    #     + ( 8 * Y22m * rhobar * taubar + 4 * sqrtof6 * Y02m * ( rhobar + rho \
    #     ) * taubar ) ) + 2 * sqrtof2 * Y12m * ( 3 * PsiCDe2 * rhobar + ( \
    #     PsiCDe2bar * rho + ( 2 * mu * rhobar * ( - rhobar + rho ) -2 * ( 2 * \
    #     rhobar + rho ) * ( taubar*taubar ) ) ) ) )
    #     )
    #     Ethphnndag=np.sum(sqrtof2 * Gqlm * rho * ( 3 * Yn12m * ( rho*rho ) + 2 * \
    #     taubar * ( - sqrtof3 * Y02m * rho + Y12m * taubar ) )
    #     )
        
    #     Thornhnm=np.sum( sqrtof2 * Y12m * rho * ( ThornThornGqlm + ( ThornGqlm * \
    #     ( rhobar + 2 * rho ) + Gqlm * ( ( rhobar*rhobar ) + ( rhobar * rho + \
    #     2 * ( rho*rho ) ) ) ) ) + Y22m * ( ThornThornGqlm * ( pi - taubar ) \
    #     + ThornGqlm * ( 2 * pi * rho - ( rhobar + rho ) * taubar ) ) 
    #     )
        
    #     Ethhnm=np.sum( 2 * Y22m * rhobar * rho * ( ThornGqlm + Gqlm * ( rhobar + \
    #     rho ) ) + ( - sqrtof2 * Y12m * ( ThornGqlm * ( 2 * rhobar + rho ) + \
    #     Gqlm * ( rhobar + rho ) * ( rhobar + 2 * rho ) ) * taubar + 0.5 * \
    #     ThornGqlm * Y22m * ( 1/( rhobar ) ) * ( 1/( rho ) ) * ( 2 * mu * \
    #     rhobar * ( ( - rhobar + rho )*( - rhobar + rho ) ) + ( PsiCDe2bar * \
    #     rho * ( 3 * rhobar + rho ) + ( PsiCDe2 * rhobar * ( rhobar + 3 * rho \
    #     ) -2 * ( -2 * ( rhobar*rhobar ) + ( 3 * rhobar * rho + ( rho*rho ) ) \
    #     ) * ( taubar*taubar ) ) ) ) ) 
    #     )
    #     Ethphnm= -1/( rhobar ) * np.sum( sqrtof6 * Y02m * rhobar * ( rho*rho ) \
    #     * ( ThornGqlm + Gqlm * ( rhobar + rho ) ) + ( sqrtof2 * Y12m * rho \
    #     * ( - rhobar + 2 * rho ) * ( ThornGqlm + Gqlm * ( rhobar + rho ) ) \
    #     * taubar + ThornGqlm * Y22m * taubar * ( 2 * pi * rho + ( rhobar * \
    #     taubar - rho * taubar ) ) ) 
    #     )

    #     Thornhmm=- np.sum(Y22m * ( ThornThornThornGqlm + 2 * rhobar * ( \
    #     ThornThornGqlm + ThornGqlm * rhobar ) )
    #     )
        
    #     Ethhmm=np.sum(2 * Y22m * rhobar * ( ThornThornGqlm + ThornGqlm * rhobar ) \
    #     * ( 1/( rho ) ) * taubar
    #     )
    #     Ethphmm=np.sum( sqrtof2 * Y12m * ( ThornThornGqlm + 2 * ThornGqlm * \
    #     rhobar ) * rho + 2 * Y22m * ( 1/( rhobar ) ) * ( ThornThornGqlm * \
    #     rho + ThornGqlm * rhobar * ( - rhobar + 2 * rho ) ) * taubar 
    #     )

    #     for m in [-2,-1,0,1,2]:
    #         Thornphnndag+=( complex( 0,-1 ) * Gqlm * a * m * rhobar * ( rho*rho ) \
    #     * ( sqrtof6 * Y02m[m+2] * rho -2 * sqrtof2 * Y12m[m+2] * taubar ) + ( sqrtof2 * \
    #     Y12m[m+2] * ( 1/( rho ) ) * ( -2 * ThornGqlm * mu * rho * taubar + Gqlm \
    #     * ( PsiCDe2 * ( pibar * rhobar - rho * taubar ) + 2 * rho * taubar * \
    #     ( -4 * mu * rhobar + ( 6 * gammabar * rho + ( mu * rho + taubar * ( \
    #     pi + taubar ) ) ) ) ) ) + sqrtof6 * Y02m[m+2] * rho * ( ThornGqlm * mu + \
    #     Gqlm * ( PsiCDe2 + ( 2 * ( -2 * ( gammabar + gamma ) + mu ) * rho + \
    #     ( PsiCDe2bar * ( 1/( rhobar ) ) * rho + 2 * ( pi + taubar ) * tau ) ) \
    #     ) ) ) )
    #         Thornphnm+=( ThornGqlm * ( 4 * Y22m[m+2] * gammabar + complex( 0,1 ) * \
    #     Y22m[m+2] * a * m * rhobar * rho ) * ( pi - taubar ) + ( sqrtof2 * ( \
    #     complex( 0,1 ) * Y12m[m+2] * a * m * rhobar * ( rho*rho ) * ( ThornGqlm + \
    #     Gqlm * ( rhobar + rho ) ) - Y12m[m+2] * ( 1/( rhobar ) ) * ( 1/( ( - \
    #     rhobar + rho ) ) ) * ( Gqlm * ( rhobar * ( -3 * mu * ( \
    #     rhobar*rhobar*rhobar ) + ( 2 * ( 3 * gammabar + mu ) * ( \
    #     rhobar*rhobar ) * rho + ( ( rho*rho ) * ( PsiCDe2 -6 * gammabar * rho \
    #     ) + rhobar * rho * ( - PsiCDe2bar + mu * rho ) ) ) ) -2 * ( - rhobar \
    #     + rho ) * ( ( rhobar + rho )*( rhobar + rho ) ) * ( taubar*taubar ) ) \
    #     + ( - rhobar + rho ) * ( ThornThornGqlm * mu * rhobar + ThornGqlm * \
    #     ( PsiCDe2bar * rho + ( rhobar * ( PsiCDe2 + ( 3 * mu * rhobar -6 * \
    #     gammabar * rho ) ) -2 * ( rhobar + rho ) * ( taubar*taubar ) ) ) ) ) \
    #     ) + Y22m[m+2] * ( ThornGqlm * ( 2 * mubar * taubar - mu * ( pi + taubar ) \
    #     ) + ( pi - taubar ) * ( ThornThornGqlm * ( mubar - mu ) * ( 1/( ( - \
    #     rhobar + rho ) ) ) -0.5 * ThornGqlm * ( 1/( rhobar ) ) * ( 1/( rho ) \
    #     ) * ( PsiCDe2 * rhobar + ( PsiCDe2bar * rho + 2 * ( rhobar + rho ) * \
    #     taubar * tau ) ) ) ) ) )
    #         Thornphmm+=( complex( 0,-1 ) * Y22m[m+2] * a * m * rhobar * ( \
    #     ThornThornGqlm + 2 * ThornGqlm * rhobar ) * rho + Y22m[m+2] * ( 1/( \
    #     rhobar ) ) * ( 1/( rho ) ) * ( ThornThornThornGqlm * mu * rhobar + ( \
    #     2 * ThornGqlm * rhobar * ( PsiCDe2bar * rho + ( rhobar * ( PsiCDe2 + \
    #     ( mu * rhobar -4 * gammabar * rho ) ) -2 * ( rhobar + rho ) * ( \
    #     taubar*taubar ) ) ) + ThornThornGqlm * ( PsiCDe2bar * rho + ( rhobar \
    #     * ( PsiCDe2 + ( 2 * mu * rhobar -4 * gammabar * rho ) ) -2 * ( rhobar \
    #     + rho ) * ( taubar*taubar ) ) ) ) ) )
            

    #     Thornhnn=np.conj(Thornhnndag)
    #     # Thornhnmb=np.conj(Thornhnm)
    #     Thornhmbmb=np.conj(Thornhmm)

    #     Thornphnn=np.conj(Thornphnndag)
    #     Thornphnmb=np.conj(Thornphnm)
    #     Thornphmbmb=np.conj(Thornphmm)

    #     Ethhnn=np.conj(Ethphnndag)
    #     Ethhnmb=np.conj(Ethphnm)
    #     Ethhmbmb=np.conj(Ethphmm)

    #     Ethphnn=np.conj(Ethhnndag)
    #     # Ethphnmb=np.conj(Ethhnm)
    #     Ethphmbmb=np.conj(Ethhmm)

    #     # PDhnn_m= 
    #     # PDhnmb_m= 
    #     # PDhmbmb_m= 



    #     dChr=np.zeros((4,4,4),np.complex128) # complex numbers are being cast as real
    #     # Some terms have explicit q-dependence, but all boost-weight dependence vanishes in the final summation.
    #     dChr[0,0,1]=0.5 * ( Thornhnn + Thornhnndag) - ( hnmb * pibar + hnm * pi )
    #     dChr[0,0,2]=0.5 * ( Thornhnm +  hnm * rhobar - hmm * ( pi + taubar ) )
    #     dChr[0,0,3]=np.conj(dChr[0,0,2]) # 0.5 * ( Thornhnmb + ( hnmb * rho - hmbmb * ( pibar + tau ) ) )
    #     dChr[0,1,0]=dChr[0,0,1]
    #     dChr[0,1,1]=0.5 * ( Thornphnn + Thornphnndag )
    #     dChr[0,1,2]=0.5 * ( Ethhnn + Ethhnndag ) - hnm * mu
    #     dChr[0,1,3]=np.conj(dChr[0,1,2])
    #     dChr[0,2,0]=dChr[0,0,2]
    #     dChr[0,2,1]=dChr[0,1,2]
    #     dChr[0,2,2]=( Ethhnm -0.5 * Thornphmm - ( hmm * mu + hnm * tau ) )
    #     dChr[0,2,3]=0.5 * ( Ethhnmb + Ethphnm + ( hnn + hnndag ) * ( rhobar + rho ) - ( hnm * taubar + hnmb * tau ) )
    #     dChr[0,3,0]=dChr[0,0,3]
    #     dChr[0,3,1]=dChr[0,1,3]
    #     dChr[0,3,2]=dChr[0,2,3]
    #     dChr[0,3,3]=np.conj(dChr[0,2,2]) # ( Ethphnmb + ( -0.5 * Thornphmbmb + ( 0 - hmbmb * mubar - hnmb * taubar ) ) )

    #     dChr[1,1,1]=( 0.5 * ( 0 -  Thornhnn - Thornhnndag ) + ( hnm * ( pi + taubar ) + hnmb * ( pibar + tau ) ) )
    #     dChr[1,1,2]=0.5 * ( 0 - Thornhnm + ( hnm * rhobar + hmm * ( pi + taubar ) ) )
    #     dChr[1,1,3]=np.conj(dChr[1,1,2])

    #     dChr[1,2,1]=dChr[1,1,2]
    #     dChr[1,2,2]=( -0.5 * Thornhmm + hmm * rhobar )

    #     dChr[1,3,1]=dChr[1,1,3]

    #     dChr[1,3,3]=np.conj(dChr[1,2,2])

    #     dChr[2,0,1]=dChr[1,1,3] - hmbmb*tau # 0.5 * ( 0 - Thornhnmb + ( hnmb * rho + hmbmb * ( pibar - tau ) ) )

    #     dChr[2,0,3]=-0.5 * Thornhmbmb
    #     dChr[2,1,0]=dChr[2,0,1]
    #     dChr[2,1,1]=0.5 * ( Ethphnn + Ethphnndag ) - ( Thornphnmb + hnmb * mubar + ( hnn + hnndag ) * taubar )
    #     dChr[2,1,2]=0.5 * ( Ethphnm - Ethhnmb + ( hnn + hnndag ) * ( rho - rhobar ) - ( hnm * taubar + hnmb * tau ) )
    #     dChr[2,1,3]=( -0.5 * Thornphmbmb - hnmb * taubar )

    #     dChr[2,2,1]=dChr[2,1,2] # 0.5 * ( 0 - Ethhnmb + ( Ethphnm + ( 0 - hnn * rhobar + ( 0 - hnndag * rhobar + ( hnn * rho + ( hnndag * rho + ( 0 - hnm * taubar - hnmb * tau ) ) ) ) ) ) )
    #     dChr[2,2,2]=( 0.5 * Ethphmm + hnm * ( rho - rhobar ) )
    #     dChr[2,2,3]=( -0.5 * Ethhmbmb - hnmb * rhobar )
    #     dChr[2,3,0]=dChr[2,0,3]
    #     dChr[2,3,1]=dChr[2,1,3]
    #     dChr[2,3,2]=dChr[2,2,3]
    #     dChr[2,3,3]=-0.5 * Ethphmbmb

    #     dChr[3,0,1]=np.conj(dChr[2,0,1]) #dChr[1,1,2] -hmm*taubar # 0.5 * ( 0 - Thornhnm + ( hnm * rhobar + hmm * ( pi - taubar ) ) )
    #     dChr[3,0,2]=-0.5 * Thornhmm # np.conj(dChr[2,0,3])

    #     dChr[3,1,0]=dChr[3,0,1]
    #     dChr[3,1,1]=np.conj(dChr[2,1,1])
    #     dChr[3,1,2]=np.conj(dChr[2,1,3]) # ( -0.5 * Thornphmm - hnm * tau )
    #     dChr[3,1,3]=np.conj(dChr[2,1,2])
    #     dChr[3,2,0]=-0.5 * Thornhmm
    #     dChr[3,2,1]=dChr[3,1,2]
    #     dChr[3,2,2]=-0.5 * Ethhmm
    #     dChr[3,2,3]=np.conj(dChr[2,3,2]) # ( -0.5 * Ethphmm - hnm * rho )

    #     dChr[3,3,1]=dChr[3,1,3]
    #     dChr[3,3,2]=dChr[3,2,3]
    #     dChr[3,3,3]=np.conj(dChr[2,2,2]) # ( 0.5 * Ethhmbmb + hnmb * ( rhobar - rho ) )

    #     # def etaouter(){

    #     # }

    #     # transform to coordinate basis
    #     # dChrKN=np.zeros((4,4,4),np.complex128)
    #     # for ii in range(4): 
    #     #     for jj in range(4): 
    #     #         for kk in range(4): 
    #     #             lnis=[1,0,3,2] # apply eta matrix for (-+++), change ">1"->"<2" for (+---) or give overall minus
    #     #             signi = 1 #if ii>1 else -1 
    #     #             signj = 1 if jj<2 else -1 
    #     #             signk = 1 if kk<2 else -1 
    #     #             dChrKN += signi*signj*signk*np.einsum('i,j,k',Tet[ii],InvTet[lnis[jj]],InvTet[lnis[kk]]) * dChr[ii, jj, kk]
    #     dChrKN = np.einsum('ia,jb,kc,ijk->abc',Tet,InvTet,InvTet, dChr)
    #     return dChrKN

    M,a,radius,theta,phi_star,z_array=1.0,.1,8.91, np.pi/2, np.pi/4,zs#np.array([2.-2.j,-1.+1.j,1,1.+1.j,2.+2.j])#zs
    dGamma(M,a,radius,theta,phi_star,z_array)
    # print(dGamma(M,a,radius,theta,phi_star,z_array))#-dGamma2(M,a,radius,theta,phi_star,z_array))