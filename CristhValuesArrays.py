import numpy as np
import spherical
import quaternionic

t=0.1
radius=4
theta=1.3
phi=0

M=1
a=.9

zs=[0,0,0,0,0]

sqrtof2=np.sqrt(2)
sqrtof3=np.sqrt(3)
sqrtof5=np.sqrt(5)
sqrtof6=sqrtof2*sqrtof3

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

# Kinnersley Tetrad in Kerr-Newman Coordinates
# lKN=np.array([0,1,0,0])
# nKN=np.array([(radius*radius+a*a),-Delta/1,0,a])*rho*rhobar
# mKN=-rhobar/sqrtof2*np.array([1.j*a*SinTH,0,0,0.j/SinTH])
# mbKN=-rho/sqrtof2*np.array([-1.j*a*SinTH,0,0,-1.j/SinTH])
# lKNd=np.transpose(lKNlowerred)
# nKNd=np.transpose(nKNlowerred)
# mKNd=np.transpose(mKNlowerred)
# mbKNd=np.transpose(mbKNlowerred)

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
sY2m=swsh_array(2, theta, phi).reshape((2*ell+1,2*ell+1))
Yn22m=sY2m[:,0]
Yn12m=sY2m[:,1]
Y02m=sY2m[:,2]
Y12m=sY2m[:,3]
Y22m=sY2m[:,4]


# #   of tetrad vectors
# Dell=
# Delm=
q=0
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
# Each Thornp of h contains explicit m-dependence, so I initialize them to be summed by a loop.
Thornphnndag=0
Thornphnm=0
Thornphmm=0
# if ell==2
if q==0:
    hnndag=np.sum(2 * zGqlm * ( 2 * Y12m * pi * rhobar - sqrtof6 * Y02m * ( rho*rho ) )
    )
    hnm=np.sum(( 2 * Y12m * rho * ( ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ThornzGqlm * Y22m * ( pi - taubar ) )
    )
    hmm=np.sum(0 - Y22m * ( ThornThornzGqlm + 2 * ThornzGqlm * rhobar )
    )

    Thornhnndag=np.sum(2 * rho * ( 0 - sqrtof6 * Y02m * rho * ( ThornzGqlm + 2 \
        * zGqlm * rho ) + 2 * Y12m * ( ThornzGqlm + zGqlm * ( rhobar + 2 * \
        rho ) ) * taubar )
    )
    #Thornphnndag=
    for m in [-2,-1,0,1,2]:
        Thornphnndag+= -2 * ( 1/( rhobar ) ) * ( complex( 0,1 ) * zGqlm * a * m \
            * ( rhobar*rhobar ) * ( rho*rho ) * ( sqrtof6 * Y02m[m+2] * rho -2 * \
            Y12m[m+2] * taubar ) + ( 0 - sqrtof6 * Y02m[m+2] * rho * ( ThornzGqlm * mu * \
            rhobar + zGqlm * ( PsiCDe2 * rhobar + ( PsiCDe2bar * rho + ( 2 * mu * \
            rhobar * rho + ( -2 * rhobar * ( taubar*taubar ) -2 * rho * ( \
            taubar*taubar ) ) ) ) ) ) + Y12m[m+2] * taubar * ( 2 * ThornzGqlm * mu * \
            rhobar - zGqlm * ( 0 - PsiCDe2 * rhobar + ( 0 - PsiCDe2bar * \
            rho + ( -2 * mu * rhobar * ( 2 * rhobar + rho ) + 2 * ( rhobar + rho \
            ) * ( taubar*taubar ) ) ) ) ) ) )
    Ethhnndag=np.sum(2 * zGqlm * ( rho * ( -6 * Y12m * rhobar * rho + ( 4 * Y22m \
        * rhobar * taubar + 2 * sqrtof6 * Y02m * ( rhobar + rho ) * taubar ) \
        ) + Y12m * ( 3 * PsiCDe2 * rhobar + ( PsiCDe2bar * rho + ( 2 * mu * \
        rhobar * ( 0 - rhobar + rho ) -2 * ( 2 * rhobar + rho ) * ( \
        taubar*taubar ) ) ) ) )
    )
    Ethphnndag=np.sum(4 * zGqlm * rho * ( 3 * Yn12m * ( rho*rho ) + taubar * ( \
        0 - sqrtof6 * Y02m * rho + Y12m * taubar ) )
    )
    Thornhnm=np.sum(( 2 * Y12m * rho * ( ThornThornzGqlm + ( ThornzGqlm * ( \
        rhobar + 2 * rho ) + zGqlm * ( ( rhobar*rhobar ) + rho * ( rhobar + 2 \
        * rho ) ) ) ) + Y22m * ( ThornThornzGqlm * ( pi - taubar ) + \
        ThornzGqlm * ( 2 * pi * rho - ( rhobar + rho ) * taubar ) ) )
    )
    # Thornphnm=
    for m in [-2,-1,0,1,2]:
        Thornphnm+=( complex( 0,2 ) * Y12m[m+2] * a * m * rhobar * ( rho*rho ) * ( \
            ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ( complex( 0,1 ) * \
            ThornzGqlm * Y22m[m+2] * a * m * rho * ( 0 - rhobar + rho ) * taubar + ( \
            ThornzGqlm * Y22m[m+2] * ( 2 * mubar * taubar - mu * ( pi + taubar ) \
            ) + ( Y12m[m+2] * ( 1/( rhobar ) ) * ( 1/( ( 0 - rhobar + rho ) ) ) * ( 2 \
            * zGqlm * ( PsiCDe2bar * ( rhobar*rhobar ) * rho + ( 0 - PsiCDe2 * \
            rhobar * ( rho*rho ) + ( mu * rhobar * ( ( rhobar*rhobar*rhobar ) + ( \
            rhobar * ( rho*rho ) -2 * ( rho*rho*rho ) ) ) + 2 * ( 0 - rhobar + \
            rho ) * ( ( rhobar + rho )*( rhobar + rho ) ) * ( taubar*taubar ) ) ) \
            ) + ( 0 - rhobar + rho ) * ( -2 * ThornThornzGqlm * mu * rhobar + \
            ThornzGqlm * ( -2 * PsiCDe2 * rhobar + ( -2 * PsiCDe2bar * rho + ( -2 \
            * mu * rhobar * ( rhobar + 2 * rho ) + 4 * ( rhobar + rho ) * ( \
            taubar*taubar ) ) ) ) ) ) + Y22m[m+2] * ( pi - taubar ) * ( \
            ThornThornzGqlm * ( mubar - mu ) * ( 1/( ( 0 - rhobar + rho ) ) \
            ) -0.5 * ThornzGqlm * ( 1/( rhobar ) ) * ( 1/( rho ) ) * ( PsiCDe2 \
            * rhobar + ( PsiCDe2bar * rho + 2 * ( rhobar + rho ) * taubar * tau ) \
            ) ) ) ) ) )
    Ethhnm=np.sum(( 4 * Y22m * rhobar * rho * ( ThornzGqlm + zGqlm * ( rhobar + \
        rho ) ) + ( -2 * Y12m * ( ThornzGqlm * ( 2 * rhobar + rho ) + zGqlm * \
        ( rhobar + rho ) * ( rhobar + 2 * rho ) ) * taubar + 0.5 * ThornzGqlm \
        * Y22m * ( 1/( rhobar ) ) * ( 1/( rho ) ) * ( 2 * mu * rhobar * ( ( \
        0 - rhobar + rho )*( 0 - rhobar + rho ) ) + ( PsiCDe2bar * rho * ( \
        3 * rhobar + rho ) + ( PsiCDe2 * rhobar * ( rhobar + 3 * rho ) -2 * \
        ( -2 * ( rhobar*rhobar ) + ( 3 * rhobar * rho + ( rho*rho ) ) ) * ( \
        taubar*taubar ) ) ) ) ) )
    )
    Ethphnm=np.sum(( 1/( rhobar ) ) * ( -2 * sqrtof6 * Y02m * rhobar * ( rho*rho \
        ) * ( ThornzGqlm + zGqlm * ( rhobar + rho ) ) + ( -2 * Y12m * rho * ( \
        0 - rhobar + 2 * rho ) * ( ThornzGqlm + zGqlm * ( rhobar + rho ) ) * \
        taubar + ThornzGqlm * Y22m * taubar * ( -2 * pi * rho + ( 0 - rhobar \
        + rho ) * taubar ) ) )
    )
    Thornhmm=np.sum(0 - Y22m * ( ThornThornThornzGqlm + 2 * rhobar * ( \
        ThornThornzGqlm + ThornzGqlm * rhobar ) )
    )
    # Thornphmm=
    for m in [-2,-1,0,1,2]:
        Thornphmm+= ( complex( 0,-1 ) * Y22m[m+2] * a * m * rhobar * ( \
            ThornThornzGqlm + 2 * ThornzGqlm * rhobar ) * rho + Y22m[m+2] * ( 1/( \
            rhobar ) ) * ( 1/( rho ) ) * ( ThornThornThornzGqlm * mu * rhobar + ( \
            2 * ThornzGqlm * rhobar * ( PsiCDe2 * rhobar + ( mu * ( rhobar*rhobar \
            ) + ( PsiCDe2bar * rho + ( -2 * rhobar * ( taubar*taubar ) -2 * rho \
            * ( taubar*taubar ) ) ) ) ) + ThornThornzGqlm * ( PsiCDe2 * rhobar + \
            ( 2 * mu * ( rhobar*rhobar ) + ( PsiCDe2bar * rho + ( -2 * rhobar * ( \
            taubar*taubar ) -2 * rho * ( taubar*taubar ) ) ) ) ) ) ) )
    Ethhmm=np.sum(2 * Y22m * rhobar * ( ThornThornzGqlm + ThornzGqlm * rhobar ) \
        * ( 1/( rho ) ) * taubar
    )
    Ethphmm=np.sum(( 2 * Y12m * ( ThornThornzGqlm + 2 * ThornzGqlm * rhobar ) * \
        rho + 2 * Y22m * ( 1/( rhobar ) ) * ( ThornThornzGqlm * rho + \
        ThornzGqlm * rhobar * ( 0 - rhobar + 2 * rho ) ) * taubar )
    )

# else:
#     hnndag=
#     hnm=
#     hmm=

#     Thornhnndag=
#     Thornhnm=
#     Thornhmm=

#     Thornphnndag=
#     Thornphnm=
#     Thornphmm=

#     Ethhnndag=
#     Ethhnm=
#     Ethhmm=

#     Ethphnndag=
#     Ethphnm=
#     Ethphmm=

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
dChr[0,2,0]=dChr[0,2,0]
dChr[0,2,1]=dChr[0,1,2]
dChr[0,2,2]=( Ethhnm -0.5 * Thornphmm - ( hmm * mu + hnm * tau ) )
dChr[0,2,3]=0.5 * ( Ethhnmb + Ethphnm + ( hnn + hnndag ) * ( rhobar + rho ) - ( hnm * taubar + hnmb * tau ) )
dChr[0,3,0]=dChr[0,0,3]
dChr[0,3,1]=dChr[0,1,3]
dChr[0,3,2]=dChr[0,2,3]
dChr[0,3,3]=np.conj(dChr[0,2,2]) # ( Ethphnmb + ( -0.5 * Thornphmbmb + ( 0 - hmbmb * mubar - hnmb * taubar ) ) )

dChr[1,1,1]=( 0.5 * ( 0 -  Thornhnn + Thornhnndag ) + ( hnm * ( pi + taubar ) + hnmb * ( pibar + tau ) ) )
dChr[1,1,2]=0.5 * ( 0 - Thornhnm + ( hnm * rhobar + hmm * ( pi + taubar ) ) )
dChr[1,1,3]=np.conj(dChr[1,1,2])

dChr[1,2,1]=dChr[1,1,2]
dChr[1,2,2]=( -0.5 * Thornhmm + hmm * rhobar )

dChr[1,3,1]=dChr[1,1,3]

dChr[1,3,3]=np.conj(dChr[1,2,2])

dChr[2,0,1]=dChr[1,1,3] # 0.5 * ( 0 - Thornhnmb + ( hnmb * rho + hmbmb * ( pibar - tau ) ) )

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

dChr[3,0,1]=dChr[1,1,2] # np.conj(dChr[2,0,1]) # 0.5 * ( 0 - Thornhnm + ( hnm * rhobar + hmm * ( pi - taubar ) ) )
dChr[3,0,2]=-0.5 * Thornhmm # np.conj(dChr[2,0,3])

dChr[3,1,0]=np.conj(dChr[3,0,1])
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

print(dChr[0,0:3,:])
# Delll=
# Dellmml=
# Delmm=

# dGammahnn_m= 
# dGammahnmb_m= 
# dGammahmbmb_m= 