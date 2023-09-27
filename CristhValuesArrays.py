import numpy as np
from scri import sph

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
# sqrtofPi=np.sqrt(np.pi)
# sqrt5overPi=1.26156626101008
# sqrt5over2Pi=0.892062058076386
# sqrt5over3Pi=0.7283656203947194
# sqrt5over6Pi=0.5150322693642528
# sqrt10over3Pi=1.030064538728506

CosTH=np.cos(theta)
SinTH=np.sin(theta)

rho=(-1/(radius - 1.j*a*CosTH))
rhobar=np.conj(rho) # (-1/(radius + 1.j*a*CosTH))

PsiCDe2=M*rho*rho*rho

Delta= a*a - 2*M*radius + radius*radius
Deltap= 2*(radius-M)
Deltapp= 2

Sigma= 1/rho/rhobar

# Kinnersley Tetrad in Kerr-Newman Coordinates
lKN=np.array([0,0,0,0])
nKN=np.array([(radius*radius+a*a),-Delta/1,0,a])*rho*rhobar
mKN=-rhobar/sqrtof2*np.array([1.j*a*SinTH,0,0,0.j/SinTH])
mbKN=-rho/sqrtof2*np.array([-1.j*a*SinTH,0,0,-1.j/SinTH])
lKNd=np.transpose(lKNlowerred)
nKNd=np.transpose(nKNlowerred)
mKNd=np.transpose(mKNlowerred)
mbKNd=np.transpose(mbKNlowerred)

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
Yn22m=
Yn12m=
Y02m=
Y12m=
Y22m=

# Common Derivatives
#   of Spin Coefs
# D\[Rho] -> \[Rho]^2
# D\[Rho]bar -> \[Rho]bar^2
# D\[Tau] -> (\[Rho] + \[Rho]bar) \[Tau]
# D\[Tau]bar -> (\[Rho] + \[Rho]bar) \[Tau]bar
# D\[CurlyPi] -> 2 \[CurlyPi] \[Rho]
# D\[CurlyPi]bar -> 2 \[CurlyPi]bar \[Rho]bar
# \[CapitalDelta]\[Rho] -> -\[Mu] \[Rho]
# \[CapitalDelta]\[Rho]bar -> -\[Mu]bar \[Rho]bar
# \[CapitalDelta]\[Tau] -> \[Mu]bar (\[CurlyPi] - \[Tau])
# \[CapitalDelta]\[Tau]bar -> -\[Mu]bar (\[CurlyPi] - \[Tau])
# \[CapitalDelta]\[CurlyPi] -> -2 \[Mu] \[CurlyPi]
# \[CapitalDelta]\[CurlyPi]bar -> -2 \[Mu]bar \[CurlyPi]bar
# \[Delta]\[Rho] -> -\[CurlyPi] \[Rho]bar
# \[Delta]\[Rho]bar -> -\[CurlyPi]bar \[Rho]bar
# \[Delta]\[Tau] -> \[Tau] (-2 \[Beta] + \[CurlyPi]bar + \[Tau])
# \[Delta]\[Tau]bar -> (2 \[Beta] - \[CurlyPi]bar + \[Tau]) \[Tau]bar
# \[Delta]\[CurlyPi] -> \[Tau] (-2 \[Beta]bar + \[CurlyPi] - \[Tau]bar)
# \[Delta]\[CurlyPi]bar -> -\[CurlyPi]bar^2 - (-2 \[Beta]bar + \[CurlyPi]) \[Tau]
# \[Delta]S\[Rho] -> -\[CurlyPi] \[Rho]
# \[Delta]S\[Rho]bar -> -\[CurlyPi]bar \[Rho]
# \[Delta]S\[Tau] -> \[Tau] (4 \[Beta]bar - 2 \[CurlyPi] + \[Tau]bar)
# \[Delta]S\[Tau]bar -> \[Tau]bar^2
# \[Delta]S\[CurlyPi] -> -\[CurlyPi]^2
# \[Delta]S\[CurlyPi]bar -> (-4 \[Beta] + 2 \[CurlyPi]bar - \[Tau]) \[Tau]bar

# #   of tetrad vectors
# Dell=
# Delm=
q=0

zG_qlm=
ThornzG_qlm=

ThornpzG_qlm=
EthzG_qlm=
EthpzG_qlm=np.conj(EthzG_qlm)
dzG_qlm=

# Pieces of dGamma
hnn_m= -Delta*Delta*(sqrtof6*rhobar*rhobar*Y02m+2*sqrtof2*rhobar*tau*Yn12m)
hnmb_m= sqrtof2*rhobar*(rho+rhobar)*Delta*Delta*Yn12m + 2*Deltap*Delta*(sqrtof2*rhobar*Yn12m + (tau-pibar)*Yn22m)
hmbmb_m= -(2*Deltap^2+Delta*(4+4*rho*Deltap))*Yn22m

# PDhnn_m= 
# PDhnmb_m= 
# PDhmbmb_m= 


dChr=np.zeros((4,4,4))

dChr[0,0,1]=( 0.5 * Thornhnn + ( -1 * hnm * pi + ( 0.5 * \
np.conjugate( Thornhnn ) + -1 * np.conjugate( hnm ) * np.conjugate( \
pi ) ) ) )
dChr[0,0,2]=( 0.5 * Thornhnm + ( -0.5 * hmm * pi + ( 0.5 * hnm * \
np.conjugate( rho ) + -0.5 * hmm * np.conjugate( tau ) ) ) )
dChr[0,0,3]=( -0.5 * tau * np.conjugate( hmm ) + ( 0.5 * rho * \
np.conjugate( hnm ) + ( 0.5 * np.conjugate( Thornhnm ) + -0.5 * \
np.conjugate( hmm ) * np.conjugate( pi ) ) ) )
dChr[0,1,0]=dChr[0,0,1]
dChr[0,1,1]=( 0.5 * Thornphnn + 0.5 * np.conjugate( Thornphnn ) )
dChr[0,1,2]=( 0.5 * Ethhnn + ( -1 * hnm * mu + 0.5 * np.conjugate( \
Ethphnn ) ) )
dChr[0,1,3]=( 0.5 * Ethphnn + ( 0.5 * np.conjugate( Ethhnn ) + -1 * \
np.conjugate( hnm ) * np.conjugate( mu ) ) )
dChr[0,2,0]=dChr[0,0,2]
dChr[0,2,1]=dChr[0,1,2]
dChr[0,2,2]=( Ethhnm + ( -0.5 * Thornphmm + ( -1 * hmm * mu + -1 * \
hnm * tau ) ) )
dChr[0,2,3]=( 0.5 * Ethphnm + ( 0.5 * hnn * rho + ( 0.5 * \
np.conjugate( Ethphnm ) + ( -0.5 * tau * np.conjugate( hnm ) + ( 0.5 \
* rho * np.conjugate( hnn ) + ( 0.5 * hnn * np.conjugate( rho ) + ( \
0.5 * np.conjugate( hnn ) * np.conjugate( rho ) + -0.5 * hnm * \
np.conjugate( tau ) ) ) ) ) ) ) )
dChr[0,3,0]=dChr[0,0,3]
dChr[0,3,1]=dChr[0,1,3]
dChr[0,3,2]=dChr[0,2,3]
dChr[0,3,3]=( np.conjugate( Ethhnm ) + ( -0.5 * np.conjugate( \
Thornphmm ) + ( -1 * np.conjugate( hmm ) * np.conjugate( mu ) + -1 * \
np.conjugate( hnm ) * np.conjugate( tau ) ) ) )

dChr[1,1,1]=( -0.5 * Thornhnn + ( hnm * pi + ( tau * np.conjugate( \
hnm ) + ( -0.5 * np.conjugate( Thornhnn ) + ( np.conjugate( hnm ) * \
np.conjugate( pi ) + hnm * np.conjugate( tau ) ) ) ) ) )
dChr[1,1,2]=( -0.5 * Thornhnm + ( 0.5 * hmm * pi + ( 0.5 * hnm * \
np.conjugate( rho ) + 0.5 * hmm * np.conjugate( tau ) ) ) )
dChr[1,1,3]=( 0.5 * tau * np.conjugate( hmm ) + ( 0.5 * rho * \
np.conjugate( hnm ) + ( -0.5 * np.conjugate( Thornhnm ) + 0.5 * \
np.conjugate( hmm ) * np.conjugate( pi ) ) ) )

dChr[1,2,1]=dChr[1,1,2]
dChr[1,2,2]=( -0.5 * Thornhmm + hmm * np.conjugate( rho ) )

dChr[1,3,1]=dChr[1,1,3]

dChr[1,3,3]=( rho * np.conjugate( hmm ) + -0.5 * np.conjugate( \
Thornhmm ) )

dChr[2,0,1]=( -0.5 * tau * np.conjugate( hmm ) + ( 0.5 * rho * \
np.conjugate( hnm ) + ( -0.5 * np.conjugate( Thornhnm ) + 0.5 * \
np.conjugate( hmm ) * np.conjugate( pi ) ) ) )

dChr[2,0,3]=-0.5 * np.conjugate( Thornhmm )
dChr[2,1,0]=dChr[2,0,1]
dChr[2,1,1]=( 0.5 * Ethphnn + ( 0.5 * np.conjugate( Ethhnn ) + ( -1 * \
np.conjugate( Thornphnm ) + ( -1 * np.conjugate( hnm ) * \
np.conjugate( mu ) + ( -1 * hnn * np.conjugate( tau ) + -1 * \
np.conjugate( hnn ) * np.conjugate( tau ) ) ) ) ) )
dChr[2,1,2]=( 0.5 * Ethphnm + ( 0.5 * hnn * rho + ( -0.5 * \
np.conjugate( Ethphnm ) + ( -0.5 * tau * np.conjugate( hnm ) + ( 0.5 \
* rho * np.conjugate( hnn ) + ( -0.5 * hnn * np.conjugate( rho ) + ( \
-0.5 * np.conjugate( hnn ) * np.conjugate( rho ) + -0.5 * hnm * \
np.conjugate( tau ) ) ) ) ) ) ) )
dChr[2,1,3]=( -0.5 * np.conjugate( Thornphmm ) + -1 * np.conjugate( \
hnm ) * np.conjugate( tau ) )

dChr[2,2,1]=dChr[2,1,2]
dChr[2,2,2]=( 0.5 * Ethphmm + ( hnm * rho + -1 * hnm * np.conjugate( \
rho ) ) )
dChr[2,2,3]=( -0.5 * np.conjugate( Ethphmm ) + -1 * np.conjugate( hnm \
) * np.conjugate( rho ) )
dChr[2,3,0]=dChr[2,0,3]
dChr[2,3,1]=dChr[2,1,3]
dChr[2,3,2]=dChr[2,2,3]
dChr[2,3,3]=-0.5 * np.conjugate( Ethhmm )

dChr[3,0,1]=( -0.5 * Thornhnm + ( 0.5 * hmm * pi + ( 0.5 * hnm * \
np.conjugate( rho ) + -0.5 * hmm * np.conjugate( tau ) ) ) )
dChr[3,0,2]=-0.5 * Thornhmm

dChr[3,1,0]=dChr[3,0,1]
dChr[3,1,1]=( 0.5 * Ethhnn + ( -1 * Thornphnm + ( -1 * hnm * mu + ( \
-1 * hnn * tau + ( 0.5 * np.conjugate( Ethphnn ) + -1 * tau * \
np.conjugate( hnn ) ) ) ) ) )
dChr[3,1,2]=( -0.5 * Thornphmm + -1 * hnm * tau )
dChr[3,1,3]=( -0.5 * Ethphnm + ( -0.5 * hnn * rho + ( 0.5 * \
np.conjugate( Ethphnm ) + ( -0.5 * tau * np.conjugate( hnm ) + ( -0.5 \
* rho * np.conjugate( hnn ) + ( 0.5 * hnn * np.conjugate( rho ) + ( \
0.5 * np.conjugate( hnn ) * np.conjugate( rho ) + -0.5 * hnm * \
np.conjugate( tau ) ) ) ) ) ) ) )
dChr[3,2,0]=dChr[3,0,2]
dChr[3,2,1]=dChr[3,1,2]
dChr[3,2,2]=-0.5 * Ethhmm
dChr[3,2,3]=( -0.5 * Ethphmm + -1 * hnm * rho )

dChr[3,3,1]=dChr[3,1,3]
dChr[3,3,2]=dChr[3,2,3]
dChr[3,3,3]=( 0.5 * np.conjugate( Ethphmm ) + ( -1 * rho * \
np.conjugate( hnm ) + np.conjugate( hnm ) * np.conjugate( rho ) ) )


# Delll=
# Dellmml=
# Delmm=

# dGammahnn_m= 
# dGammahnmb_m= 
# dGammahmbmb_m= 