#KerrB.py
import math
# for Schwarzhild conditions
def schwarzEnergy(a,p,e,x):
    if(e==0):return (-2 + p)/math.sqrt((-3 + p)* p)# circular orbit
    return  math.sqrt((-4 *e*e + (-2 + p)*(-2 + p))/(p* (-3 - e*e + p)))

def schwarzMomentum(a,p,e,x):return (p*x)/math.sqrt(-3 - e*e + p)
def schwarzQ(a,p,e,x):
    if(e==0):return 0-p*p*(-1 + x*x)/(-3 + p)# Circular orbit. The zero is necessary for the compiler, for some reason.
    return p*p *(-1 + x*x)/(3 + e*e - p)


# Equatorial orbits, not circular
def equatEnergy(a,p,e,x):
    return math.sqrt(1 - \
                     ((1 - e*e)* \
                                            (1 + \
                       ((-1 + e*e)* \
                                                (a*a *(1 + 3 *e*e + p) + p* \
                                                                          (-3 - e*e + p - 2* x* \
                                                                                                     math.sqrt((a**6 *(-1 + e*e)*(-1 + e*e) + a*a *(-4 *e*e + (-2 + p)*(-2 + p)) *p*p + 2 *a*a*a*a *p *(-2 + p + e*e *(2 + p))
                              )/(p*p*p *x*x)
                          )
                         )
                        )
                       )/(-4 *a*a *(-1 + e*e)*(-1 + e*e) + (3 + e*e - p)*(3 + e*e - p) *p)
                      )
                     )/p
                    )

def equatMomentum(a,p,e,x):
    return p*x* \
                        math.sqrt(
            (a*a *(1 + 3 *e*e + p) + p * \
                          (-3 - e*e + p - 2 *x * \
                                         math.sqrt(\
                  (a**6 *(-1 + e*e)*(-1 + e*e) + a*a *(-4 *e*e + (-2 + p)*(-2 + p)) *p*p + 2 *a*a*a*a *p *(-2 + p + e*e *(2 + p))\
                  )/(x*x *p*p*p)\
              )\
             )\
            )/(\
                (-4 *a*a *(-1 + e*e)*(-1 + e*e) + (3 + e*e - p)*(3 + e*e - p) *p)\
                *x*x\
            )\
        )\
        + a* \
                            math.sqrt(1 - \
                  ((1 - e*e) * \
                                      (1 + \
                    ((-1 + e*e) * \
                                          (a*a *(1 + 3 *e*e + p) + p * \
                                                                 (-3 - e*e + p - 2 *x * \
                                                                                         math.sqrt( \
                           (a**6 *(-1 + e*e)*(-1 + e*e) + a*a *(-4 *e*e + (-2 + p)*(-2 + p)) *p*p + 2 *a*a*a*a *p *(-2 + p + e*e *(2 + p)) \
                           )/(p*p*p *x*x) \
                       )\
                      )\
                     )\
                    )/(-4 *a*a *(-1 + e*e)*(-1 + e*e) + (3 + e*e - p)*(3 + e*e - p) *p) \
                   )\
                  )/p \
                 )\

def equatQ(a,p,e,x):return 0

# Equatorial and Circular
def circEquatEnergy(a,p,e,x):return ((-2 + p)* math.sqrt(p) + a/x) \
                                  /math.sqrt(2* a/x *p**(3/2) + (-3 + p)* p*p)

def circEquatMomentum(a,p,e,x):return (a*a - 2 *a/x *math.sqrt(p) + p*p) \
                                    /math.sqrt(2 *a/x + (-3 + p)* math.sqrt(p)) \
                                    /p**(3/4)

def circEquatQ(a,p,e,x):return 0

# Polar orbits, not spherical
def polEnergy(a,p,e,x):
    return math.sqrt( \
        -( \
            (p * \
                          (a*a*a*a *(-1 + e*e)*(-1 + e*e) + (-4 *e*e + (-2 + p)*(-2 + p)) *p*p + 2 *a*a *p *(-2 + p + e*e *(2 + p)) \
             ) \
            ) \
                /(a*a*a*a *(-1 + e*e)*(-1 + e*e) *(-1 + e*e - p) + (3 + e*e - p) *p*p*p*p - 2 *a*a *p*p *(-1 - e*e*e*e + p + e*e *(2 + p)) \
                 ) \
        ) \
    ) 

def polMomentum(a,p,e,x):return 0
def polQ(a,p,e,x):
    return -(p*p) * \
                            (a*a*a*a *(-1 + e*e)*(-1 + e*e) + p*p*p*p + 2 *a*a *p *(-2 + p + e*e *(2 + p))) \
          /(
              a*a*a*a *(-1 + e*e)*(-1 + e*e) *(-1 + e*e - p) + (3 + e*e - p) *p*p*p*p - 2 *a*a *p*p *(-1 - e*e*e*e + p + e*e *(2 + p))
           )


# Polar and Spherical
def polSpherEnergy(a,p,e,x):
    return  math.sqrt(
        (p* (a*a - 2* p + p*p)*(a*a - 2* p + p*p))/(
            (a*a + p*p) *(a*a + a*a *p - 3 *p*p + p*p*p)
        )
    )

def polSpherMomentum(a,p,e,x):return 0
def polSpherQ(a,p,e,x):
    return p*p *(a*a*a*a + 2 *a*a *(-2 + p) *p + p*p*p*p) \
        /(a*a + p*p) \
        /((-3 + p) *p*p + a*a *(1 + p))


# Other Spherical orbits
def spherEnergy(a,p,e,x):
    return math.sqrt(
        (
            (-3 + p)* (-2 + p)*(-2 + p) *p**5 - 2 *a**5 *x *(-1 + x*x)* \
            math.sqrt(p*p*p + a*a *p *(-1 + x*x)) + a*a*a*a *p*p *(-1 + x*x) * \
            (4 - 5* p* (-1 + x*x) + 3* p*p *(-1 + x*x)) - a**6 *(-1 + x*x)*(-1 + x*x) * \
            (x*x + p*p *(-1 + x*x) - p *(1 + 2 *x*x))+ a*a *p*p*p * \
            (4 - 4 *x*x + p *(12 - 7 *x*x) - 3 *p*p*p *(-1 + x*x) + p*p *(-13 + 10 *x*x))+ a* \
            (-2* p**(9/2) *x * \
                 math.sqrt(p*p + a*a *(-1 + x*x))
                 + 4* p*p*p *x * \
                math.sqrt(p*p*p + a*a *p *(-1 + x*x))
            )
            + 2 *a*a*a * \
            (2 *p *x *(-1 + x*x) * \
                math.sqrt(p*p*p + a*a *p *(-1 + x*x))
                 - x*x*x * \
                math.sqrt(p**7 + a*a *p**5 * (-1 + x*x))
            )
        )/(
            (p*p - a*a *(-1 + x*x)) * \
            ((-3 + p)*(-3 + p) *p*p*p*p - 2 *a*a *p*p * \
                (3 + 2 *p - 3 *x*x + p*p *(-1 + x*x))
                 + a*a*a*a *(-1 + x*x) * \
                (-1 + x*x + p*p *(-1 + x*x) - 2 *p *(1 + x*x))
            )
          )
    )

def spherMomentum(a,p,e,x,En):
    g = 2 *a *p
    d = (a*a + (-2 + p) *p) *(p*p - a*a *(-1 + x*x))
    h = ((-2 + p) *p - a*a *(-1 + x*x))/(x*x)
    f = p*p*p*p + a*a *(p *(2 + p) - (a*a + (-2 + p) *p) *(-1 + x*x))
    return (-En *g + x *math.sqrt((-d *h + En*En *(g*g + f *h))/(x*x)))/h

# def spherQ(a,p,e,x):

# (*Generic orbits*)
def genE(a,p,e,x):# I somehow broke this when trying to speed it up.
    r1 = p/(1 - e)
    r2 = p/(1 + e)
    zm = math.sqrt(1 - x*x)
    Delta1= r1*r1 - 2* r1 + a*a
    Delta2= r2*r2 - 2* r2 + a*a
    f1 = r1*r1*r1*r1 + a*a *(r1 *(r1 + 2) + zm*zm * Delta1)
    f2 = r2*r2*r2*r2 + a*a *(r2 *(r2 + 2) + zm*zm * Delta2)
    g1 = 2* a* r1
    g2 = 2* a* r2
    h1 = r1*(r1 - 2) + zm*zm/(1 - zm*zm)*Delta1
    h2 = r2*(r2 - 2) + zm*zm/(1 - zm*zm)*Delta2
    d1 = (r1*r1 + a*a * zm*zm)* Delta1
    d2 = (r2*r2 + a*a * zm*zm)* Delta2
    Kappa = d1* h2 - h1* d2
    Epsilon = d1* g2 - g1* d2
    Rho = f1* h2 - h1* f2
    Eta = f1* g2 - g1* f2
    Sigma = g1* h2 - h1* g2
    return math.sqrt((Kappa * Rho + 2* Epsilon * Sigma - x * 2* math.sqrt(Sigma * (Sigma * Epsilon*Epsilon + Rho * Epsilon * Kappa - Eta * Kappa*Kappa)/(x*x)))/(Rho*Rho + 4* Eta * Sigma))

def genL(a,p,e,x,En):
    r = p/(1 - e)
    zm = math.sqrt(1 - x*x)
    Delta =r*r - 2* r + a*a
    f = r*r*r*r + a*a *(r* (r + 2) + zm*zm * Delta)
    g = 2* a* r
    h = r*(r - 2) + zm*zm/(1 - zm*zm)*Delta
    d = (r*r + a*a * zm*zm)* Delta
    return (-En*g + x*math.sqrt((-d*h + En*En *(g*g + f*h))/(x*x)))/h

def genQ(a,p,e,x,En,L):
    zm = math.sqrt(1 - x*x)
    return zm*zm *(a*a * (1 - En*En) + L*L/(1 - zm*zm))


# Switchboard defs for the Constants of Motion
def KerrGeoEnergy(a, p, e, x) :
    if(a==0):return schwarzEnergy(a,p,e,x)
    if(e==0):
        if((x*x) == 1):return circEquatEnergy(a,p,e,x)
        if(x == 0):return polSpherEnergy(a,p,e,x)
        return spherEnergy(a,p,e,x)
               
    if(x==0):return polEnergy(a,p,e,x)
    if((x*x) == 1):return equatEnergy(a,p,e,x)

    return genE(a,p,e,x)



def KerrGeoAngularMomentum(a, p, e, x, En= None):
    if(x==0):return 0# p may be large
    if(a==0):return schwarzMomentum(a,p,e,x)
    if((x*x) == 1):
        if(e==0):return circEquatMomentum(a,p,e,x)
        return equatMomentum(a,p,e,x)

    if(En == None): En = KerrGeoEnergy(a, p, e, x)# just in case
    if(e==0):return spherMomentum(a,p,e,x,En)

    return genL(a,p,e,x,En) # not currently calculated by Controller



def KerrGeoCarterConstant(a, p, e, x, En= None, L= None):
    if((x*x) == 1) :return 0# I'm using 1 rather than ~1 for case of large p
    if(a==0):return schwarzQ(a,p,e,x)
    if(x==0):
        if(e==0):return polSpherQ(a,p,e,x)
        return polQ(a,p,e,x)

    if(En == None): En = KerrGeoEnergy(a, p, e, x)
    if(L == None):L = KerrGeoAngularMomentum(a, p, e, x, En)

    return genQ(a,p,e,x,En,L) # not currently calculated by Controller



def KerrGeoConstantsOfMotion(a, p, e, x):
    En = KerrGeoEnergy(a, p, e, x)
    L = KerrGeoAngularMomentum(a, p, e, x, En)
    Q = KerrGeoCarterConstant(a, p, e, x, En, L)
    return [En, L, Q]


# (*Roots of the radial and polar equations*)(*Returns the roots of the
# radial equation,as given by Fujita and Hikida*)
def KerrGeoRadialRoots(a, p, e, x, En= None, Q= None):
    if(En == None): En = KerrGeoEnergy(a, p, e, x)
    if(Q == None): Q = KerrGeoCarterConstant(a, p, e, x)
    M=1
    r1 = p/(1 - e)
    r2 = p/(1 + e)
    AplusB = (2*M)/(1 - En*En) - (r1 + r2)# (*Eq.(11)*)
    AB = (a*a *Q)/((1 - En*En) *r1 *r2)# (*Eq.(11)*)
    r3 = (AplusB + math.sqrt(AplusB*AplusB - 4*AB))/2# (*Eq.(11)*)
    r4 = AB/r3
    return [r1, r2, r3, r4]


def KerrGeoPolarRoots(a, p, e, x):
    En, L, Q = KerrGeoConstantsOfMotion(a, p, e, x)
    zm = math.sqrt(1 - x*x)
    zp = math.sqrt((a*a *(1 - En*En) + L*L/(1 - zm*zm)))
    return [zp, zm]

# These Root defs don't get called, but seem important somehow.

## Work on this if we need the Separatrix ##

# # (*Photon Sphere, Needed for IBSO*)
# def KerrGeoPhotonSphereRadius(a, x0):#  Abs[x0] <= 1] /  Precision[:a1, x0] != \[Infinity] :=
#   # Approximately Schwarzschild
#     if (a<=0.001) :
#         #  console.log('a==0')
#         return 3
    
#     # (*Radius of photon sphere for equatorial orbits from Bardeen,Press,Teukolsky ApJ,178,p347 (1972),Eq.2.18*)
#     if ((x0*x0)>=.98) :# Abs(x0)>~.99
#         #  console.log('(x0*x0)>.98')
#         return 2*(1 + math.cos(2/3 *math.acos(-x0*a)))
    

#     # (*For polar orbits the radius was given by E.Teo,General Relativity and Gravitation,v.35,Issue 11,p.1909-1926 (2003),Eq.(14)*)
#     if ((x0*x0) < 0.001) :# Abs(x0)<~0.0316
#         #  console.log('(x0*x0)<.01')
#         return 1+2*math.sqrt(1 - a*a /3)*math.cos(1/3 *math.acos((1 - a*a)/(1 - 1/3 *a*a)**(3/2)))
    

#     if(a==1):
#         if(x0<(math.sqrt(3)-1)):
#             return 1+math.sqrt(2*(1-x0))-x0
#         else :
#             return 1
        
#     # This case will not be allowed.

#     # General Case
#       # M = 1# redundant
#       #  prec = Precision[:a1, x0]
#       # sign = x0/math.abs(x0)
#     #  console.log('Going to call KGPSR(' + a + ',' + sign + ')')
#     req = KerrGeoPhotonSphereRadius(a, math.sign(x0))
#     #  console.log('Going to call KGPSR(' + a + ',' + 0 + ')')
#     rpolar = KerrGeoPhotonSphereRadius(a, 0)
    
#     def Phi(r,a):return -(r*r*r - 3 *M *r*r + a*a *r + a*a *M)/(a *(r - M))
#     def Q(r,a) :return -(r*r*r *(r*r*r - 6 *M *r*r + 9 *M*M *r - 4 *a*a *M))
#         /(a*a *(r - M)*(r - M))
#     def u0Sq(r,a) :return(
#         (a*a - Q(r,a) - Phi(r,a)**2)
#             + math.sqrt(
#                 (a*a - Q(r,a) - Phi(r,a)**2)**2 + 4 *a*a *Q(r,a)
#             )
#     )
#                         /(2 *a*a)
#     def f(r):return 1-u0Sq(r,a)-(x0*x0)
#     # let ans=bisection(f,math.min(req,rpolar),math.max(req,rpolar))
#     # if(ans==0):let too=newtonraphson(f,math.min(req,rpolar),math.max(req,rpolar),(req + rpolar)/2)
#     # return too
#     # console.log("bisection "+ans)
#     # return ans
#     return newtonraphson(f,math.min(req,rpolar),math.max(req,rpolar),(req + rpolar)/2)



# # Innermost Bound Spherical Orbits
# def KerrGeoIBSO(a, x0, E=None ):#  / Abs[x1] <= 1] / # Precision[:a1, x1] != \[Infinity] :=
#     if(a==0):return 4# Approximately Schwarzschild
#     if((x0*x0)==1):return 2-x0*a+2*math.sqrt(1-x0*a)# (*Equatorial IBSO results from Bardeen,Press,Teukolsky 1972*)
#     if(x0==0):if(a==1)# case will not be allowed
#               :return (3 + (54 - 6*math.sqrt(33))**(1/3) + (6* (9 + math.sqrt(33)))**(1/3))/3
#               let delt=26*a*a*a*a-8*a**6+3*math.sqrt(3*(27*a**8-16*a**10))
#               return 1 +
#               math.sqrt(12 - 4* a*a -
#                         (6* math.sqrt(6)* (-2 + a*a))
#                         /math.sqrt(6 - 2* a*a + 4* a*a*a*a/delt**(1/3) + delt**(1/3))
#                         - 4* a*a*a*a/delt**(1/3) - delt**(1/3)
#                        )/math.sqrt(6)
#               + math.sqrt(6 - 2* a*a + 4* a*a*a*a/delt**(1/3) + delt**(1/3))
#               /math.sqrt(6)
#              # Polar case

#     # prec = Precision[:a1, x]

#     # General case
#     rph = KerrGeoPhotonSphereRadius(a, x0)
#     E0 = ru => (E==None)? KerrGeoEnergy(a, ru, 0, x0):E(ru)# E0 could be given if it has already been calculated
#     let n=0
#     while(E0(rph + 10**(-n)) < 1): n++
#     # let f=ru=>E0(ru)-1
#     let ans=bisection(def(ru):return E0(ru)-1,rph+10**(-n),10)
#     return ans
#     # console.log("bisection "+ans)
#     # return newtonraphson(f, rph + 10**(-n), 10, rph+2,15)


# def KerrGeoSeparatrix(a, e1, x):#  / Abs[x1] <= 1] /
#     # Precision[:a1, e1, x1] != \[Infinity] :=
#     if(a==0):return 6+2*e1# Schwarzschild
#     #  spin equal to 1 will not be an allowed case, but here it is anyway.# /
#     if(a==1):if(x==1):return 1+e1# From Glampedakis and Kennefick arXiv : gr - qc/0203086, for a = M we have Subscript[p, s] = 1 + e
#              if(x==0 && e1==0):return 1 + math.sqrt(3) + math.sqrt(3 + 2* Sqrt(3))# Polar ISSO in extremal case
            
#     # # # 
#     if(e1>0.99):return 2*KerrGeoIBSO(a, x)
#     if((x*x)==1):
#         def f(ru,a,e1,x):
#             return e1-
#             (0-ru*ru + 6* ru - 8* x*a* math.sqrt(ru) + 3* a*a)
#                 /(ru*ru - 2* ru + a*a)
        
#         let ribco=2-x*a+2*math.sqrt(1-x*a)# radius of inner-most bound circular orbit. From Levin and Perez-Giz arXiv:0811.3814
#         let rup=bisection(def(ru):return f(ru,a,e1,x),ribco,10)
#         #  console.log(f(ribco)+" "+rup)
#         return (4* rup* (math.sqrt(rup) - a)**2)/(rup*p - 2* rup + a*a)
    
#     C0 = ru =>  KerrGeoConstantsOfMotion(a, ru, 0, x)# :E0, L0, Q0->C0(ru).En .L .Q
#     Beta =  ru => (-1 + C0(ru).En**2) # (-1 + E0^2)
#     ra2 =  ru => 2 *(-a *C0(ru).En + C0(ru).L)**2
#         + 2 *C0(ru).Q
#         + ru*ru * \        (-1 - ru *Beta(ru) +
#          math.sqrt(1 + Beta(ru) * \                   (C0(ru).L**2 + C0(ru).Q - a*a *Beta(ru) - 2 *ru * \                    (1 + ru *Beta(ru))
#                    )
#                   )
#         )
#     Beta2 = (ru) => :let B= Beta(ru)
#                          return ru * \                         (2 + ru *B - 2 * \                          math.sqrt(1 +B* \                                    (C0(ru).L**2 + C0(ru).Q - 2 *ru - a*a *B - 2 *ru*ru *B
#                                     )
#                                    )
#                          )# I simplified this to reduce the number of def calls.
#     e2 = ru => (ra2(ru) - ru *Beta2(ru))/(ra2(ru) + ru *Beta2(ru))
#     #  prec = Precision[:a, e1, x]
#     E =r=>C0(r).En
#     r1 = KerrGeoIBSO(a, x, E)
#     let sol=ru=>e2(ru)-e1
#     ru0 = bisection(sol,r1,10)# newtonraphson(sol, r1, 10, (r1 + 10)/2, 20, undefined, restrain=true)
#     #  console.log(ru0)
#     #  console.log((2 *ra2(ru0) *ru0)/(ra2(ru0) + ru0 *Beta2(ru0)))
#     return (2 *ra2(ru0) *ru0)/(ra2(ru0) + ru0 *Beta2(ru0))
#     # (*This method is an extension of the method in arXiv:1108.1819.See N.Warburton's notes for details.*)


# def KerrGeoBoundOrbitQ(a, p, e, x):return p >= KerrGeoSeparatrix(a, e, x)
