#EoMs

# acceleration^mu = - (Gamma + dGamma)^mu_ab v^a v^b

# rungeKutta to for vdot and v
def vdot_func(TotalGamma):
    def vdot(zeta,vel):
        return np.einsum('ijk,j,k->i',TotalGamma[i, j, k],vel[j],vel[k])
    return vdot

# rk4 again for xdot and x
def xdot_func():
    def xdot(zeta,xdot0):
        vel=Integrator.rungeKutta(Integrator,xdot0, vdot_func, zeta, 1)
        return vel
    return xdot