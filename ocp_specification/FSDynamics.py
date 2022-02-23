import casadi as cas
import numpy as np


class FSDynamics:
    def __init__(self):
        self.indR = range(9)
        self.indp = range(9, 12)

    def geo_integrator_tra(R_t, p_obj, u, h):
        """Integrate invariants over interval h starting from a current state (object pose + moving frames)"""
        # Define a geometric integrator for eFSI,
        # (meaning rigid-body motion is perfectly integrated assuming constant invariants)

        # object translation speed
        # curvature speed translational Frenet-Serret
        # torsion speed translational Frenet-Serret

        i4 = u[0]
        i5 = u[1]
        i6 = u[2]

        omega = cas.vertcat(i6, i5, 0)
        omega_norm = cas.norm_2(omega)
        v = cas.vertcat(i4, 0, 0)

        deltaR = np.eye(3) + cas.sin(omega_norm @ h)/omega_norm*cas.skew(omega) + \
            (1-cas.cos(omega_norm @ h))/omega_norm**2 * \
            cas.mtimes(cas.skew(omega), cas.skew(omega))
        deltaP = (np.eye(3)-deltaR) @ cas.skew(omega) @ v / \
            omega_norm**2 + omega @ omega.T @ v/omega_norm**2*h

        R_t_plus1 = R_t @ deltaR
        p_obj_plus1 = R_t @ deltaP + p_obj

        return (R_t_plus1, p_obj_plus1)

    def dynamics(self, uk, xk, dt):
        Rk = xk[self.indR].reshape(3, 3).T
        Rkp1, pkp1 = self.geo_integrator_tra(Rk, uk, dt)
        xkp1 = cas.SX.zeros(12, 1)
        xkp1[self.indR] = (Rkp1.T).reshape(9, 1)
        xkp1[self.indp] = pkp1
        return xkp1
