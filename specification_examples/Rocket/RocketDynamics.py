from casadi import *
class RockDyns:
   def __init__(self):
       self.ind_pos = list(range(0,2))
       self.ind_rot = list(range(2,6))
       self.ind_vel = list(range(6,8))
       self.ind_omega = [8]
   def dynamics(self, uk, xk, dt = 0.1):
       # split up state into its components position, invariants and orientation (cf. above)
       nx = 9
       self.xk_pos = xk[self.ind_pos]
       self.xk_rot = reshape(xk[self.ind_rot], 2,2)
       self.xk_vel = xk[self.ind_vel]
       self.xk_omega = xk[self.ind_omega]
       xkp1 = MX.zeros(nx)
       drot = MX.zeros(2,2)
       drot[0,0] =cos(self.xk_omega*dt)
       drot[0,1] = -sin(self.xk_omega*dt)
       drot[1,0] = sin(self.xk_omega*dt)
       drot[1,1] = cos(self.xk_omega*dt)
       xkp1[self.ind_rot] = reshape(drot@self.xk_rot, 4,1)
       xkp1[self.ind_pos] = self.xk_pos + self.xk_vel*dt
       xkp1[self.ind_vel] = self.xk_vel + self.xk_rot[:,1]*uk[0]*dt
       xkp1[self.ind_omega] = self.xk_omega + uk[1]*dt
       return xkp1
