import math

import numpy as np
import numpy.linalg as nl

import aircraft_funcs as af
import aircraft_data 

from scipy.optimize import fsolve as fsolve

# Data
S_sa = 0.166                                                                   # Side area - Booms, fuselage - m2
Stot = 1.02*(S_sa + Sv/2)                                                      # Total projected side area - m2

x_sa_ca = 0.6153                                                               # Centro geometrico S_sa - m
x_ca = (Stot*x_sa_ca+(Sv/2)*2.029)/(Stot + Sv/2)                               # Centro geometrico Stot - m

# Posiciones - Con respecto a la nariz
h_0x_position = [0.47695, 0.51249]                                             # Posiciones del CG - dx_h0_list
dx_cg_h_list = [1.4783, 1.4425]                                                # Brazo de palanca - Most forward, most aft CG - m

dc = x_ca-h_0x_position[1]                                                     # Distancia entre el CG y el centro geometrico - m

dσ_dB = 0                                                                      # Vertical Tail Sidewash Gradient
Cy0 = 0                                                                        # Simetrico
Cn0 = 0                                                                        # Simetrico

Sr_max = 24                                                                    # Rudder def. max. - deg

rho = 1.225
Vs = af.Vel(mtow*9.81, CL_max, f_flaps_list[0], rho, Sw)                       # Stall velocity - m/s
Vw = 1                                                                         # Crosswind velocity - m/s
Vapp = 1.1*Vs                                                                  # Approach velocity - m/s
Vt = ((Vapp**2)+(Vw**2))**0.5                                                  # Total velocity - m/s

Cdy = 0.5                                                                      # Side drag coef.
Vv = 0.05                                                                      # Vertical volume Coef.
br_bv, cr_cv = 1, 0.42                                                         # Table 12.3 - Cessna 182 (pg. 12.20)
Kf1, Kf2 = 0.65, 1.3                                                           # Sect. 12.6.2.3


# Calculo del Rudder - Mohammed 12.8.3 - pg. 738
# Criterio Crosswing Landing
class Rudder():
    def __init__(self, rho, Vapp, Vw, Vt, Stot, Sv, Sw, bw, br_bv, bv, cr_cv, cv, dσ_dB, Vv, nv, Cdy, Cn0, Cy0,
                 Kf1, Kf2, Sr_max, av, dx_cg_h_list, dc):
        self.rho = rho
        self.Vapp = Vapp
        self.Vw = Vw
        self.Vt = Vt
        
        self.Stot = Stot
        self.Sv = Sv
        self.Sw = Sw
        self.bw = bw
        
        self.br_bv = br_bv
        self.bv = bv
        self.cr_cv = cr_cv
        self.cv = cv
        
        self.dσ_dB = dσ_dB
        self.Vv = Vv
        self.nv = nv
        
        self.Cdy = Cdy
        self.Cn0 = Cn0
        self.Cy0 = Cy0
        
        self.Kf1 = Kf1
        self.Kf2 = Kf2
        self.Sr_max = Sr_max
        self.av = av
        
        self.dx_cg_h_list = dx_cg_h_list
        self.dc = dc
 
    # Side wind force
    def SideWindForce(self):
        return 0.5*self.rho*(self.Vw**2)*self.Stot*self.Cdy                                  
    
    # Sideslip angle
    def SideSlipAngle(self):
        return math.degrees(math.atan(self.Vw/self.Vapp))                             
    
    # Sideslip derivative
    def CNBDeriv(self):
        return self.Kf1*self.av*(1 - self.dσ_dB)*self.nv*(self.dx_cg_h_list[1]*self.Sv/(self.bw*self.Sw))      
    
    # Sideslip derivative
    def CYBDeriv(self):
        return -self.Kf2*self.av*(1 - self.dσ_dB)*self.nv*(self.Sv/self.Sw)                         
    
    # Control derivative
    def CYSRDeriv(self):
        return self.av*self.nv*self.tau_r*self.br_bv*(self.Sv/self.Sw)                              
    
    # Control derivatives
    def CNSRDeriv(self):
        return -self.av*self.Vv*self.nv*self.tau_r*self.br_bv                                   
    
    # Nonlinear System of Equations - x = Sr, y = σ 
    def hg(self, xy):
        x = xy[0]
        y = xy[1]
        
        h = self.elem1*self.Cn0 + self.elem1*self.Cnb*(math.radians(self.β)) - self.elem1*self.Cnb*x + self.elem1*self.Cn_Sr*y + self.Fw*self.dc*np.cos(x)
        g = self.elem2*self.Cy0 +  self.elem2*self.Cyb*(math.radians(self.β)) - self.elem2*self.Cyb*x + self.elem2*self.Cy_Sr*y - self.Fw
        return np.array([h, g])
    
    def RudderSolver(self):
        self.Fw = self.SideWindForce()
        self.β = self.SideSlipAngle()
        self.Cnb = self.CNBDeriv()
        self.Cyb = self.CYBDeriv()
        self.tau_r = af.TauForCxCy(self.cr_cv, None)
        self.Cy_Sr = self.CYSRDeriv()
        self.Cn_Sr = self.CNSRDeriv()
        
        self.elem1 = 0.5*self.rho*(self.Vt**2)*self.Sw*self.bw
        self.elem2 = self.elem1 / self.bw 
        
        xy0 = np.array([1, 1])
        xy = fsolve(self.hg, xy0)                                                           
        xy = [n*57.3 for n in xy]
        
        print(xy)
        print(self.Cy_Sr, self.Cn_Sr)
        
        if xy[0] < self.Sr_max:
            cr = self.cr_cv*self.cv
            br = self.br_bv*self.bv / 2
            Sr = br*cr / 2 
            return [cr, br, Sr]
        
        return print("SR max. has been surpassed")
    




x = Rudder(rho, Vapp, Vw, Vt, Stot, Sv, Sw, bw, br_bv, bv, cr_cv, cv, dσ_dB, Vv, nv, Cdy, Cn0, Cy0, Kf1, Kf2, Sr_max, av, dx_cg_h_list, dc)
rudder_data = x.RudderSolver()
print(rudder_data)


