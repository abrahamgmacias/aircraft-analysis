import math
import aircraft_data 
import aircraft_funcs as af

# Matplotlib
import matplotlib.pyplot as plt

# Aircraft Data
lf = 0.9        # Longitud del fuselaje - m
h_cg = 0.2      # Altura del CG al suelo - m

# Aceleraciones - m / s2
a_TO = 2.5035
a_LA = -0.2912
μ = 0.05

# Posiciones - Con respecto a la nariz
margen = 0.05                               # Margen de movimiento del CG
h_0x_position = [0.476, 0.520]              # Posiciones del CG
ng_x = 0.1                                  # Posicion del Nose Gear
Nose_100 = 15                               # Porcentaje de carga total por el Nose Gear

# FOR SECOND FUNCTION
# Posiciones - Con respecto a la nariz
dx_cg_h_list = [1.48, 1.436]                                                   # Brazo de palanca - Most forward, most aft CG - m
dx_acw = 0.4535                                                                # Centro aerodinamico del ala - m 
dx_ach = 2.0115                                                                # Centro aerodinamico del empenaje - m

# Posiciones - Con respecto al centroide del fuselaje
dy_w_ac = 0.1                                                                  # Centro aerodinamico del ala - m  
dy_h_ac = dy_w_ac                                                              # Centro aerodinamico del empenaje
dy_h0 = 0.045                                                                  # Centro de gravedad                                           

# Empennage
be_bh = 1                                                                      # Span-to-tail ratio - From Table 12.3
Se_max = -20                                                                   # Deflexion max. - deg - From Table 12.3

delta_alpha_he = 8.5                                                           # Table 
graph_actuator = 1

# Yada
alpha_s = 14                                                                   # Wing stall angle - deg 
alpha_s_to = alpha_s - iw                                                      # Aircraft stall angle - deg 
alpha_TO = alpha_s_to - 2                                                      # Aircraft max. TO angle - deg 


# -------------------------------------------- Classes / Functions ------------------------------------------ #
class AircraftSubsystems():
    def __init__(self, mtow, g_accl, lf, h_cg, a_TO, a_LA, margen, ng_x, Nose_100):
        self.mtow = mtow
        self.g_accl = g_accl
        self.lf = lf
        self.h_cg = h_cg
        self.a_TO = a_TO
        self.a_LA = a_LA
        self.margen = margen
        self.h_0x_position = h_0x_position
        self.ng_x = ng_x
        self.Nose_100 = Nose_100

    def LandingGear(self):
        # Ex. 9.3 - Aircraft Design: A Systems Engineering Approach
        # Distancia entre el Nose Gear (Bn) y el CG
        Bn_positions = [h_0x - self.ng_x for h_0x in h_0x_position]           
        Bn_min, Bn_max = min(Bn_positions), max(Bn_positions)          
                                                    
        # Wheelbase
        B = Bn_min*(100/self.Nose_100)/(100/self.Nose_100-1)
        
        # Distancia entre el CG y el Main Gear (Bm)
        Bm_positions = [B - Bn for Bn in Bn_positions]
        Bm_min, Bm_max = min(Bn_positions), max(Bn_positions)          
        
        # Carga Max en Bn - Landing Brake
        Fn = self.mtow*self.g_accl*(Bm_max / B) + ((self.mtow*self.g_accl*abs(self.a_LA)*self.h_cg) / self.g_accl*B)
        Load_N = ((Fn / self.g_accl) / self.mtow)*100
        
        # Carga Max en Bm - Take-Off
        Fm = self.mtow*self.g_accl*(Bn_max / B) + ((self.mtow*self.g_accl*self.a_TO*self.h_cg) / self.g_accl*B)
        Load_M = ((Fm / self.g_accl) / self.mtow)*100
        
        # Porcentage de Carga sobre Nose/Main Gear
        Load_perc = [Load_N, Load_M]
        
        # Validacion
        l_margen = self.lf - (self.ng_x + B) 
        if l_margen < 0:
            print("Main Gear exceeds the limits of the fuselage")
            
        return Bn_positions, Bm_positions, round(B,4), Load_perc

class LongitudinalControlSurfaces():
    def __init__(self, g_accl, Iyy, mtow, f_flaps, CL_max, CL_min, rho, μ, aw, Sw, cw, iw, ARw, ew, CL0, CD0, CM0,
                 Tmax, ah, Sh, bh, ih, nh, Vh, B, h_cg, ng_x, dy_h0, dx_acw, dx_ach, Bm_positions, h_0x_position, Se_max, be_bh,
                 lf, a_TO, a_LA, margen, Nose_100, delta_alpha_he, alpha_TO, alpha_hs, graph_actuator):
        self.g_accl = g_accl
        self.Iyy = Iyy
        self.mtow = mtow
        
        self.f_flaps = f_flaps
        self.CL_max = CL_max
        self.CL_min = CL_min
        self.rho = rho
        self.μ = μ
    
        self.aw = aw
        self.Sw = Sw
        self.cw = cw
        self.iw = iw
        self.ARw = ARw
        self.ew = ew
        
        self.CL0 = CL0 
        self.CD0 = CD0
        self.CM0 = CM0
        self.Tmax = Tmax
        
        self.ah = ah
        self.Sh = Sh
        self.bh = bh
        self.ih = ih
        self.nh = nh
        self.Vh = Vh
        
        self.B = B                                                             #
        self.h_cg = h_cg                                                       #
        self.ng_x = ng_x                                                       #
        self.dy_h0 = dy_h0
        self.dx_ac_w = dx_acw                                                  # 
        self.dx_ac_h = dx_ach                                                  # 
        self.Bm_positions = Bm_positions                                               # Distance from CG to Main Gear
        self.h_0x_position = h_0x_position
    
        self.Se_max = Se_max
        self.be_bh = be_bh
        
        self.lf = lf
        self.a_TO = a_TO
        self.a_LA = a_LA
        self.margen = margen
        self.Nose_100 = Nose_100
        
        self.delta_alpha_he = delta_alpha_he
        self.alpha_TO = alpha_TO
        self.alpha_hs = alpha_hs
        
        self.graph_actuator = graph_actuator
        
    def Forces(self, Vs):
        # Forces - Take Off
        Vr = Vs
        K = 1 / (math.pi*self.ew*self.ARw)
        CD0_TO = self.CD0 + K*self.CL0*self.f_flaps**2
        
        D_TO = af.Lift(self.rho, CD0_TO, self.Sw, Vr)
        L_TO = af.Drag(self.rho, self.CL0*self.f_flaps, self.Sw, Vr)
        Ff = self.μ*(self.mtow*self.g_accl - L_TO)
        a = (self.Tmax - D_TO - Ff) / self.mtow
        
        return Vr, D_TO, L_TO, Ff, a
        
    def Moments(self, Vr, D_TO, L_TO, a, Cm_acw):
        # Distances for moments - Con respecto al main gear - Take Off
        dx_lw_mg = self.dx_ac_w-self.B+self.ng_x                               # Wing ac - m
        dy_lh_mg = self.dx_ac_h-self.B+self.ng_x                               # Empennage ac - m   
        dy_T_mg, dy_D_mg = self.h_cg + self.dy_h0, self.h_cg + self.dy_h0      # Motor, drag - m
        
        # Moment - Take Off - Most forward CG (N*m)
        M_acw = af.Moment(self.rho, Cm_acw, self.Sw, self.cw, Vr)              # N*m
        Mw = -self.mtow*self.g_accl*self.Bm_positions[0]                        # N*m
        Md = D_TO*dy_D_mg                                                      # N*m
        Mt = -self.Tmax*dy_T_mg                                                # N*m
        Ml_w = L_TO*dx_lw_mg                                                   # N*m
        Ma = self.mtow*a*self.h_cg   
        
        return dx_lw_mg, dy_lh_mg, M_acw, Mw, Md, Mt, Ml_w, Ma
    
    def CLiftDesired(self, Vs, dy_lh_mg, M_acw, Mw, Md, Mt, Ml_w, Ma):
        # Desired horizontal tail lift and coef. - Take Off
        Lh = (Ml_w + M_acw + Ma + Mw + Md + Mt - (self.Iyy*(self.mtow*self.g_accl/57.3))) / dy_lh_mg
        Clh = af.CLift(Lh, self.rho, self.f_flaps, Vs, self.Sh)
        
        return Lh, Clh

    def EffectiveAngle(self, Clh):
        # Figure 12.12 - Control surface angle of attack effectiveness parameter
        tau_ce_ch = [[0, 0.02, 0.03, 0.04, 0.05, 0.05, 0.06, 0.07, 0.08, 0.09, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.14, 0.15, 
                        0.16, 0.17, 0.18, 0.19, 0.2, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.27, 0.28, 0.29, 0.3, 0.31, 
                        0.32, 0.33, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.46, 
                        0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.59, 0.6, 0.61, 0.62, 
                        0.63, 0.64, 0.65, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.01, 0.02, 0.04, 0.01], [0, 0.08, 0.09, 0.11, 0.15, 0.17,
                        0.18, 0.2, 0.22, 0.23, 0.25, 0.26, 0.28, 0.29, 0.3, 0.32, 0.33, 0.34, 0.35, 0.36, 0.38, 0.39, 0.39, 0.4, 0.41, 
                        0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.49, 0.5, 0.51, 0.52, 0.53, 0.53, 0.54, 0.55, 0.55, 0.56, 0.57,
                        0.57, 0.58, 0.59, 0.6, 0.6, 0.61, 0.61, 0.62, 0.63, 0.63, 0.64, 0.65, 0.65, 0.66, 0.67, 0.67, 0.68, 0.69, 0.69, 
                        0.7, 0.7, 0.71, 0.71, 0.72, 0.73, 0.73, 0.74, 0.74, 0.75, 0.76, 0.76, 0.77, 0.77, 0.78, 0.79, 0.79, 0.8, 0.02,
                        0.06, 0.13, 0.04]]
        
        # Angle of effectiveness of the elevator
        E_0 = 2*self.CL0*57.3 / (math.pi*self.ARw)                                      # deg
        delta_E = 2*self.aw / (math.pi*self.ARw)                                  # deg/deg
        E = E_0 + delta_E*self.iw
        alpha_h = self.iw + self.ih - E                                                  # deg
    
        tau_e = ((alpha_h/57.3) + (Clh / self.ah)) / (self.Se_max / 57.3)
        ce_ch = af.Closest(tau_ce_ch[1], tau_ce_ch[0], tau_e)   
        
        # Change in tail lift coef. for elevator deflection
        delta_alpha_e = -1.15*ce_ch*self.Se_max
        
        return tau_e, delta_alpha_e, delta_E, E_0, ce_ch
    
    def EffectiveDerivs(self, v, dx_h0, dx_lw_mg, dy_lh_mg, tau_e, delta_E):
        # Elevator effectiveness derivatives
        Cm_delta_e = -self.ah*self.nh*self.Vh*self.be_bh*tau_e
        Cl_delta_e = self.ah*self.nh*(self.Sh/self.Sw)*self.be_bh*tau_e
        Clh_delta_e = self.ah*tau_e
        Cm_alpha = self.aw*((dx_lw_mg-dx_h0)/self.cw) - self.ah*self.nh*(self.Sh/self.Sw)*((dy_lh_mg+dx_h0)/self.cw)*(1-delta_E)
        
        q = af.DynPressure(self.rho, v)
        Cl_1 = af.CLift(self.mtow*self.g_accl, self.rho, self.f_flaps, v, self.Sw)
        Delta_e = -((((self.Tmax*self.dy_h0/(q*self.Sw*self.cw))+self.CM0)*self.aw + (Cl_1-self.CL0)*Cm_alpha)/(self.aw*Cm_delta_e-Cm_alpha*Cl_delta_e))*57.3

        return Delta_e
    
    def Verification(self, delta_E, E_0, ce_ch):
        # Stall check @ TO Rotation
        alpha_h_TO = self.alpha_TO*(1 - delta_E) + self.ih - E_0
        alpha_h_s_min, alpha_h_s_max = -(self.alpha_hs-self.delta_alpha_he), (self.alpha_hs-self.delta_alpha_he)
        
        # Elevator Geometry
        if abs(alpha_h_s_min) > alpha_h_TO:
            be = self.bh 
            ch = self.Sh / self.bh 
            ce = ce_ch*ch
            Se = be*ce 
            
        return be, ce, Se, alpha_h_s_min, alpha_h_s_max
    
    def Elevator(self):
        Vs = af.Vel(self.mtow*self.g_accl, self.CL_max, self.f_flaps, self.rho, self.Sw)
        Vmax = af.Vel(self.mtow*self.g_accl, self.CL_min, self.f_flaps, self.rho, self.Sw)
        V_list = [v*0.25 for v in range(int(Vs)*4, int(Vmax)*4+1)]
        Cm_acw = self.CM0
        
        Vr, D_TO, L_TO, Ff, a = self.Forces(Vs)
        dx_lw_mg, dy_lh_mg, M_acw, Mw, Md, Mt, Ml_w, Ma = self.Moments(Vr, D_TO, L_TO, a, Cm_acw)
        Lh, Clh = self.CLiftDesired(Vs, dy_lh_mg, M_acw, Mw, Md, Mt, Ml_w, Ma)
        tau_e, delta_alpha_e, delta_E, E_0, ce_ch = self.EffectiveAngle(Clh)
        
        Delta_e_list = [[], []]
        for dx_h0 in self.Bm_positions:
            for v in V_list:
                Delta_e = self.EffectiveDerivs(v, dx_h0, dx_lw_mg, dy_lh_mg, tau_e, delta_E)
                Delta_e_list[self.Bm_positions.index(dx_h0)].append(Delta_e)
            plt.plot(V_list, Delta_e_list[self.Bm_positions.index(dx_h0)], '-o')
        
        # Grafica - Elevator deflection vs. velocity   
        if self.graph_actuator == 1:
            plt.grid("On")
            plt.title('Elevator deflection vs. Velocity')
            plt.xlabel('Velocity (m/s)')
            plt.ylabel ('Elevator deflection (deg)')
            plt.show()
        
        delta_e_max = max([max(leest) for leest in Delta_e_list])
        delta_e_min = min([min(leest) for leest in Delta_e_list])
        
        be, ce, Se, alpha_h_s_min, alpha_h_s_max = self.Verification(delta_E, E_0, ce_ch)
        
        return round(be,4), round(ce,4), round(Se,4), round(delta_e_max,4), round(delta_e_min,4)


# ------------------------------------------------- Solvers ------------------------------------------------- #
x = AircraftSubsystems(mtow, g_accl, lf, h_cg, a_TO, a_LA, margen, ng_x, Nose_100)
Bn_positions, Bm_positions, B, Load_perc = x.LandingGear()

y = LongitudinalControlSurfaces(g_accl, Ifull_list[1], mtow, f_flaps_list[0], CL_max, CL_min, rho_list[0], μ, aw, Sw, cw, iw, 
                    ARw, ew, CL0, CD0, CM0, T_stat_N, ah, Sh, bh, ih, nh, Vh, B, h_cg, ng_x, dy_h0, dx_acw, dx_ach, 
                    Bm_positions, h_0x_position, Se_max, be_bh, lf, a_TO, a_LA, margen, Nose_100, delta_alpha_he, 
                    alpha_TO, alpha_hs, graph_actuator)

be, ce, Se, delta_e_max, delta_e_min = y.Elevator()
print(Bm_positions)
print("B: {a}".format(a=B))
print("be: {a}, ce: {b}, Se_max: {c}, Se_min: {d}".format(a=be,b=ce,c=delta_e_max,d=delta_e_min))
print(1-ce/ch)
