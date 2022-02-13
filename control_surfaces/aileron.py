import aircraft_data
import aircraft_funcs as af

def Qdyn(u0, rho):
    return 0.5*rho*u0**2
    
def CLP(aw, lamdaw):
    return -aw*((1+3*lamdaw)/(1+lamdaw))/12

def LPfunc(Q, u0, Sw, bw, Clp, Ixx):
    return (Q*Sw*Clp*bw**2)/(2*Ixx*u0)

def LSafunc(Cl_Sa, Q, Sw, bw, Ixx):
    return Cl_Sa*Q*Sw*bw/Ixx

def PSS(L_Sa, Lp, Sa):
    return -L_Sa*Sa/Lp
    
def CLSA(aw, tau, cw, Sw, bw, lamdaw, y):
    return (2*aw*tau*cw/(Sw*bw))*(((y**2)/2)+((lamdaw-1)/bw)*(2*y**3)/3) 
    
y_limits = [0.61, 0.95]
y_lengths = [bw*y/2 for y in y_limits]
Cl_Sa = [CLSA(aw, 0.45, cw, Sw, bw, 1, y) for y in y_lengths]
Cl_Sa = Cl_Sa[1]-Cl_Sa[0]

V = 14
Ixx = Ifull_list[0]
Q = Qdyn(V, rho_list[0])
Clp = CLP(aw, lamdaw)
Lp = LPfunc(Q, V, Sw, bw, Clp, Ixx)
tau = 1/Lp
Cl_Sa = 0.3156
Sa = 5

L_Sa = LSafunc(Cl_Sa, Q, Sw, bw, Ixx)
Pss = PSS(L_Sa, Lp, 5)/57.3


print(Clp, Q, Sw, bw, Ifull_list[0])
print(L_Sa, Lp, Cl_Sa, Pss)