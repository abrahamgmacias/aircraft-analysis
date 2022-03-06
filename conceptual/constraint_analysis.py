import math
import numpy as np
import matplotlib.pyplot as plt

# Parametros propuestos
vs = 12 #Velocidad de stall, m/s
vc = 14 #Velocidad de crucero, m/s
vv = 0.508 #Velocidad vertical en Vy

# Aircraft
mtow = 22 #Peso máximo, kg
Cdmin = 0.035 #Coeficiente de drag mínimo, tabla 3.1 GUDMUNDSSON

# Wing
AR = 7.5 #Relación de aspecto
lam = 1 #Taper, rectangular
e = 0.75 #Factor de eficiencia de Oswald                # Make formula for this 
k = 1/((math.pi)*AR*e)                                  # Make formula for this 

# Motor
p = 900 #Potencia, Watts
n = 2 #Factor de carga = 1/cos(angulo de banqueo)
ep = 0.8 #Eficiencia prop
tw_real = (ep*p)/(mtow*g*vc)

# Flight conditions
rho_sea = 1.225 #densidad del aire a nivel del mar, kg/m^3
rho_desired = 0.974 #densidad del aire en un lugar en específico, kg/m^3
mu = 0.05 #Coeficiente de fricción, tabla 17.1 Raymer, Dry concrete/asphalt
g = 9.807 #Gravedad, m/s^2

# Aerodynamic conds
CLmax = 1.2 #Coeficiente de lift máximo deseado
CLto = 1.2 #Coeficiente de lift deseado en despegue
CDto = 0.04 #Coeficiente de drag deseado en despegue
to_d = 61 #Distancia de despegue, metros

# Range data
ws = np.linspace(30,170,171)

# Exportar la clase a af
# Sacar la plot section y dejarlo en el executable
# Hacer un dict con cada uno de los elementos del init
class Constraints():
    def __init__(self,vs,vc,vv,mtow,Cdmin,AR,lam,p,n,rho_sea,rho_desired,mu,e,CLmax,CLto,CDto,to_d,k,ws,ep,tw_real):
        self.vs = vs 
        self.vc = vc 
        self.vv = vv 
        self.mtow = mtow 
        self.Cdmin = Cdmin 
        self.AR = AR 
        self.lam = lam 
        self.p = p 
        self.n = n
        self.rho_sea = rho_sea 
        self.rho_desired = rho_desired 
        self.mu = mu 
        self.g = g 
        self.e = e 
        self.CLmax = CLmax 
        self.CLto = CLto 
        self.CDto = CDto 
        self.to_d = to_d 
        self.k = k
        self.ws = ws
        self.ep = np
        self.tw_real = tw_real

    #Métodos para encontrar T/W
    def turn(self): #Constant velocity turn
        q_turn = 0.5*(self.rho_sea)*(self.vc**2)
        a = (self.Cdmin)/(self.ws)
        b = (self.n)/(q_turn)
        return q_turn*(a + (self.k)*(b**2)*(self.ws))
    
    def roc(self): #Rate of Climb
        vy = ((2/self.rho_sea)*self.ws*np.sqrt(self.k/(3*self.Cdmin)))**0.5
        q_climb = 0.5*self.rho_sea*(vy**2)
        return (self.vv/vy + (q_climb/self.ws)*self.Cdmin) + (self.k/q_climb)*self.ws

    def takeoff(self): #Desired Takeoff Distance
        vto = 1.2*np.sqrt(self.ws*(2/(self.rho_sea*self.CLmax)))
        q_takeoff = 0.5*self.rho_sea*(vto**2)
        return (vto**2)/(2*self.g*self.to_d) + (q_takeoff*self.CDto)/self.ws + self.mu*(1-(q_takeoff*self.CLto)/self.ws)

    def cruise(self): #Desired cruise Airspeed
        q_cruise = 0.5*self.rho_sea*(vc**2)
        return q_cruise*self.Cdmin*(1/self.ws) + self.k*(1/q_cruise)*self.ws
    
    def CLmax1(self):
        q_stall1 = 0.5*self.rho_sea*(self.vs**2)
        return (1/q_stall1)*self.ws

    def CLmax2(self):
        q_stall2 = 0.5*self.rho_sea*((self.vs-2)**2)
        return (1/q_stall2)*self.ws

    def CLmax3(self):
        q_stall3 = 0.5*self.rho_sea*((self.vs+2)**2)
        return (1/q_stall3)*self.ws

    def CLmax4(self):
        q_stall4 = 0.5*self.rho_sea*((self.vs+4)**2)
        return (1/q_stall4)*self.ws

    def CLmax5(self):
        q_stall5 = 0.5*self.rho_sea*((self.vs-4)**2)
        return (1/q_stall5)*self.ws


def plot(self):
    fig, ax1 = plt.subplots()
    ax1.plot(self.ws, self.turn(), color = 'black',label='Constant velocity turn')
    ax1.plot(self.ws, self.roc(), color = 'red',label='Rate of Climb')
    ax1.plot(self.ws, self.takeoff(), color = 'blue',label='Desired Takeoff Distance')
    ax1.plot(self.ws, self.cruise(), color = 'yellow',label='Desired cruise Airspeed')
    ax1.hlines(self.tw_real, 30, 170, colors='k', linestyles='dashed',label='T/W Real posible con Vcrucero')
    ax2 = ax1.twinx()
    ax2.set_ylabel('CLmax',fontsize=25)
    ax2.plot(self.ws, self.CLmax1(), color = 'k', linestyle='-.',label='Recta Vstall = 12 m/s')
    ax2.plot(self.ws, self.CLmax2(), color = 'b', linestyle='-.',label='Recta Vstall = 10 m/s')
    ax2.plot(self.ws, self.CLmax3(), color = 'r', linestyle='-.',label='Recta Vstall = 14 m/s')
    ax2.plot(self.ws, self.CLmax4(), color = 'magenta', linestyle='-.',label='Recta Vstall = 16 m/s')
    #ax2.plot(self.ws, self.CLmax5(), color = 'green', linestyle='-.',label='Recta Vstall = 8 m/s')
    plt.title('Constraint Diagram', fontsize = 25,fontweight='bold')
    ax1.set_ylabel('T/W', fontsize = 25)
    ax1.set_xlabel(' W/S (N/m^2)', fontsize = 25)
    ax1.legend()
    ax2.legend(loc='upper center')
    fig.tight_layout()
    plt.show()
