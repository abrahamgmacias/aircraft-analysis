import math

# General Aircraft Functions
class Aerodynamics():
    # Aircraft Dynamic Pressure
    @staticmethod
    def dynamicPressure(rho, v):
        return 0.5*rho*v**2

    # Reynolds Number
    @staticmethod
    def numReynolds(rho, v, MAC, miu):
        return rho*v*MAC / miu

    # Velocity / Stall Velocity / Vuelo Nivelado
    @staticmethod
    def velocity(CL, w, rho, S, f_flaps=1):
        return (2*w / (CL*f_flaps*rho*S))**0.5

    # Velocity / Stall Velocity / Giro Nivelado
    @staticmethod
    def stallVelocity(CL, n_factor, w, rho, S):
        return (2*w*n_factor / (CL*rho*S))**0.5

    # Lift Coefficient / Vuelo Niveglado / Take Off / Touchdown
    @staticmethod
    def liftCoefficient(rho, vel, w, S, f_flaps=1):
        return 2*w / (rho*f_flaps*(vel**2)*S) 

    # Lift Force / Vuelo Nivelado / Take Off / Touchdown
    @staticmethod
    def liftForce(CL, rho, S, v):
        return 0.5*rho*CL*S*(v**2)

    # Drag Coefficient / Vuelo Nivelado
    @staticmethod
    def dragCoefficient(CL, CD_0, e, AR):
        return CD_0 + (CL**2) / (math.pi*e*AR)

    # Drag Curve Slope
    @staticmethod
    def dragCurveSlope(CL, e, ARw, aw):
        return (2*CL / (math.pi*e*ARw))*aw      

    # Thrust Required to maintain Level Flight
    @staticmethod
    def thrustRequired(mass, CL, CD):
        return mass / (CL / CD)



# Aircraft class
class Aircraft():
    def __init__(self):
        self.aeroCoefficients = {}

        # Wing, empenage, landing gear, etc. 
        self.components = {}

        # Empty weight, MTOW, etc.
        self.weights = {}

    # Coefficients must be arrays
    def addCoefficients(self, **kwargs):
        self.aeroCoefficients.update(kwargs)
        return f"'{list(kwargs.keys())[-1]}' was added to aeroCoefficients"

    # Components should be classes 
    def addComponents(self, **kwargs):
        self.components.update(kwargs)
        return f"'{list(kwargs.keys())[-1]}' was added to components" 

    def addWeights(self, **kwargs):
        self.weights.update(kwargs)
        return f"'{list(kwargs.keys())[-1]}' was added to weights" 

    def getCoefficients(self):
        return self.aeroCoefficients


class Wing():
    def __init__(self, bw, cw, airfoil):
        self.geometry = {'bw': bw, 'cw': cw, 'airfoil': airfoil}
        self.wingCoefficients = {}

    def simpleWingArea(self):
        self.Sw = self.geometry['bw']*self.geometry['cw']

    def addGeometry(self, kwargs):
        self.geometry.update(kwargs)
        return f"'{list(kwargs.keys())[-1]}' parameter was added to wing geometry" 

    def addCoefficients(self, kwargs):
        self.wingCoefficients.update(kwargs)
        return f"'{list(kwargs.keys())[-1]}' was added to wing coefficients" 




sae = Aircraft()
sae.addWeights(mtow=22)
sae.addCoefficients(cl_list=[], cd_list=[], cd_0=0.03324, cl_min_d=0, cl_0=0.6096)

wing = Wing('saeWing', 3.71, 0.494, 'MH-114')
wing.addGeometry(ARw=7.5, Sw=1.831)
wing.addCoefficients(ew=0.825)