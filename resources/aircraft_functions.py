import pandas as pd
import math

# General Aircraft Functions
class Aerodynamics():
    # Aircraft Dynamic Pressure
    @staticmethod
    def dynamicPressure(rho, v):
        return 0.5*rho*v**2

    # Reynolds Number
    @staticmethod
    def numReynolds(rho, v, mac, miu):
        return rho*v*mac / miu

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

    # Power Required to maintain Level Flight
    @staticmethod
    def powerRequired():
        pass



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

    def getComponents(self, componentRequested):
        for component in self.components:
            if component == componentRequested:
                return self.components[componentRequested]
            else:
                raise NameError(f"No element named '{componentRequested}'...")


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



from inspect import getsourcefile
from os.path import abspath

def locate_dir(folder_name, root=None):
    # does not handle length error...
    curr_dir = abspath(getsourcefile(lambda:0))

    if root and type(root) == int:
        split = curr_dir.split("\\")[0:-root]
        return '\\'.join(split) + f'\\{folder_name}'

    elif root == None:
        return curr_dir + f'\\{folder_name}'

    else:
        raise TypeError('Root location must be int.')



