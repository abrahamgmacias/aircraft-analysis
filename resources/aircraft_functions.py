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


# English System, correct at the end
class Propeller():
    coefficient_location = {'thrust_coeffs': 'thrust_coefficient_', 'power_coeffs': 'power_coefficient_'}

    def __init__(self, diameter, pitch):
        self.diameter = diameter
        self.pitch = pitch
        self.pitch_diameter = pitch/diameter

    def getPropellerData(self, prefix):
        return pd.read_csv(f'{prefix}{self.pitch_diameter}.csv')

    # If input: Cp / Ct - output: J (And viceversa) - Handles Parabolas
    def getCoefficientRatio(self, value, prefix, Cp=False, Ct=False, J=False):
        coefficient_data = self.get_propeller_data(Propeller.coefficient_location[prefix])
        x_points, y_points = coefficient_data['X'], coefficient_data['Y']

        if Cp == True or Ct == True:
            residuals = [abs(value-y) for y in y_points]
            residuals_copy = residuals.copy()

            if Cp == True:
                residuals.sort(reverse=False)
                x_selected_points_position = [residuals_copy.index(y) for y in residuals[0:2]]
                xy_selected_point = [x_points[y_position] for y_position in x_selected_points_position]

            if Ct == True:
                xy_selected_point = x_points[residuals.index(min(residuals))]

        if J == True:
            residuals = [abs(value-x) for x in x_points]
            xy_selected_point = y_points[residuals.index(min(residuals))]
        
        return xy_selected_point

    def propEfficiency(self, Ct, Cp, advance_ratio):
        nProp = Ct*advance_ratio / Cp
        return nProp

    # Check 550 HP - De dónde sale?
    def cpFactored(self, Ps, air_density, rpm):
        Cp = 550*(Ps / (air_density*((rpm/60)**3)*((self.diameter/12)**5)))
        return Cp

    def ctFactored(self, thrust, air_density, rpm):
        Ct = thrust / (air_density*(rpm**2)*(self.diameter**4))
        return Ct


# rpm_voltage - rpm to voltage ratio
class Motor():
    def __init__(self, propeller, max_rpm=None, max_voltage=None, max_amperage=None, rpm_voltage=None):
        self.propeller = propeller
        self.max_rpm = max_rpm
        self.max_voltage = max_voltage
        self.max_amperage = max_amperage
        self.rpm_voltage = rpm_voltage

        self.motor_data = self.set_motor_data('motor_data.csv')

    def setMotorData(self, file):
        self.motor_data = pd.read_csv(file)

    def getMotorData(self):
        return self.motor_data

    def simpleMaxPower(self):
        simple_power = self.max_voltage*self.max_amperage
        return round(simple_power, 4)

    def complexMaxPower(self):
        complexPower = self.max_rpm*(1/self.rpm_voltage)*self.max_amperage
        return complexPower

    def thrust(self, Ct, air_density, rpm):
        T = Ct*air_density*((rpm/60)**2)*(self.propeller.diameter/12)**4
        return T

    def motorEfficiency(self, motor_input):
        coefficient_data = pd.read_csv('motor_efficiency.csv')
        x_points, y_points = coefficient_data['X'], coefficient_data['Y']
        residuals = [motor_input-x for x in x_points if motor_input-x >= 0]
        y_selected_point = y_points[residuals.index(min(residuals))]
        return y_selected_point

    def powerS(self, Cp, air_density, rpm):
        Ps = Cp*air_density*(rpm**3)*(self.propeller.diameter**5)
        return Ps

    def velocity(self, J, rpm):
        v = J*(rpm/60)*(self.propeller.diameter/12)
        return v


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



