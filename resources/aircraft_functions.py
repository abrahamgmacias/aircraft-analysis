from urllib import request
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
    def powerRequired(T, v):
        return T*v

    @staticmethod
    def rateOfClimb(PA, PR, wto):
        return (PA - PR) / wto


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

    def getCoefficients(self, *args, dict=False):
        object_to_return = []
        if dict == True:
                object_to_return = {}

        if args == ():
            if dict == False:
                for element in self.aeroCoefficients:
                    object_to_return += [self.aeroCoefficients[element]]
                return object_to_return

            else:
                return self.aeroCoefficients

        else:
            for arg in args:
                if arg in self.aeroCoefficients:
                    if dict == False:
                        object_to_return += [self.aeroCoefficients[arg]]
                    else:
                        object_to_return[arg] = self.aeroCoefficients[arg]
                 
            return object_to_return
        
    def getComponents(self, *args, dict=False): 
        objectToReturn = []
        if dict == True:
            objectToReturn = {}

        for arg in args:
            if arg in self.components:
                element = self.components[arg]
                if dict == False:
                    objectToReturn += [element]
                else:
                    objectToReturn[arg] = element

        return objectToReturn



# Description can be removed if a tail class is created
class Wing():
    def __init__(self, bw=None, cw=None, airfoil=None, description=None):
        self.geometry = {'bw': bw, 'cw': cw, 'airfoil': airfoil}
        self.wingCoefficients = {}
        self.description = description

    def simpleWingArea(self):
        self.Sw = self.geometry['bw']*self.geometry['cw']

    def addGeometry(self, **kwargs):
        self.geometry.update(kwargs)
        return f"'{list(kwargs.keys())[-1]}' parameter was added to wing geometry" 

    def addCoefficients(self, **kwargs):
        self.wingCoefficients.update(kwargs)
        return f"'{list(kwargs.keys())[-1]}' was added to wing coefficients"
        
    def getCoefficients(self, *args, dict=False):
        object_to_return = []
        if dict == True:
                object_to_return = {}

        if args == ():
            if dict == False:
                for element in self.wingCoefficients:
                    object_to_return += [self.wingCoefficients[element]]
                return object_to_return

            else:
                return self.wingCoefficients

        else:
            for arg in args:
                if arg in self.wingCoefficients:
                    if dict == False:
                        object_to_return += [self.wingCoefficients[arg]]
                    else:
                        object_to_return[arg] = self.wingCoefficients[arg]
                 
            return object_to_return


# English System, correct at the end
class Propeller():
    coefficient_location = {'thrust_coeffs': 'thrust_coefficient_', 'power_coeffs': 'power_coefficient_'}

    def __init__(self, diameter=None, pitch=None):
        self.propellerSpecs = {'diameter': diameter, 'pitch': pitch, 'pitch_diameter': pitch/diameter}
        
    def getPropellerData(self, prefix):
        return pd.read_csv(f'{prefix}{self.propellerSpecs["pitch_diameter"]}.csv')

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
        Cp = 550*(Ps / (air_density*((rpm/60)**3)*((self.propellerSpecs['diameter']/12)**5)))
        return Cp

    def ctFactored(self, thrust, air_density, rpm):
        Ct = thrust / (air_density*(rpm**2)*(self.propellerSpecs['diameter']**4))
        return Ct

    def addPropellerSpecs(self, **kwargs):
        self.propellerSpecs.update(kwargs)
        return f"'{list(kwargs.keys())[-1]}' was added to propellerSpecs"


class Motor():
    def __init__(self, propeller, max_rpm=None, max_voltage=None, max_amperage=None, rpm_voltage=None):
        self.propeller = propeller

        self.motorSpecs = {'max_rpm': max_rpm, 'max_voltage': max_voltage,
                           'max_amperage': max_amperage, 'rpm_voltage': rpm_voltage}

    def setMotorData(self, file):
        self.motor_data = pd.read_csv(file)

    def addMotorSpecs(self, **kwargs):
        self.motorSpecs.update(kwargs)
        return f"'{list(kwargs.keys())[-1]}' was added to motorSpecs"

    def simpleMaxPower(self):
        simple_power = self.motorSpecs['max_voltage']*self.motorSpecs['max_voltage']
        return round(simple_power, 4)

    def complexMaxPower(self):
        complexPower = self.motorSpecs['max_rpm']*(1/self.motorSpecs['rpm_voltage'])*self.motorSpecs['max_voltage']
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

    def getPropeller(self):
        return self.propeller    


# Placeholder class to append methods to 
class Atmospheric():
    air_densities = [1.255]

    def __init__(self, curr_density, curr_temperature=None, curr_height=None):
        self.air_conditions = {'curr_temperature': curr_temperature, 'curr_density': curr_density,
                               'curr_height': curr_height}

    def atmosphericRatio(self):
        return Atmospheric.air_densities[0] / self.air_conditions['curr_density']



class LongitudinalStaticStability():
    def __init__(self, aircraft):
        self.ac = aircraft
        self.acCoefficients = aircraft.aeroCoefficients

    def setTail(self, tailName):
        self.acTail = self.ac.getComponents(tailName)[0]
        self.acTailCoefficients = self.acTail.wingCoefficients
        return 'Tail object has been set properly.'

    def setWing(self, wingName):
        self.acWing = self.ac.getComponents(wingName)[0]
        self.acWingCoefficients = self.acWing.wingCoefficients
        return 'Wing object has been set properly.'

    def setMotor(self, motorName):
        self.acMotor = self.ac.getComponents(motorName)[0]
        self.acMotorCoefficients = self.acMotor.wingCoefficients
        return 'Wing object has been set properly.'
    
    # Cornell
    def cmAlpha(self):
        staticMargin = self.ac_coefficients['x_cg_cw'] - self.ac_coefficients['x_ac_cw']
        cm_a  = staticMargin*self.ac_coefficients['aw'] - self.ac['aircraft_specs']['n_ef']*self.ac['aircraft_specs']['Vh']*self.ac['general_coefficients']['at']*(1-self.ac['general_coefficients']['epsilon_a'])
        return cm_a

    # Epsilon @ AoA = 0
    def epsilonZero(self): 
        cl_0w, arw = self.acWingCoefficients.getCoefficients('cl0w', 'arw')
        epsilon_0 = 2*cl_0w / (arw*math.pi) 
        return round(epsilon_0, 4)

    # Epsilon In Function of AoA - dE/alpha - Downwash - 1/rad
    def epsilonAlpha(self):
        aw, arw = self.acWingCoefficients.getCoefficients('aw', 'arw')
        epsilon_alpha = 2*aw / (arw*math.pi)
        return round(epsilon_alpha, 4)

    def liftCoefZero(self):
        CL_0 = self.ac['general_coefficients']['aw']*(self.ac['aircraft_specs']['iw'] - self.ac['aircraft_specs']['alpha_0w']) + self.ac['aircraft_specs']['n_ef']*(self.ac['aircraft_specs']['St']/self.ac['aircraft_specs']['Sw'])*self.ac['general_coefficients']['at']*(self.ac['aircraft_specs']['it'] - self.ac['general_coefficients']['epsilon_0'])
        return round(CL_0, 4)

    def alphaZero(self):
        alpha_0 = -self.ac['general_coefficients']['CL_0'] / self.ac['general_coefficients']['aw']
        return round(alpha_0, 4)

    # Cornell - Eq. 3.17
    def cmZero(self, alpha_0):
        cm_0 = self.ac['general_coefficients']['Cm_0w'] - self.ac['aircraft_specs']['n_ef']*self.ac['aircraft_specs']['Vh']*self.ac['general_coefficients']['at']*(self.ac['aircraft_specs']['it'] - self.ac['general_coefficients']['epsilon_0'] + (1-self.ac['general_coefficients']['epsilon_a'])*alpha_0)
        return round(cm_0, 4)

    # E and R
    def alphaEquilibrium(self):
        alpha_e = (-self.ac['general_coefficients']['Cm_0'] / self.ac['general_coefficients']['Cm_a']) - abs(self.ac['aircraft_specs']['alpha_0'])
        return round(alpha_e, 4)




# class StaticStability():
#     def __init__(self, aircraft):
#         self.ac = aircraft
#         self.ac_coefficients = aircraft.aeroCoefficients

#     # Cornell
#     def cmAlpha(self):
#         cm_a  = (self.ac['aircraft_specs']['x_cg_cw'] - self.ac['aircraft_specs']['x_ac_cw'])*self.ac['general_coefficients']['aw'] - self.ac['aircraft_specs']['n_ef']*self.ac['aircraft_specs']['Vh']*self.ac['general_coefficients']['at']*(1-self.ac['general_coefficients']['epsilon_a'])
#         return cm_a

#     # Epsilon @ AoA = 0
#     def epsilonZero(self): 
#         epsilon_0 = 2*self.ac['general_coefficients']['CL_0w'] / (self.ac['aircraft_specs']['ARw']*math.pi) 
#         return round(epsilon_0, 4)

#     # Epsilon In Function of AoA - dE/alpha - Downwash - 1/rad
#     def epsilonAlpha(self):
#         epsilon_alpha = 2*self.ac['general_coefficients']['aw'] / (self.ac['aircraft_specs']['ARw']*math.pi)
#         return round(epsilon_alpha, 4)

#     def liftCoefZero(self):
#         CL_0 = self.ac['general_coefficients']['aw']*(self.ac['aircraft_specs']['iw'] - self.ac['aircraft_specs']['alpha_0w']) + self.ac['aircraft_specs']['n_ef']*(self.ac['aircraft_specs']['St']/self.ac['aircraft_specs']['Sw'])*self.ac['general_coefficients']['at']*(self.ac['aircraft_specs']['it'] - self.ac['general_coefficients']['epsilon_0'])
#         return round(CL_0, 4)

#     def alphaZero(self):
#         alpha_0 = -self.ac['general_coefficients']['CL_0'] / self.ac['general_coefficients']['aw']
#         return round(alpha_0, 4)

#     # Cornell - Eq. 3.17
#     def cmZero(self, alpha_0):
#         cm_0 = self.ac['general_coefficients']['Cm_0w'] - self.ac['aircraft_specs']['n_ef']*self.ac['aircraft_specs']['Vh']*self.ac['general_coefficients']['at']*(self.ac['aircraft_specs']['it'] - self.ac['general_coefficients']['epsilon_0'] + (1-self.ac['general_coefficients']['epsilon_a'])*alpha_0)
#         return round(cm_0, 4)

#     # E and R
#     def alphaEquilibrium(self):
#         alpha_e = (-self.ac['general_coefficients']['Cm_0'] / self.ac['general_coefficients']['Cm_a']) - abs(self.ac['aircraft_specs']['alpha_0'])
#         return round(alpha_e, 4)