import pandas as pd
import math

# MANY METHODS SHOULD BE ADDED THE @staticmethod DECORATOR

# Aircraft Class
class Aircraft():
    def __init__(self, propulsion_system=None, general_aircraft_data=None, general_coefficients=None, longitudinal_dynamic_data=None):
        self.propulsion = propulsion_system
        self.aircraft_specs = general_aircraft_data
        self.general_coefficients = general_coefficients
        self.longitudinal_dynamic_data = longitudinal_dynamic_data

        # self.flight_envelope = {}
        self.temp_data = {}

    def addPropulsionSystem(self, propulsion_system):
        self.propulsion = propulsion_system

    def addSpeed(self, speed):
        self.flight_envelope['speed'] = speed

    def addAlpha(self, alpha):
        self.flight_envelope['alpha'] = math.radians(alpha)

    # def addAirDensity(self, air_density):
    #     self.flight_envelope['air_density'] = air_density

    def addSpec(self, spec, name):
        self.aircraft_specs[name] = spec 

    def addCoefficient(self, coeff, name):
        if type(coeff) == list:
            for coefficient, title in zip(coeff, name):
                self.general_coefficients[title] = coefficient
        else:
            self.general_coefficients[name] = coeff

    def mergeToCoefficients(self, new_dict):
        self.general_coefficients.update(new_dict)

    def addTemp(self, temp, name):
        if type(temp) == list:
            for coefficient, title in zip(temp, name):
                self.temp_data[title] = coefficient
        else:
            self.temp_data[name] = temp


# General Aircraft Functions
class Aerodynamics():
    def __init__(self, aircraft=None):
        self.ac = aircraft.__dict__
        self.g = 9.81

    def changeUnits(self):
        self.g = int(input('Gravity: '))

    def addAircraft(self, aircraft):
        self.ac = aircraft

    def aircraftCheck(self):
        if self.ac != None:
            return True
        else:
            return False

    def showData(self):
        while self.aircraftCheck():
            print(self.ac)
    
    # Aircraft Dynamic Pressure
    def dynamicPressure(self, rho=None, v=None):
        while self.aircraftCheck():
            return 0.5*self.ac['flight_envelope']['air_density']*self.ac['flight_envelope']['speed']**2

        return 0.5*rho*v**2

    # Aircraft Dynamic Pressure
    def dynamicPressure(self, rho=None, v=None):
        while self.aircraftCheck():
            return 0.5*self.ac['flight_envelope']['air_density']*self.ac['flight_envelope']['speed']**2

        return 0.5*rho*v**2

    # Reynolds Number
    def numReynolds(self, rho, v, MAC, miu):
        return rho*v*MAC / miu

    # Velocity / Stall Velocity / Vuelo Nivelado
    def velocity(self, CL, f_flaps, w=None, rho=None, S=None):
        while self.aircraftCheck():
            return (2*self.ac['aircraft_specs']['mass']*9.81 / (CL*f_flaps*self.ac['flight_envelope']['air_density']*self.ac['aircraft_specs']['Sw']))**0.5
        
        return (2*w / (CL*f_flaps*rho*S))**0.5

    # Velocity / Stall Velocity / Giro Nivelado
    def stallVelocity(self, CL, n_factor, w=None, rho=None, S=None):
        while self.aircraftCheck():
            return (2*self.ac['aircraft_specs']['mass']*9.81*n_factor / (CL*self.ac['flight_envelope']['air_density']*self.ac['aircraft_specs']['Sw']))**0.5

        return (2*w*n_factor / (CL*rho*S))**0.5

    # Lift Coefficient / Vuelo Niveglado / Take Off / Touchdown
    def liftCoefficient(self, f_flaps, rho, vel, w=None, S=None):
        while self.aircraftCheck():
            return 2*self.ac['aircraft_specs']['mass']*self.g / (rho*f_flaps*(vel**2)*self.ac['aircraft_specs']['Sw'])

        return 2*w / (rho*f_flaps*(vel**2)*S) 

    # Lift Force / Vuelo Nivelado / Take Off / Touchdown
    def liftForce(self, CL, rho=None, S=None, v=None):
        while self.aircraftCheck():
            return 0.5*self.ac['flight_envelope']['air_density']*CL*self.ac['aircraft_specs']['Sw']*(self.ac['flight_envelope']['speed']**2)

        return 0.5*rho*CL*S*(v**2)

    # Drag Coefficient / Vuelo Nivelado
    def dragCoefficient(self, CL, CD_0=None, e=None, AR=None):
        while self.aircraftCheck():
            return self.ac['general_coefficients']['CD_0'] + ((CL**2) / (math.pi*self.ac['general_coefficients']['ew']*self.ac['aircraft_specs']['ARw']))

        return CD_0 + (CL**2) / (math.pi*e*AR)

    # Drag Curve Slope
    def dragCurveSlope(self, CL, e=None, ARw=None, aw=None):
        while self.aircraftCheck():
            return self.ac['general_coefficients']['ew']
            #return (2*CL / (math.pi*self.ac['general_coefficients']['ew']*self.ac['aircraft_specs']['ARw']))*self.ac['general_coefficients']['aw']

        return (2*CL / (math.pi*e*ARw))*aw      

    # Thrust Required to maintain Level Flight
    def thrustRequired(self, CL, CD):
        T = self.ac['aircraft_specs']['mass'] / (CL / CD)
        return round(T, 4)

    # Aircraft Speed given RPM and Thrust
    def aircraftVelocity(self, rpm, thrust, air_density):
        try:
            Ct = self.ac['propulsion'].propeller.ct_factored(thrust, rpm, air_density)
        except AttributeError:
            print('Aircraft object has no Propulsion System.')
            return
        else:
            J_coeff = self.ac['propulsion'].propeller.get_coefficient_ratio(Ct, prefix='thrust_coeffs', Ct=True)
            vel = self.ac['propulsion'].velocity(J_coeff, rpm)   
            return vel


# rpm_voltage - rpm to voltage ratio
class Motor():
    def __init__(self, propeller, max_rpm=None, max_voltage=None, max_amperage=None, rpm_voltage=None):
        self.propeller = propeller
        self.max_rpm = max_rpm
        self.max_voltage = max_voltage
        self.max_amperage = max_amperage
        self.rpm_voltage = rpm_voltage

        self.motor_data = self.set_motor_data('motor_data.csv')

    def set_motor_data(self, file):
        self.motor_data = pd.read_csv(file)

    def get_motor_data(self):
        return self.motor_data

    def simple_max_power(self):
        simple_power = self.max_voltage*self.max_amperage
        return round(simple_power, 4)

    def thrust(self, Ct, air_density, rpm):
        T = Ct*air_density*((rpm/60)**2)*(self.propeller.diameter/12)**4
        return round(T, 4)

    def motor_efficiency(self, motor_input):
        coefficient_data = pd.read_csv('motor_efficiency.csv')
        x_points, y_points = coefficient_data['X'], coefficient_data['Y']
        residuals = [motor_input-x for x in x_points if motor_input-x >= 0]
        y_selected_point = y_points[residuals.index(min(residuals))]
        return y_selected_point

    def power_s(self, Cp, air_density, rpm):
        Ps = Cp*air_density*(rpm**3)*(self.propeller.diameter**5)
        return round(Ps, 4)

    def velocity(self, J, rpm):
        v = J*(rpm/60)*(self.propeller.diameter/12)
        return round(v, 4)


# English System, correct at the end
class Propeller():
    coefficient_location = {'thrust_coeffs': 'thrust_coefficient_', 'power_coeffs': 'power_coefficient_'}

    def __init__(self, diameter, pitch):
        self.diameter = diameter
        self.pitch = pitch
        self.pitch_diameter = pitch/diameter

    def get_propeller_data(self, prefix):
        return pd.read_csv(f'{prefix}{self.pitch_diameter}.csv')

    # If input: Cp / Ct - output: J (And viceversa) - Handles Parabolas
    def get_coefficient_ratio(self, value, prefix, Cp=False, Ct=False, J=False):
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