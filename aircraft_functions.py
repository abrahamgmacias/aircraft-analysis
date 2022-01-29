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

    




