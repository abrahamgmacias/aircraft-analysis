import sys, os
sys.path.append(os.path.realpath('.'))

import matplotlib.pyplot as plt
from config import config_data
from resources.aircraft_functions import Aerodynamics as aero

# Config data
trim = config_data['analyses']['trim_analysis']
parameters = trim['parameters']

# Data extraction section
aircraft = config_data['aircraft']
wto = aircraft.weights['mtow']                       # What if weight varies?
cd_0 = aircraft.aeroCoefficients['cd_0']
cl_min_d = aircraft.aeroCoefficients['cl_min_d']

wing = aircraft.getComponents('wing')
sw, arw = wing.geometry['Sw'], wing.geometry['ARw']
ew = wing.wingCoefficients['ew']

motor = aircraft.getComponents('motor')
pa_max = motor.motorSpecs['pa_max']

propeller = motor.getPropeller()
n_eff = propeller.propellerSpecs['n_eff']

rho = parameters['atmospheric_conditions']

# Execution section
results = {'power_required': [], 'power_available': []}         # Missing other values

for v in parameters['velocity_range']:
    coefLift = aero.liftCoefficient(wto, rho, v, sw)
    coefDrag = aero.dragCoefficient(cd_0, coefLift, cl_min_d, ew, arw)

    # Def these functions below...
    # alpha = aero.Closest(CL_list, Angles_list, CL)

    aeroEfficiency = coefLift / coefDrag
    thrustReq = aero.thrustRequired(wto, coefDrag, coefLift)
    powerReq = aero.powerRequired(thrustReq, v)

    powerAvailable = pa_max*n_eff*(rho.atmosphericRatio())
    rClimb = aero.rateOfClimb(powerAvailable, powerReq, wto)

    # Appending section


# Plotting section
if trim['plotting'] == True:
    pass
