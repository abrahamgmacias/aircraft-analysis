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
cd0 = aircraft.aeroCoefficients['cd_0']
clMinD = aircraft.aeroCoefficients['cl_min_d']

wing = aircraft.getComponents('wing')
sw, arw = wing.geometry['Sw'], wing.geometry['ARw']
ew = wing.wingCoefficients['ew']

motor = aircraft.getComponents('motor')
paMax = motor.motorSpecs['pa_max']

propeller = motor.getPropeller()
nEff = propeller.propellerSpecs['n_eff']

atmosphere = parameters['atmospheric_conditions']
rho = atmosphere.air_conditions['curr_density']

# Execution section
trim_results = {}                               # Missing other values

for v in parameters['velocity_range']:
    coefLift = aero.liftCoefficient(wto, rho, v, sw)
    coefDrag = aero.dragCoefficient(coefLift, cd0, ew, arw)

    aeroEfficiency = coefLift / coefDrag
    try:
        thrustReq = aero.thrustRequired(wto, coefDrag, coefLift)
    except ZeroDivisionError:
        pass
    else:
        powerReq = aero.powerRequired(thrustReq, v)
        powerAvail = paMax*nEff*atmosphere.atmosphericRatio()
        rClimb = aero.rateOfClimb(powerAvail, powerReq, wto)

        # Appending section
        trim_results[v] = (coefLift, coefDrag, thrustReq, powerReq, powerAvail, rClimb)

if trim['print_results']:
    print(trim_results)

# Plotting section
if trim['plotting'] == True:
    pass



# Handle division by zeroes given null velocity