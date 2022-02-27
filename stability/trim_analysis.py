import sys, os

from matplotlib.style import available
sys.path.append(os.path.realpath('.'))

import matplotlib.pyplot as plt
from config import config_data
from resources.aircraft_functions import Aerodynamics as aero

# Config data
trim = config_data['analyses']['trim_analysis']
parameters = trim['parameters']
velocityRange = parameters['velocity_range']

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
trim_results = {}                                           # Missing other values
required_list, available_list, errorList = [], [], []

if trim['print_steps'] == True:
    print('Initializing trim analysis...\n')

for v in velocityRange:
    if trim['print_steps'] == True:
        print(f'Computing trim @ v = {v} m/s...')

    coefLift = aero.liftCoefficient(wto, rho, v, sw)
    coefDrag = aero.dragCoefficient(coefLift, cd0, ew, arw)
    aeroEfficiency = coefLift / coefDrag

    try:
        thrustReq = aero.thrustRequired(wto, coefDrag, coefLift)

    except ZeroDivisionError:
        errorList += [v]

        if trim['print_steps'] == True:
            print(f'v = {v} m/s produces computing error...')

    else:
        powerReq = aero.powerRequired(thrustReq, v)
        powerAvail = paMax*nEff*atmosphere.atmosphericRatio()
        rClimb = aero.rateOfClimb(powerAvail, powerReq, wto)

        # Appending section
        trim_results[v] = (coefLift, coefDrag, thrustReq, powerReq, powerAvail, rClimb)
        required_list += [powerReq]
        available_list += [powerAvail]

for error in errorList:
    index = velocityRange.index(error)
    velocityRange.pop(index)

if (len(trim_results.keys()) > 0) and (trim['print_steps'] == True):
    print('\nAnalysis finished, check results...')


# Print results section
if trim['print_results']:
    print('')
    print(trim_results)


# Plotting section
if trim['plotting'] == True:
    # Power required / power available as a function of trim velocity
    plt.title('PR, PA vs. Trim Velocity')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('PR, PA (Watts')

    plt.plot(parameters['velocity_range'], required_list)
    plt.plot(parameters['velocity_range'], available_list)
    plt.show()
