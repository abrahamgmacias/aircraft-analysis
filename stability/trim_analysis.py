import sys, os
sys.path.append(os.path.realpath('.'))

import matplotlib.pyplot as plt
from config import config_data
from main import general_aircraft, wing, motor, propeller
from resources.aircraft_functions import Aerodynamics as aero

# Extract data from wing
Sw, ARw = wing.geometry['Sw'], wing.geometry['ARw']
ew = wing.wingCoefficients['ew']

# Extract data from aircraft
wto = general_aircraft.weights['mtow']                       # What if weight varies?
cd_0 = general_aircraft.aeroCoefficients['cd_0']
cl_min_d = general_aircraft.aeroCoefficients['cl_min_d']

# Extract data from motor and propeller
pa_max = motor.motorSpecs['pa_max']
n_eff = propeller.propellerSpecs['n_eff']

results = {'power_required': [], 'power_available': []}         # Missing other values

for v in config_data['analyses']['trim_analysis']['parameters']['velocity_range']:
    coefLift = aero.liftCoefficient(wto, rho, v, Sw)
    coefDrag = aero.dragCoefficient(cd_0, coefLift, cl_min_d, ew, ARw)

    # Def these functions below...
    # alpha = aero.Closest(CL_list, Angles_list, CL)

    aeroEfficiency = coefLift / coefDrag
    thrustReq = aero.thrustRequired(wto, coefDrag, coefLift)
    powerReq = aero.powerRequired(thrustReq, v)

    powerAvailable = pa_max*n_eff*(rho / rho_list[0])
    rClimb = aero.rateOfClimb(powerAvailable, powerReq, wto)

    # Appending section


# Plotting section
if config_data['analyses']['trim_analysis']['plotting'] == True:
    
