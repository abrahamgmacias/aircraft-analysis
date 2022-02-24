import sys, os
sys.path.append(os.path.realpath('.'))

from main import sae, wing, motor
from resources.aircraft_functions import Aerodynamics as aero


# Extract data from wing
Sw, ARw = wing.geometry['Sw'], wing.geometry['ARw']
ew = wing.wingCoefficients['ew']

# Extract data from aircraft
wto = sae.weights['mtow']                       # What if weight varies?
cd_0 = wing.wingCoefficients['cd_0']
cl_min_d = wing.wingCoefficients['cl_min_d']

# Extract data from motor
# motor

# Non class variables
rho = 1.225

rangoDeVelocidades = []
for v in rangoDeVelocidades:
    coefLift = aero.liftCoefficient(wto, rho, v, Sw)
    coefDrag = aero.dragCoefficient(cd_0, coefLift, cl_min_d, ew, ARw)

    # Def these functions below...
    # alpha = aero.Closest(CL_list, Angles_list, CL)

    aeroEfficiency = coefLift / coefDrag
    thrustReq = aero.thrustRequired(wto, coefDrag, coefLift)
    powerReq = aero.powerRequired(thrustReq, v)

    pAvailable = PA_max*n_ef*(rho / rho_list[0])
    rClimb = af.R_C(PA, PR, wto)

    # Appending section

# Plotting section
