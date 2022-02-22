import sys, os
sys.path.append(os.path.realpath('.'))

from main import sae
from resources.aircraft_functions import Aerodynamics as aero

rangoDeVelocidades = []

# Extract data from sae class...
print(sae.getComponents('motor'))

for v in rangoDeVelocidades:
    coefLift = aero.liftCoefficient(wto, rho, v, Sw)
    coefDrag = aero.dragCoefficient(CD0, coefLift, CL_min_D, ew, ARw)

    # Def these functions below...
    alpha = aero.Closest(CL_list, Angles_list, CL)

    aeroEfficiency = coefLift / coefDrag
    thrustReq = aero.thrustRequired(wto, coefDrag, coefLift)
    powerReq = aero.powerRequired(thrustReq, v)

    pAvailable = PA_max*n_ef*(rho / rho_list[0])
    rClimb = af.R_C(PA, PR, wto)

    # Appending section

# Plotting section
