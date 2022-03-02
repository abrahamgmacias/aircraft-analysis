import math
import sys, os
sys.path.append(os.path.realpath('.'))
from config import config_data
from resources.aircraft_functions import LongitudinalStaticStability as staticLong

# Config data
trim = config_data['analyses']['static_analysis']
parameters = trim['parameters']

# Data extraction section
aircraft = config_data['aircraft']
longitudinal = staticLong(aircraft)

# Can be removed if standardized...
longitudinal.setTail('h_stabilizer')
longitudinal.setWing('wing')
longitudinal.setMotor('motor')

# # Compute for Cm_0 
epsilon0 = longitudinal.epsilonZero()
epsilonA = longitudinal.epsilonAlpha()

cl0 = longitudinal.liftCoefZero()

# alpha0 = longitudinal.alphaZero()
# aircraft.addCoefficients(alpha0=alpha0)

# cm0 = longitudinal.cmZero(alpha0)
# aircraft.addCoefficients(cm0=cm0)

# # Compute for Cm_alpha
# cmA = longitudinal.cmAlpha() 

# # Compute Equilibrium Angle - deg
# aE = longitudinal.alphaEquilibrium()*57.3
# aircraft.addCoefficients(cmA=cmA, aE=aE)