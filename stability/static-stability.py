import math
import sys, os
sys.path.append(os.path.realpath('.'))
from config import config_data
from resources.aircraft_functions import LongitudinalStaticStability as staticLong

# Config data
stability = config_data['analyses']['static_analysis']
parameters = stability['parameters']

# Data extraction section
aircraft = config_data['aircraft']
longitudinal = staticLong(aircraft)

if stability['print_steps'] == True:
    print('Intializing static stability analysis...\n')

# Can be removed if standardized...
longitudinal.setTail('h_stabilizer')
longitudinal.setWing('wing')
longitudinal.setMotor('motor')

# Compute for Cm_0 
epsilon0 = longitudinal.epsilonZero()
epsilonAlpha = longitudinal.epsilonAlpha()

cl0 = longitudinal.liftCoefZero()
alpha0 = longitudinal.alphaZero()
cm0 = longitudinal.cmZero()

# Compute for Cm_alpha
cmAlpha = longitudinal.cmAlpha() 

# Compute Equilibrium Angle - deg
angleEq = longitudinal.alphaEquilibrium()*57.3

if stability['print_steps'] == True:
    print('Analysis finished, check results...')

if stability['print_results'] == True:
    stability_results = {'epsilon0': epsilon0, 'epsilonAlpha': epsilonAlpha, 'cl0': cl0, 'alpha0': alpha0,
                         'cm0': cm0, 'cmAlpha': cmAlpha, 'angleEq': angleEq}
    print(stability_results )

# Missing plot chart
if stability['plotting'] == True:
    pass