import math
import sys, os
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.path.realpath('.'))
from config import config_data
from resources.aircraft_functions import Constraints as con

# Import section
constraints = config_data['analyses']['constraint_analysis']
parameters = constraints['parameters']
optional = parameters['optional']

# Import parameters
atmospheric = parameters['atmospheric_conditions']
aerodynamics = parameters['aerodynamics']
performance = parameters['performance']
propulsion = parameters['propulsion']
aircraft = parameters['aircraft']
wing = parameters['wing']

densitySeaLevel = atmospheric['densitySeaLevel']

arw = wing['arw']
kFactor = wing['kFactor']
oswaldSpan = wing['oswaldSpan']

gravity = performance['gravity']
loadFactor = performance['loadAtBanking']
takeOffDistance = performance['takeOffDistance']
groundFrictionCoefficient = performance['groundFrictionCoefficient']
vVertical = performance['vVertical']
vCruise = performance['vCruise']
vStall = performance['vStall']

cdMin = aerodynamics['cdMin']
clMax = aerodynamics['clMax']
clTakeOff = aerodynamics['clTakeOff']
cdTakeOff = aerodynamics['cdTakeOff']

# Execution section
if parameters['print_steps']:
    print('Intializing constraint analysis...\n')

results = {'turn': [], 'rateOfClimb': [], 'takeOff': [], 'cruise': []}
if optional['status']:
    results['deltaRangeResults'] = []

for ws in parameters['ws_range']:
    if parameters['print_steps']:
        print(f'Computing T/W for W/S = {ws}')

    # Computing section
    try: 
        turn = con.turn(ws, densitySeaLevel, vCruise, cdMin, loadFactor, kFactor)
        rateOfClimb = con.rateOfClimb(ws, densitySeaLevel, vVertical, cdMin, kFactor)
        takeOff = con.takeoff(ws, densitySeaLevel, takeOffDistance, clMax, clTakeOff, cdTakeOff, groundFrictionCoefficient, gravity)
        cruise = con.cruise(ws, densitySeaLevel, vCruise, cdMin, kFactor)
    except:
        print(f'W/S = {ws} produces computing error...')

    # Optional constraint section
    if optional['status'] == True:
        deltaIndividual = []

        for delta in optional['deltaRange']:
            clMax = con.clMax(ws, densitySeaLevel, vStall, delta)
            deltaIndividual += [clMax]

        results['deltaRangeResults'] += [deltaIndividual]

    # Appending section
    results['rateOfClimb'] += [rateOfClimb]
    results['takeOff'] += [takeOff]
    results['cruise'] += [cruise]
    results['turn'] += [turn]

if parameters['print_results']:
    print('\nAnalysis finished, check results... \n')

if parameters['print_results']:
    print(results)

if parameters['plotting']:
    pass
    

    
