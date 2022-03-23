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
wsRange = parameters['ws_range']
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
    for delta in optional['deltaRange']:
        results[f'deltaRange{delta}'] = []

for ws in wsRange:
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
        for delta in optional['deltaRange']:
            clMax = con.clMax(ws, densitySeaLevel, vStall, delta)
            results[f'deltaRange{delta}'] += [clMax]

    # Appending section
    results['rateOfClimb'] += [rateOfClimb]
    results['takeOff'] += [takeOff]
    results['cruise'] += [cruise]
    results['turn'] += [turn]

if parameters['print_steps']:
    print('\nAnalysis finished, check results... \n')

if parameters['print_results']:
    print(results)

# Plotting section
if parameters['plotting']:
    plt.plot(wsRange, results['turn'])
    plt.plot(wsRange, results['rateOfClimb'])
    plt.plot(wsRange, results['takeOff'])
    plt.plot(wsRange, results['cruise'])

    if optional['status']:
        for delta in optional['deltaRange']:
            plt.plot(wsRange, results[f'deltaRange{delta}'])

    plt.show()


    
# def plot(self):
#         fig, ax1 = plt.subplots()
#         ax1.plot(self.ws, self.turn(), color = 'black',label='Constant velocity turn')
#         ax1.plot(self.ws, self.roc(), color = 'red',label='Rate of Climb')
#         ax1.plot(self.ws, self.takeoff(), color = 'blue',label='Desired Takeoff Distance')
#         ax1.plot(self.ws, self.cruise(), color = 'yellow',label='Desired cruise Airspeed')
#         ax1.hlines(self.tw_real, 30, 170, colors='k', linestyles='dashed',label='T/W Real posible con Vcrucero')
#         ax2 = ax1.twinx()
#         ax2.set_ylabel('CLmax',fontsize=25)
#         ax2.plot(self.ws, self.CLmax1(), color = 'k', linestyle='-.',label='Recta Vstall = 12 m/s')
#         ax2.plot(self.ws, self.CLmax2(), color = 'b', linestyle='-.',label='Recta Vstall = 10 m/s')
#         ax2.plot(self.ws, self.CLmax3(), color = 'r', linestyle='-.',label='Recta Vstall = 14 m/s')
#         ax2.plot(self.ws, self.CLmax4(), color = 'magenta', linestyle='-.',label='Recta Vstall = 16 m/s')
#         #ax2.plot(self.ws, self.CLmax5(), color = 'green', linestyle='-.',label='Recta Vstall = 8 m/s')
#         plt.title('Constraint Diagram', fontsize = 25,fontweight='bold')
#         ax1.set_ylabel('T/W', fontsize = 25)
#         ax1.set_xlabel(' W/S (N/m^2)', fontsize = 25)
#         ax1.legend()
#         ax2.legend(loc='upper center')
#         fig.tight_layout()
#         plt.show()
    
