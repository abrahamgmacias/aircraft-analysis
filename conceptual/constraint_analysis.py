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

cdMin = aerodynamics['cdMin']
clMax = aerodynamics['clMax']
clTakeOff = aerodynamics['clTakeOff']
cdTakeOff = aerodynamics['cdTakeOff']

# Execution section
results = {'turn': [], 'rateOfClimb': [], 'takeOff': [], 'cruise': []}

for ws in parameters['ws_range']:
    # Computing section
    turn = con.turn(ws, densitySeaLevel, vCruise, cdMin, loadFactor, kFactor)
    rateOfClimb = con.rateOfClimb(ws, densitySeaLevel, vVertical, cdMin, kFactor)
    takeOff = con.takeoff(ws, densitySeaLevel, takeOffDistance, clMax, clTakeOff, cdTakeOff, groundFrictionCoefficient, gravity)
    cruise = con.cruise(ws, densitySeaLevel, vCruise, cdMin, kFactor)

    # Optional constraint section

    # Appending section
    results['rateOfClimb'] += [rateOfClimb]
    results['takeOff'] += [takeOff]
    results['cruise'] += [cruise]
    results['turn'] += [turn]

    


# Exportar la clase a af
# Sacar la plot section y dejarlo en el executable
# Hacer un dict con cada uno de los elementos del init

# Add to aircraft
# conceptualAircraft.addComponents(wing=conceptualWing, motor=conceptualMotor, )

# conceptualAtmosphere = Atmospheric()

# Parametros propuestos
# vs = 12 #Velocidad de stall, m/s
# vc = 14 #Velocidad de crucero, m/s
# vv = 0.508 #Velocidad vertical en Vy

# Aircraft
# mtow = 22 #Peso máximo, kg
# Cdmin = 0.035 #Coeficiente de drag mínimo, tabla 3.1 GUDMUNDSSON

# # Wing
# AR = 7.5 #Relación de aspecto
# lam = 1 #Taper, rectangular
# e = 0.75 #Factor de eficiencia de Oswald                # Make formula for this 
# k = 1/((math.pi)*AR*e)                                  # Make formula for this 

# Motor
# p = 900 #Potencia, Watts
# n = 2 #Factor de carga = 1/cos(angulo de banqueo)
# ep = 0.8 #Eficiencia prop

# tw_real = (ep*p)/(mtow*g*vc)          UNCOMMENT 


# Flight conditions
# rho_sea = 1.225 #densidad del aire a nivel del mar, kg/m^3
# rho_desired = 0.974 #densidad del aire en un lugar en específico, kg/m^3
# mu = 0.05 #Coeficiente de fricción, tabla 17.1 Raymer, Dry concrete/asphalt
# g = 9.807 #Gravedad, m/s^2

# # Aerodynamic conds
# CLmax = 1.2 #Coeficiente de lift máximo deseado
# CLto = 1.2 #Coeficiente de lift deseado en despegue
# CDto = 0.04 #Coeficiente de drag deseado en despegue


# # Add to aircraft
# conceptualAircraft.addComponents(wing=conceptualWing, motor=conceptualMotor)

# # Create objects and populate it
# atmospheric = parameters['atmospheric_conditions']
# aerodynamics = parameters['aerodynamics']
# performance = parameters['performance']
# propulsion = parameters['propulsion']
# velocities = parameters['velocities']
# aircraft = parameters['aircraft']
# wing = parameters['wing']

# conceptualAtmosphere = Atmospheric(atmospheric['densityDesiredLevel'])
# conceptualAtmosphere.setGravity(performance['gravity'])

# if atmospheric['units'] == 'metric':
#     conceptualAtmosphere.setSeaLevel()
# else:
#     conceptualAtmosphere.setSeaLevel(metric=False, imperial=True)

# conceptualAircraft = Aircraft()
# conceptualAircraft.addWeights(mtow=aircraft['mtow'])

# conceptualAircraft.addCoefficients(cdMin=aerodynamics['cdMin'],
#                                    clMax=aerodynamics['clMax'],
#                                    clTakeOff=aerodynamics['clTakeOff'],
#                                    cdTakeOff=aerodynamics['cdTakeOff'])

# conceptualAircraft.addVelocities(vStall=velocities['vStall'],
#                                  vCruise=velocities['vCruise'],
#                                  vVertical=velocities['vVertical'])

# conceptualAircraft.setPerformanceData(loadAtBanking=performance['loadAtBanking'],
#                                       takeOffDistance=performance['takeOffDistance'],
#                                       groundFrictionCoefficient=performance['groundFrictionCoefficient'])

# conceptualWing = Wing()
# conceptualWing.addGeometry(arw=wing['arw'],
#                            lamda=wing['lamda'])
                        
# conceptualWing.addCoefficients(oswaldSpan=wing['oswaldSpan'],
#                                kFactor=wing['kFactor'])

# conceptualPropeller = Propeller()
# conceptualPropeller.addPropellerSpecs(propellerEff=propulsion['propellerEff'])

# conceptualMotor = Motor(conceptualPropeller)
# conceptualMotor.addMotorSpecs(maxPower=propulsion['maxPower'],
#                               thrustToWeight=propulsion['thrustToWeight'])
