from urllib import request
import pandas as pd
import numpy as np
import math

# Remove boilerplate code by using decorators?

# General Aircraft Functions
class Aerodynamics():
    # Aircraft Dynamic Pressure
    @staticmethod
    def dynamicPressure(rho, v):
        return 0.5*rho*v**2

    # Reynolds Number
    @staticmethod
    def numReynolds(rho, v, mac, miu):
        return rho*v*mac / miu

    # Velocity / Stall Velocity / Vuelo Nivelado
    @staticmethod
    def velocity(CL, w, rho, S, f_flaps=1):
        return (2*w / (CL*f_flaps*rho*S))**0.5

    # Velocity / Stall Velocity / Giro Nivelado
    @staticmethod
    def stallVelocity(CL, n_factor, w, rho, S):
        return (2*w*n_factor / (CL*rho*S))**0.5

    # Lift Coefficient / Vuelo Niveglado / Take Off / Touchdown
    @staticmethod
    def liftCoefficient(rho, vel, w, S, f_flaps=1):
        return 2*w / (rho*f_flaps*(vel**2)*S) 

    # Lift Force / Vuelo Nivelado / Take Off / Touchdown
    @staticmethod
    def liftForce(CL, rho, S, v):
        return 0.5*rho*CL*S*(v**2)

    # Drag Coefficient / Vuelo Nivelado
    @staticmethod
    def dragCoefficient(CL, CD_0, e, AR):
        return CD_0 + (CL**2) / (math.pi*e*AR)

    # Drag Curve Slope
    @staticmethod
    def dragCurveSlope(CL, e, ARw, aw):
        return (2*CL / (math.pi*e*ARw))*aw      

    # Thrust Required to maintain Level Flight
    @staticmethod
    def thrustRequired(mass, CL, CD):
        return mass / (CL / CD)

    # Power Required to maintain Level Flight
    @staticmethod
    def powerRequired(T, v):
        return T*v

    @staticmethod
    def rateOfClimb(PA, PR, wto):
        return (PA - PR) / wto


class MetaClass():
    def __init__(self):
        self.aeroCoefficients = {}

    def setData(*args):
        objectToSet, dictToSet, arguments = args
        dictToSet.update(arguments)

    def setAeroCoefficients(self, **kwargs):
        self.setData(self.aeroCoefficients, kwargs)

    def getAeroCoefficients(self, *args, dict=False):
        return self.getData(self.aeroCoefficients, args, dict)
        
    def getData(*args):
        objectToSearch, dictToSearch, arguments, returnDict = args
    
        objectToReturn = []
        if returnDict == True:
            objectToReturn = {}

        if arguments == ():
            if returnDict == False:
                for element in dictToSearch:
                    objectToReturn += [dictToSearch[element]]
                return objectToReturn

            else:
                return dictToSearch

        else:
            for arg in arguments:
                if arg in dictToSearch:
                    if returnDict == False:
                        objectToReturn += [dictToSearch[arg]]
                    else:
                        objectToReturn[arg] = dictToSearch[arg]
                 
            return objectToReturn


# Aircraft class
class Aircraft(MetaClass):
    def __init__(self):
        self.performanceData = {}
        self.components = {}
        self.velocities = {}
        self.weights = {}
        super().__init__()
    
    def setVelocities(self, **kwargs):
        self.setData(self.velocities, kwargs)

    # Components should be classes 
    def setComponents(self, **kwargs):
        self.setData(self.components, kwargs)

    def setPerformanceData(self, **kwargs):
        self.setData(self.performanceData, kwargs)

    def setWeights(self, **kwargs):
        self.setData(self.weights, kwargs) 
        
    def getComponents(self, *args, dict=False): 
        return self.getData(self.components, args, dict)

    def getWeights(self, *args, dict=False): 
        return self.getData(self.weights, args, dict)
    

# Description can be removed if a tail class is created
class Wing(MetaClass):
    def __init__(self, bw=None, cw=None, airfoil=None, description=None):
        self.geometry = {'bw': bw, 'cw': cw, 'airfoil': airfoil}
        self.aeroCoefficients = {}
        self.description = description

    def simpleWingArea(self):
        self.Sw = self.geometry['bw']*self.geometry['cw']

    def setGeometry(self, **kwargs):
         self.setData(self.geometry, kwargs)

    def getGeometry(self, *args, dict=False):
        return self.getData(self.geometry, args, dict)


# English System, correct at the end
class Propeller(MetaClass):
    coefficient_location = {'thrust_coeffs': 'thrust_coefficient_', 'power_coeffs': 'power_coefficient_'}

    def __init__(self, diameter=None, pitch=None):
        self.specs = {'diameter': diameter, 'pitch': pitch}

    def setSpecs(self, **kwargs):
        self.setData(self.specs, kwargs)

    def getSpecs(self, *args, dict=False):
        return self.getData(self.specs, args, dict)
    
    def getPitchDiameterRatio(self):
        pitchDiameterRatio = self.propellerSpecs['pitch'] / self.propellerSpecs['diameter']
        self.propellerSpecs['pitchDiameter'] = pitchDiameterRatio
        return pitchDiameterRatio
        
    # def getData(self, prefix):
    #     return pd.read_csv(f'{prefix}{self.propellerSpecs["pitch_diameter"]}.csv')

    # If input: Cp / Ct - output: J (And viceversa) - Handles Parabolas
    def getCoefficientRatio(self, value, prefix, Cp=False, Ct=False, J=False):
        coefficient_data = self.get_propeller_data(Propeller.coefficient_location[prefix])
        x_points, y_points = coefficient_data['X'], coefficient_data['Y']

        if Cp == True or Ct == True:
            residuals = [abs(value-y) for y in y_points]
            residuals_copy = residuals.copy()

            if Cp == True:
                residuals.sort(reverse=False)
                x_selected_points_position = [residuals_copy.index(y) for y in residuals[0:2]]
                xy_selected_point = [x_points[y_position] for y_position in x_selected_points_position]

            if Ct == True:
                xy_selected_point = x_points[residuals.index(min(residuals))]

        if J == True:
            residuals = [abs(value-x) for x in x_points]
            xy_selected_point = y_points[residuals.index(min(residuals))]
        
        return xy_selected_point

    def propEfficiency(self, Ct, Cp, advance_ratio):
        nProp = Ct*advance_ratio / Cp
        return nProp

    # Check 550 HP - De d??nde sale?
    def cpFactored(self, Ps, air_density, rpm):
        Cp = 550*(Ps / (air_density*((rpm/60)**3)*((self.propellerSpecs['diameter']/12)**5)))
        return Cp

    def ctFactored(self, thrust, air_density, rpm):
        Ct = thrust / (air_density*(rpm**2)*(self.propellerSpecs['diameter']**4))
        return Ct


class Motor(MetaClass):
    def __init__(self, propeller, max_rpm=None, max_voltage=None, max_amperage=None, rpm_voltage=None):
        self.propeller = propeller
        self.specs = {'max_rpm': max_rpm, 'max_voltage': max_voltage,
                      'max_amperage': max_amperage, 'rpm_voltage': rpm_voltage}

    def setData(self, file):
        self.motor_data = pd.read_csv(file)

    def setSpecs(self, **kwargs):
        self.specs.update(kwargs)

    def getPropeller(self):
        return self.propeller 

    def getSpecs(self, *args, dict=False):
        return self.getData(self.specs, args, dict)
    
    def simpleMaxPower(self):
        simple_power = self.specs['max_voltage']*self.specs['max_voltage']
        return round(simple_power, 4)

    def complexMaxPower(self):
        complexPower = self.specs['max_rpm']*(1/self.specs['rpm_voltage'])*self.specs['max_voltage']
        return complexPower

    def thrust(self, Ct, air_density, rpm):
        T = Ct*air_density*((rpm/60)**2)*(self.propeller.diameter/12)**4
        return T

    def motorEfficiency(self, motor_input):
        coefficient_data = pd.read_csv('motor_efficiency.csv')
        x_points, y_points = coefficient_data['X'], coefficient_data['Y']
        residuals = [motor_input-x for x in x_points if motor_input-x >= 0]
        y_selected_point = y_points[residuals.index(min(residuals))]
        return y_selected_point

    def powerS(self, Cp, air_density, rpm):
        Ps = Cp*air_density*(rpm**3)*(self.propeller.diameter**5)
        return Ps

    def velocity(self, J, rpm):
        v = J*(rpm/60)*(self.propeller.diameter/12)
        return v


class Atmospheric(MetaClass):
    def __init__(self, curr_density, curr_temperature=None, curr_height=None):
        self.airConditions = {'currTemperature': curr_temperature, 'currDensity': curr_density,
                              'currHeight': curr_height}

    def setUnits(self, metric=True, imperial=False):
        if metric == True:
            self.units = 'metric'
            return

        else:
            self.units = 'imperial'

    def atmosphericRatio(self):
        if self.units == 'metric':
            return 1.255 / self.airConditions['currDensity']

    # Imperial logic missing
    def setSeaLevel(self, metric=True, imperial=False):
        if metric == True:
            self.seaLevelConditions = {'seaTemperature': None, 'seaDensity': 1.255,
                                       'seaHeight': 0}
            return self.seaLevelConditions

        if imperial == True:
            pass

    def setGravity(self, gravity):
        self.airConditions['gravity'] = gravity
        return 'Gravity set successfully.'

    def getCurrentAirConditions(self, *args, dict=False):
        return self.getData(self.airConditions, args, dict)


class LongitudinalStaticStability():
    def __init__(self, aircraft):
        self.ac = aircraft
        self.acCoefficients = aircraft.aeroCoefficients

    def setTail(self, tailName):
        self.acTail = self.ac.getComponents(tailName)[0]
        self.acTailCoefficients = self.acTail.aeroCoefficients
        self.acWingGeometry = self.acTail.geometry
        return 'Tail object has been set properly.'

    def setWing(self, wingName):
        self.acWing = self.ac.getComponents(wingName)[0]
        self.acWingCoefficients = self.acWing.aeroCoefficients
        self.acWingGeometry = self.acWing.geometry
        return 'Wing object has been set properly.'

    def setMotor(self, motorName):
        self.acMotor = self.ac.getComponents(motorName)[0]
        self.acMotorCoefficients = self.acMotor.specs
        return 'Wing object has been set properly.'
    
    # Cornell
    def cmAlpha(self):
        x_cg_cw, x_ac_cw = self.ac.getAeroCoefficients('x_cg_cw', 'x_ac_cw')
        epsilonAlpha = self.ac.getAeroCoefficients('epsilonAlpha')[0]
        nEff = self.acMotor.propeller.getSpecs('nEff')[0]
        vh, at = self.acTail.getAeroCoefficients('vh', 'at')
        aw = self.acWing.getAeroCoefficients('aw')[0]

        staticMargin = x_cg_cw - x_ac_cw
        cmAlpha  = staticMargin*aw - nEff*vh*at*(1-epsilonAlpha)
        self.ac.setAeroCoefficients(cmAlpha=cmAlpha)
        return round(cmAlpha, 4)

    # Epsilon @ AoA = 0
    def epsilonZero(self): 
        cl_0w = self.acWingCoefficients['cl0w']

        epsilon0 = 2*cl_0w / (self.acWingGeometry['arw']*math.pi) 
        self.ac.setAeroCoefficients(epsilon0=epsilon0)
        return round(epsilon0, 4)

    # Epsilon In Function of AoA - dE/alpha - Downwash - 1/rad
    def epsilonAlpha(self):
        aw = self.acWing.getAeroCoefficients('aw')[0]

        epsilonAlpha = 2*aw / (self.acWingGeometry['arw']*math.pi)
        self.ac.setAeroCoefficients(epsilonAlpha=epsilonAlpha)
        return round(epsilonAlpha, 4)

    def liftCoefZero(self):
        aw, alpha_0w = self.acWing.getAeroCoefficients('aw', 'alpha0w')
        nEff = self.acMotor.propeller.getSpecs('nEff')[0]
        epsilon0 = self.ac.getAeroCoefficients('epsilon0')[0]
        iw, sw = self.acWing.getGeometry('iw', 'sw')
        st, it = self.acTail.getGeometry('st', 'it')
        at = self.acTail.getAeroCoefficients('at')[0]

        cl0 = aw*(iw - alpha_0w) + nEff*(st/sw)*at*(it - epsilon0)
        self.ac.setAeroCoefficients(cl0=cl0)
        return round(cl0, 4)

    def alphaZero(self):
        cl0 = self.ac.getAeroCoefficients('cl0')[0]
        aw = self.acWing.getAeroCoefficients('aw')[0]
        
        alpha0 = -cl0 / aw
        self.ac.setAeroCoefficients(alpha0=alpha0)
        return round(alpha0, 4)

    # Cornell - Eq. 3.17
    def cmZero(self):
        epsilon0, epsilonAlpha = self.ac.getAeroCoefficients('epsilon0', 'epsilonAlpha')
        nEff = self.acMotor.propeller.getSpecs('nEff')[0]
        vh, at = self.acTail.getAeroCoefficients('vh', 'at')
        cm0w = self.acWing.getAeroCoefficients('cm0w')[0]
        alpha0 = self.ac.getAeroCoefficients('alpha0')[0]
        it = self.acTail.getGeometry('it')[0]

        cm0 = cm0w - nEff*vh*at*(it - epsilon0 + (1-epsilonAlpha)*alpha0)
        self.ac.setAeroCoefficients(cm0=cm0)
        return round(cm0, 4)

    # E and R
    def alphaEquilibrium(self):
        cm0, cmAlpha, alpha0 = self.ac.getAeroCoefficients('cm0', 'cmAlpha', 'alpha0')

        alphaEq = (-cm0 / cmAlpha) - abs(alpha0)
        self.ac.setAeroCoefficients(alphaEq=alphaEq)
        return round(alphaEq, 4)


class Constraints():
    # M??todos para encontrar T/W
    # Constant velocity turn
    @staticmethod
    def turn(ws, densitySeaLevel, vCruise, cdMin, loadFactor, kFactor):  
        q_turn = 0.5*(densitySeaLevel)*(vCruise**2)
        firstTerm = (cdMin)/(ws)
        secondTerm = (loadFactor)/(q_turn)
        return q_turn*(firstTerm + kFactor*(secondTerm**2)*ws)
    
    # Rate of Climb
    @staticmethod
    def rateOfClimb(ws, densitySeaLevel, vVertical, cdMin, kFactor): #Rate of Climb
        rateOfClimbSpeed = np.sqrt((2/densitySeaLevel)*ws*np.sqrt(kFactor/(3*cdMin)))
        qClimb = 0.5*densitySeaLevel*(rateOfClimbSpeed**2)
        return (vVertical/rateOfClimbSpeed + (qClimb/ws)*cdMin) + (kFactor/qClimb)*ws

    # Desired Takeoff Distance
    @staticmethod
    def takeoff(ws, densitySeaLevel, takeOffDistance, clMax, clTO, cdTO, groundFrictionCoefficient, gravity): 
        vto = 1.2*np.sqrt(ws*(2/(densitySeaLevel*clMax)))
        qTakeOff = 0.5*densitySeaLevel*(vto**2)
        return (vto**2)/(2*gravity*takeOffDistance) + (qTakeOff*cdTO)/ws + groundFrictionCoefficient*(1-(qTakeOff*clTO)/ws)

    @staticmethod
    def cruise(ws, densitySeaLevel, vCruise, cdMin, kFactor): #Desired cruise Airspeed
        qCruise = 0.5*densitySeaLevel*(vCruise**2)
        return qCruise*cdMin*(1/ws) + kFactor*(1/qCruise)*ws

    # Optional method
    @staticmethod
    def clMax(ws, densitySeaLevel, vStall, vDelta):
        qStall = 0.5*densitySeaLevel*pow(vStall+vDelta, 2)
        return (1/qStall)*ws
        



    # def plot(self):
    #     fig, ax1 = plt.subplots()
    #     ax1.plot(ws, self.turn(), color = 'black',label='Constant velocity turn')
    #     ax1.plot(ws, self.roc(), color = 'red',label='Rate of Climb')
    #     ax1.plot(ws, self.takeoff(), color = 'blue',label='Desired Takeoff Distance')
    #     ax1.plot(ws, self.cruise(), color = 'yellow',label='Desired cruise Airspeed')
    #     ax1.hlines(self.tw_real, 30, 170, colors='k', linestyles='dashed',label='T/W Real posible con Vcrucero')

    #     ax2 = ax1.twinx()
    #     ax2.set_ylabel('CLmax',fontsize=25)
    #     ax2.plot(ws, self.CLmax1(), color = 'k', linestyle='-.',label='Recta Vstall = 12 m/s')
    #     ax2.plot(ws, self.CLmax2(), color = 'b', linestyle='-.',label='Recta Vstall = 10 m/s')
    #     ax2.plot(ws, self.CLmax3(), color = 'r', linestyle='-.',label='Recta Vstall = 14 m/s')
    #     ax2.plot(ws, self.CLmax4(), color = 'magenta', linestyle='-.',label='Recta Vstall = 16 m/s')
    #     #ax2.plot(ws, self.CLmax5(), color = 'green', linestyle='-.',label='Recta Vstall = 8 m/s')
    #     plt.title('Constraint Diagram', fontsize = 25,fontweight='bold')

    #     ax1.set_ylabel('T/W', fontsize = 25)
    #     ax1.set_xlabel(' W/S (N/m^2)', fontsize = 25)
    #     ax1.legend()
    #     ax2.legend(loc='upper center')
    #     fig.tight_layout()
    #     plt.show()
