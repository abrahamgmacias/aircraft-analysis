from re import M
from main import *


# Executor soon-to-be JSON portion
trim_analysis = True
input_parameters = None

config_data = {
            'aircraft': general_aircraft,
            'analyses': {
                            'trim_analysis': {
                                                'parameters': {
                                                                'atmospheric_conditions': atmosphere,
                                                                'velocity_range': list(range(0, 8))
                                                },
                                                'plotting': False,
                                                'print_results': False,
                                                'print_steps': False
                            },

                            'static_analysis': {
                                                 'parameters': {},
                                                 'plotting': False,
                                                 'print_results': False,
                                                 'print_steps': False
                            },

                            'constraint_analysis': {
                                                 'parameters': {
                                                                 'velocities': {
                                                                                'vStall': None,
                                                                                'vCruise': None,
                                                                                'vVertical': None,
                                                                 }, 

                                                                 'atmospheric_conditions': {
                                                                           'units': 'metric',
                                                                           'densityDesiredLevel': None
                                                                 }, 

                                                                 'wing': {
                                                                           'arw': None,
                                                                           'lamda': None,
                                                                           'oswaldSpan': None,
                                                                           'kFactor': None
                                                                 },

                                                                 'propulsion': {
                                                                           'maxPower': None,
                                                                           'propellerEff': None,
                                                                           'thrustToWeight': None
                                                                 },

                                                                 'aircraft': {
                                                                           'mtow': None
                                                                 },

                                                                 'performance': {
                                                                           'loadAtBanking': None,
                                                                           'takeOffDistance': None,
                                                                           'groundFrictionCoefficient': None,
                                                                           'gravity': 9.807
                                                                 }, 

                                                                 'aerodynamics': {
                                                                           'cdMin': None,
                                                                           'clMax': None,
                                                                           'clTakeOff': None,
                                                                           'cdTakeOff': None,
                                                                 }
                                                 },

                                                 'plotting': False,
                                                 'print_results': False,
                                                 'print_steps': False
                            }
            }
}

