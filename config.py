from re import M
from main import *
import numpy as np

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
                                                                 'atmospheric_conditions': {
                                                                           'densitySeaLevel': 1.225,
                                                                           'densityDesiredLevel': None
                                                                 }, 

                                                                 'wing': {
                                                                           'arw': None,
                                                                           'lamda': None,
                                                                           'kFactor': None,
                                                                           'oswaldSpan': None
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
                                                                           'vStall': None,
                                                                           'vCruise': None,
                                                                           'vVertical': None,
                                                                           'gravity': 9.807,
                                                                           'loadAtBanking': None,
                                                                           'takeOffDistance': None,
                                                                           'groundFrictionCoefficient': None
                                                                 }, 

                                                                 'aerodynamics': {
                                                                           'cdMin': None,
                                                                           'clMax': None,
                                                                           'clTakeOff': None,
                                                                           'cdTakeOff': None,
                                                                 },

                                                                 'ws_range': np.linspace(30, 170, 171),
                                                 },

                                                 'plotting': False,
                                                 'print_results': False,
                                                 'print_steps': False
                            }
            }
}