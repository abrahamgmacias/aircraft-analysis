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
                                                                           'densityDesiredLevel': 0.974
                                                                 }, 

                                                                 'wing': {
                                                                           'arw': 7.5,
                                                                           'lamda': 1,
                                                                           'kFactor': 0.056,
                                                                           'oswaldSpan': 0.75
                                                                 },

                                                                 'propulsion': {
                                                                           'maxPower': 900,
                                                                           'propellerEff': 0.8,
                                                                           'thrustToWeight': 0.24
                                                                 },

                                                                 'aircraft': {
                                                                           'mtow': 22
                                                                 },

                                                                 'performance': {
                                                                           'vStall': 12,
                                                                           'vCruise': 14,
                                                                           'vVertical': 0.508,
                                                                           'gravity': 9.807,
                                                                           'loadAtBanking': 2,
                                                                           'takeOffDistance': 61,
                                                                           'groundFrictionCoefficient': 0.05
                                                                 }, 

                                                                 'aerodynamics': {
                                                                           'cdMin': 0.035,
                                                                           'clMax': 1.2,
                                                                           'clTakeOff': 1.2,
                                                                           'cdTakeOff': 0.04
                                                                 },

                                                                 'ws_range': range(30, 170),
                                                                 'optional': {
                                                                            'status': True,
                                                                            'deltaRange': [(0, 'green', 'Recta Vstall = 8 m/s'),
                                                                                           (-2, 'k', 'Recta Vstall = 12 m/s'),
                                                                                           (2, 'b', 'Recta Vstall = 10 m/s'),
                                                                                           (-4, 'r', 'Recta Vstall = 14 m/s'),
                                                                                           (4, 'magenta', 'Recta Vstall = 16 m/s')]
                                                 },

                                                 'plotting': True,
                                                 'print_results': False,
                                                 'print_steps': False
                            }
                        }
            }
}