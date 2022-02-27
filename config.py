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
                                                 'print_results': False,
                                                 'print_steps': False
                            }
            }
}

