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
                                                                'velocity_range': []
                                                },
                                                'plotting': False
                            }
            }
}

