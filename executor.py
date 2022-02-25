from main import general_aircraft


# Executor soon-to-be JSON portion
trim_analysis = True
input_parameters = None

config = {
            'aircraft': general_aircraft,
            'analyses': {
                            'trim_analysis': {
                                                'execution': False,
                                                'air_density': False,
                                                'plotting': False
                            }
            }
}



# Algo
if config['analyses']['trim_analysis']['execution'] == True:
    if config['analyses']['trim_analysis']['plotting'] == True:
        pass

    # else:
    #     from stability.trim_analysis import 
