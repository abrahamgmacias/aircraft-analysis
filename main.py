import resources.aircraft_functions as af

# General Aircraft Data / general_aircraft_data / aircraft_specs
uas = {'mass': 22, 'ARw': 7.5, 'Sw': 1.831, 'cw': 0.494, 'iw': round(1/57.3, 4), 'St': 0.37, 'lt': 1.48, 
       'it': round(-2.5/57.3, 4), 'n_ef': 0.95, 'Vh': 0.6, 'x_cg_cw': 0.2834, 'x_ac_cw': 0.25,
       'alpha_0w': round(-7.8/57.3, 4)}

# Longitudinal Coefficients / general_coefficients 
long_coeffs = {'aw': 4.40637, 'at': 3.63855, 'ew': 0.825, 'Cm_0w': -0.208, 'CL_0w': 0.62, 'CL_0': 0.6096, 'CD_0': 0.03324}


# --------- Extra Computation ---------- # 
# Class Instantiation
sae = af.Aircraft()
sae.addWeights(mtow=22)
sae.addCoefficients(cl_list=[], cd_list=[])      # Pull from excel file...
sae.addCoefficients(cd_0=0.03324, cl_min_d=0, cl_0=0.6096)

wing = af.Wing(3.71, 0.494, 'MH-114')
wing.addGeometry(ARw=7.5, Sw=1.831)
wing.addCoefficients(ew=0.825)

propeller = af.Propeller(12, 12)   # Change data...
motor = af.Motor(propeller)
sae.addComponents(motor=motor)