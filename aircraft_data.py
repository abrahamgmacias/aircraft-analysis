import sys
# Cambiar dependiendo de la ubicacion de sus archivos
sys.path.append("D:\Scripts\Python Scripts\SAE-New-Structure\Lib")
import aircraft_functions as af


# General Aircraft Data / general_aircraft_data / aircraft_specs
uas = {'mass': 22, 'ARw': 7.5, 'Sw': 1.831, 'cw': 0.494, 'iw': round(1/57.3, 4), 'St': 0.37, 'lt': 1.48, 
       'it': round(-2.5/57.3, 4), 'n_ef': 0.95, 'Vh': 0.6, 'x_cg_cw': 0.2834, 'x_ac_cw': 0.25,
       'alpha_0w': round(-7.8/57.3, 4)}

# Longitudinal Coefficients / general_coefficients 
long_coeffs = {'aw': 4.40637, 'at': 3.63855, 'ew': 0.825, 'Cm_0w': -0.208, 'CL_0w': 0.62, 'CL_0': 0.6096, 'CD_0': 0.03324}

# Lateral Coefficients / general_coefficients
lat_coeffs = {}


# --------- Extra Computation ---------- # 
# Class Instantiation
sae = af.Aircraft(uas, long_coeffs)
aerodynamics = af.Aerodynamics(sae.__dict__)

#CD_alpha = aerodynamics.dragCurveSlope(0.16)                   # Revisar ; Datos incorrectos - Test
#sae.addCoefficient(CD_alpha, 'CD_alpha')