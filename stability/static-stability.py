import sys
import math 
# Cambiar para CWD
sys.path.append("D:\Scripts\Python Scripts\SAE-New-Structure\Lib")   
import aircraft_functions as af 
import aircraft_data as ad


# Class Instantiation
sae = ad.sae
sae_longitudinal = Static_Stability(sae)

# Compute for Cm_0 
epsilon_0 = sae_longitudinal.epsilonZero()
epsilon_a = sae_longitudinal.epsilonAlpha()
sae.addCoefficient(epsilon_0, 'epsilon_0')
sae.addCoefficient(epsilon_a, 'epsilon_a')

epsilon_0 = sae_longitudinal.epsilonZero()
epsilon_a = sae_longitudinal.epsilonAlpha()
sae.addCoefficient(epsilon_0, 'epsilon_0')
sae.addCoefficient(epsilon_a, 'epsilon_a')

CL_0 = sae_longitudinal.liftCoefZero()
sae.addCoefficient(CL_0, 'CL_0')

alpha_0 = sae_longitudinal.alphaZero()
sae.addSpec(alpha_0, 'alpha_0')
cm_0 = sae_longitudinal.cmZero(alpha_0)
sae.addCoefficient(cm_0, 'Cm_0')

# Compute for Cm_alpha
cm_a = sae_longitudinal.cmAlpha() 
sae.addCoefficient(cm_a, 'Cm_a')

# Compute Equilibrium Angle - deg
a_e = sae_longitudinal.alphaEquilibrium()*57.3