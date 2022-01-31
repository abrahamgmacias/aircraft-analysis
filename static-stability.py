import sys
import math 
# Cambiar para CWD
sys.path.append("D:\Scripts\Python Scripts\SAE-New-Structure\Lib")   
import aircraft_functions as af 
import aircraft_data as ad


# ------- Classes ------- # 
class Static_Stability(af.Aircraft):
    def __init__(self, aircraft_obj):
        super().__init__(ad.uas, ad.long_coeffs)
        self.ac = aircraft_obj.__dict__

    # Cornell
    def cmAlpha(self):
        cm_a  = (self.ac['aircraft_specs']['x_cg_cw'] - self.ac['aircraft_specs']['x_ac_cw'])*self.ac['general_coefficients']['aw'] - self.ac['aircraft_specs']['n_ef']*self.ac['aircraft_specs']['Vh']*self.ac['general_coefficients']['at']*(1-self.ac['general_coefficients']['epsilon_a'])
        return cm_a

    # Epsilon @ AoA = 0
    def epsilonZero(self): 
        epsilon_0 = 2*self.ac['general_coefficients']['CL_0w'] / (self.ac['aircraft_specs']['ARw']*math.pi) 
        return round(epsilon_0, 4)

    # Epsilon In Function of AoA - dE/alpha - Downwash - 1/rad
    def epsilonAlpha(self):
        epsilon_alpha = 2*self.ac['general_coefficients']['aw'] / (self.ac['aircraft_specs']['ARw']*math.pi)
        return round(epsilon_alpha, 4)

    def liftCoefZero(self):
        CL_0 = self.ac['general_coefficients']['aw']*(self.ac['aircraft_specs']['iw'] - self.ac['aircraft_specs']['alpha_0w']) + self.ac['aircraft_specs']['n_ef']*(self.ac['aircraft_specs']['St']/self.ac['aircraft_specs']['Sw'])*self.ac['general_coefficients']['at']*(self.ac['aircraft_specs']['it'] - self.ac['general_coefficients']['epsilon_0'])
        return round(CL_0, 4)

    def alphaZero(self):
        alpha_0 = -self.ac['general_coefficients']['CL_0'] / self.ac['general_coefficients']['aw']
        return round(alpha_0, 4)

    # Cornell - Eq. 3.17
    def cmZero(self, alpha_0):
        cm_0 = self.ac['general_coefficients']['Cm_0w'] - self.ac['aircraft_specs']['n_ef']*self.ac['aircraft_specs']['Vh']*self.ac['general_coefficients']['at']*(self.ac['aircraft_specs']['it'] - self.ac['general_coefficients']['epsilon_0'] + (1-self.ac['general_coefficients']['epsilon_a'])*alpha_0)
        return round(cm_0, 4)

    # E and R
    def alphaEquilibrium(self):
        alpha_e = (-self.ac['general_coefficients']['Cm_0'] / self.ac['general_coefficients']['Cm_a']) - abs(self.ac['aircraft_specs']['alpha_0'])
        return round(alpha_e, 4)

    def test(self):
        return self.ac['general_coefficients']['epsilon_0'] + self.ac['general_coefficients']['epsilon_a']
