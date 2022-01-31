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
