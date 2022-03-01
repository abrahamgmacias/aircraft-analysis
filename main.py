from resources import aircraft_functions as af

# --------- Extra Computation ---------- # 
# Class Instantiation
general_aircraft = af.Aircraft()
general_aircraft.addWeights(mtow=22)

# Data appending for trim_analysis.py
general_aircraft.addCoefficients(cl_list=[], cd_list=[])      # Pull from excel file...
general_aircraft.addCoefficients(cd_0=0.03324, cl_min_d=0, cl_0=0.6096)

wing = af.Wing(3.71, 0.494, 'MH-114')
wing.addGeometry(arw=7.5, sw=1.831, cw=0.494, iw=round(1/57.3, 4))
wing.addCoefficients(ew=0.825)

propeller = af.Propeller(12, 12)   # placeholder data
motor = af.Motor(propeller)

motor.addMotorSpecs(pa_max=32)   # placeholder pa_max
propeller.addPropellerSpecs(n_eff=0.95)

atmosphere = af.Atmospheric(1.225)

# Data appending for static_stability.py
# Must create a tail class... can modify the wing class to morph into a tail
horizontal_stab = af.Wing(description='Horizontal Stabilizer')
horizontal_stab.addGeometry(st=0.37, lt=1.48, it=round(-2.5/57.3, 4))
horizontal_stab.addCoefficients(vh=0.6, at=3.6385)

wing.addCoefficients(aw=4.4063, ew=0.825, cm_0w=-0.208, cl_0w=0.62, alpha_0w=round(-7.8/57.3, 4))

general_aircraft.addCoefficients(x_cg_cw=0.2834, x_ac_cw=0.25, cl_0=0.6096, cd_0=0.03324)
general_aircraft.addComponents(wing=wing, motor=motor, h_stabilizer=horizontal_stab)
