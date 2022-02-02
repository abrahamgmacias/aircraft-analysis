from trim import Aerodynamics as aero

rangoDeVelocidades = []

# Se extraen los valores del CL y CD del aeronave completa via X5
# 

for v in rangoDeVelocidades:
    cLift = aero.liftCoefficient(wto, rho, v, Sw)
    cDrag = aero.dragCoefficient(CD0, CL, CL_min_D, ew, ARw)

    # Def these functions below...
    alpha = aero.Closest(CL_list, Angles_list, CL)

    aeroEfficiency = 
    tRequired = aero.TReq(wto, CD, CL)
    pRequired = aero.PReq(TR, v)

    pAvailable = PA_max*n_ef*(rho / rho_list[0])
    rClimb = af.R_C(PA, PR, wto)

    # Appending section

# Plotting section

