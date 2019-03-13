#import
from math import exp, sqrt
from globals import GEAR, CONHIS, SPOIL, u, x, V, uInc, tuHis, deluHis, TrimHist, SMI, MODEL, RUNNING
from Atoms import Atoms


# global ??

    
def AeroModelUser(x,u,Mach,alphar,betar,V):
    '''
        T-X Trainer Aerodynamic Coefficients, Thrust Model,
        and Geometric and Inertial Properties for FLIGHT.m
        Low-Angle-of-Attack, Mach-Dependent Model

        Called by:
                EoM.py
                EoMQ.py
    '''

    global m
    global Ixx
    global Iyy
    global Izz
    global Ixz
    global S
    global b
    global cBar

    # Mass, Inertial, and Reference Properties
    m               =   4800
    Ixx             =   20950
    Iyy             =   49675
    Izz             =   62525
    Ixz             =   -1710
    cBar            =   3.03
    b               =   10
    S               =   27.77
    lHT             =   5.2
    lVT             =   3.9
    StaticThrust    =   49000

    (airDens,airPres,temp,soundSpeed) = Atoms(-x[5])
    Thrust   =   u[3]*StaticThrust*(airDens/1.225)*(1 - exp((-x[5] - 18000)/2000))
    
    # Mach Number Effect on All Incompressible Coefficients
    Prandtl  = 1 / sqrt(1 - Mach**2)	 # Prandtl Fac
    
    # Current Longitudinal Characteristics
    # ====================================
    # Lift Coefficient
    CLo     =   0
    CLa     =   4.92
    CLqhat  =   2.49
    CLq     =	CLqhat*cBar/(2*V)
    CLdE    =   0.72
    CLdS    =   CLdE

    # Total Lift Coefficient, w/Mach Correction
    CL      =	(CLo + CLa*alphar + CLq*x[7] + CLdS*u[6] + CLdE*u[0])*Prandtl

    # Drag Coefficient
    CDo     =   0.019
    Epsilon =   0.093*Prandtl
        
    # Total Drag Coefficient, w/Mach Correction
    CD      =	CDo*Prandtl + (Epsilon * CL**2)
                
    # Pitching Moment Coefficient
    StaticMargin    =   0.2
    Cmo     =   0
    Cma 	=	-CLa*StaticMargin
    Cmqhat  =   -4.3
    Cmq     =	Cmqhat*cBar/(2*V)
    CmV     =   0
    CmdE	=	-1.25
    CmdS 	=	CmdE
        
    # Total Pitching Moment Coefficient, w/Mach Correction
    Cm      =	(Cmo + Cma*alphar + Cmq*x[7] + CmdS*u[6] + CmdE*u[0])*Prandtl

    # Current Lateral-Directional Characteristics
    # ===========================================
    # Side-Force Coefficient
    CYo     =   0
    CYb     =   -0.5
    CYphat  =   0
    CYp     =   CYphat*(b/(2*V))
    CYrhat  =   0
    CYr     =   CYrhat*(b/(2*V))
    CYdA	=	0
    CYdR	=	0.04
        
    # Total Side-Force Coefficient, w/Mach Correction
    CYo     =   0
    CY      =	(CYo + CYb*betar + CYdR*u[2] + CYdA*u[1] + CYp*x[6] + CYr*x[8])*Prandtl

    # Rolling Moment Coefficient
    Clo     =   0
    Clb  	=	-0.066
    Clphat  =   -0.5
    Clp 	=	Clphat*(b/(2*V))				
    Clrhat  =   -0.5
    Clr 	=	Clrhat* (b/(2*V))
    CldA 	=	0.12
    CldR 	=	0.03
        
    # Total Rolling-Moment Coefficient, w/Mach Correction
    Cl      =	(Clo + Clb*betar + Clp * x[6] + Clr * x[8] + CldA*u[1] + CldR*u[2])* Prandtl

    # Yawing Moment Coefficient
    Cno     =   0
    CnBeta  =	0.37
    Cnphat  =   -0.06
    Cnp 	=	Cnphat*(b/(2*V))				
    Cnrhat  =   -0.5
    Cnr 	=	Cnrhat*(b/(2*V))				
    CndA 	=	0
    CndR 	=	-0.2
        
    # Total Yawing-Moment Coefficient, w/Mach Correction
    Cn      =	(Cno + CnBeta*betar + Cnp*x[6] + Cnr*x[8] + CndA*u[1] + CndR*u[2])*Prandtl


    return (CD,CL,CY,Cl,Cm,Cn,Thrust)
