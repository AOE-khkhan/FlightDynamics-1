### FLIGHT.py -- Generic 6-DOF Trim, Linear Model, and Flight Path Simulation

#####   ????? <== search for queries and doubts

#import
import numpy as np
from Atoms import Atoms


## global variables
GEAR        = 0
CONHIS      = 0
SPOIL       = 0
u           = 0
x           = 0
V           = 0
uInc        = 0
tuHis       = 0
deluHis     = 0
TrimHist    = 0
SMI         = 0
MODEL       = 0
RUNNING     = 0

print('** ======================= **')
print('** 6-DOF FLIGHT Simulation **')
print('** ======================= **\n')

#### Main File. It contain ...
#       Define initial conditions
#       Contain aerodynamic data tables(if required)
#       Calculates longitudinal trim condition
#       Calculate stability & control matrices for linear model
#       Simulate flight path using nonlinear equation of motion

#### Function used by FLIGHT:
#       AeroModelAlpha.m    High-Alpha, Low-Mach aerodynamic coefficients of the aircraft,
#                           thrust model, and geometric and inertial properties
#       AeroModelMach.m     Low-Alpha, High-Mach aerodynamic coefficients of the aircraft,
#                           thrust model, and geometric and inertial properties
#       AeroMoedelUser.m    User-defined aerodynamic corfficietns of the aircraft,
#                           thrust model, and geometric and inertial properties
#       Atoms.m             Air density, sound speed
#       DCM.m               Direction-cosine (Rotation) matrix from Euler Angles
#       RMQ.m               Direction-cosine (Rotation) matrix from Quaternions
#       EoM.m               Equations of motion for integration (Euler Angles)
#       EoMQ.m              Equations of motion for integration (Quaternions)
#       Linmodel.m          Equations of motion for defining linear-model
#                           (F & G) matrices via central differences
#       TrimCost.m          Cost function for trim solution
#       WindField.m         Wind velocity components

### DEFINITION OF THE STATE VECTOR
##   With Euler Angle DCM option (QUAT=0):
#       x[0]    =   Body-axis x in inertial velocity, ub, m/s
#       x[1]    =   Body-axis y in inertial velocity, vb, m/s
#       x[2]    =   Body-axis z in inertial velocity, wb, m/s
#       x[3]    =   North position of centre of mass w.r.t Earth, xe, m
#       x[4]    =   East position of centre of mass w.r.t Earth, ye, m
#       x[5]    =   Negative position of centre of mass w.r.t Earth, ze = -h, m
#       x[6]    =   Body-axis roll rate, pr, rad/s
#       x[7]    =   Body-axis pitch rate, qr, rad/s
#       x[8]    =   Body-axis yaw rate, rr, rad/s
#       x[9]    =   Roll angle of body w.r.t Earth, phir, rad
#       x[10]   =   Pitch angle of body w.r.t Earth, qhir, rad
#       x[11]   =   Yaw angle of body w.r.t Earth, psir, rad 
##  With Quaternions DCM option (QUAT=1):
#       ## WILL DO THIS LATER ##


### DEFINITION OF THE CONTROL VECTOR
#       u[0]    =   Elevator, dEr, rad, positive: trailing edge down
#       u[1]    =   Aileron, dAr, rad, positive: left trailing edge down
#       u[2]    =   Rudder, dRr, rad, positive: trailing edge left
#       u[3]    =   Throttle, dT
#       u[4]    =   Asymmetric Spoiler, dASr, rad
#       u[5]    =   Flap, dFr, rad
#       u[6]    =   Stabilator, dSr, rad


### ================================================================================
##  USER INPUT
##  ============
##  --  Flag define which analysis or Input conditions(I.C.)/Inputs will be engaged
##  ===============================================================
##  FLIGHT Flags (1 = ON, 0 = OFF)

MODEL   =   2       # Aerodynamic model selection
                    # 0: Incompressible flow, high angle of attack  ## WILL DO THIS LATER ##
                    # 1: Compressible flow, low angle of attack  ## WILL DO THIS LATER ##
                    # 2: User-Defined model
QUAT    =   0       # 0: Rotation matrix (DCM) from Euler Angles
                    # 1: Rotation matrix (DCM) from Quaternions  ## WILL DO THIS LATER ##
TRIM    =   1       # Trim Flag (=1 to calculate trim at Intial Condition)
LINEAR  =   1       # Linear model flag (=1 to calculate and store F and G)
SIMUL   =   1       # Flight path flag (=1 for nonlinear simulation)
CONHIS  =   1       # Control History ON(=1) or OFF(=0)
RUNNING =   0       # internal flag
GEAR    =   0       # Landing gear DOWN(=1) or UP(=0)
SPOIL   =   0       # Symmetric Spoiler DEPLOYED(=1) or CLOSED(=0)
dF      =   0       # Flap setting, deg

### Initial Altitude(ft), Indicated Airspeed(kt)
hft     =   10000   # Altitude above Sea Level, ft
VKIAS   =   150     # Indicated Airspeed, kt

hm      =   hft * 0.3048    # Altitude above Sea Level, m
VmsIAS  =   VKIAS * 0.5154  # Indicated Airspeed, m/s

print('Initial Conditions')
print('==================')
print(f'Altitude           = {str(hft)} ft,   = {str(hm)} m')
print(f'Indicated Airspeed = {str(VKIAS)} kt,   = {str(VmsIAS)} m/s')

### US Standard Atmosphere, 1976, Table lookup for I.C.
(airDens, airPres, temp, soundSpeed) = Atoms(hm) 
print(f'Air Density     = {str(airDens)} kg/m**3, Air Pressure = {str(airPres)} N/m**2')
print(f'Air Temperature = {str(temp)} K,         Sound Speed  = {str(soundSpeed)} m/s')
        
# Dynamic Pressure (N/m**2), and True Airspeed (m/s)
qBarSL  =   0.5*1.225 * VmsIAS**2      # Dynamic Pressure at sea level, N/m**2
V       =   np.sqrt(2*qBarSL/airDens)	# True Airspeed, TAS, m/s
TASms   =   V
print(f'Dynamic Pressure = {str(qBarSL)} N/m**2, True Airspeed = {str(V)} m/s')

### Alphabetical List of Initial Conditions

alpha   =	0      # Angle of attack, deg (relative to air mass)
beta    =	0      # Sideslip angle, deg (relative to air mass)
dA      =	0      # Aileron angle, deg
dAS     =	0      # Asymmetric spoiler angle, deg
dE      =	0      # Elevator angle, deg
dR      =	0      # Rudder angle, deg
dS      = 	0      # Stabilator setting, deg
dT      = 	0      # Throttle setting,   # / 100
hdot    =	0      # Altitude rate, m/s
p       =	0      # Body-axis roll rate, deg/s
phi     =	0      # Body roll angle wrt earth, deg
psi     =	0      # Body yaw angle wrt earth, deg
q       =	0      # Body-axis pitch rate, deg/sec
r       =	0      # Body-axis yaw rate, deg/s
SMI     =	0      # Static margin increment due to center-of-mass variation from reference, #/100
tf      =	100    # Final time for simulation, sec
ti      = 	0      # Initial time for simulation, sec
theta   =	alpha  # Body pitch angle wrt earth, deg [theta = alpha if hdot = 0]
xe      =	0      # Initial longitudinal position, m
ye      = 	0      # Initial lateral position, m
ze      = 	-hm    # Initial vertical position, m [h: + up, z: + down]

## Initial Conditions Depending on Prior Initial Conditions
gamma   =   57.2957795 * np.arctan(hdot / np.sqrt(V**2 - hdot**2))    ######  ?????? # Inertial Vertical Flight Path Angle, deg  # deg = 57.2957795 * rad
qbar	= 	0.5 * airDens * V**2	    # Dynamic Pressure, N/m**2
IAS		=	np.sqrt(2 * qbar / 1.225)      # Indicated Air Speed, m/s
Mach	= 	V / soundSpeed              # Mach Number
print(f'Mach number = {str(Mach)}, Flight Path Angle = {str(gamma)}, deg\n')

uInc    =   [];      ######  ??????

if MODEL == 0:
    print('<<Low-Mach, High-Alpha Model>>')
elif MODEL == 1:
    print('<<High-Mach, Low-Alpha Aerodynamic Model>>')
else:
    print('<<USer Defined AeroModel>>')
print('  ======================================')

## Initial Control Perturbation (Test Inputs: rad or 100#)			
delu	=	np.array([0,0,0,0,0,0,0])  # column vector  # delta u   in u = u0 + delta u
## Initial State Perturbation (Test Inputs: m, m/s, rad, or rad/s)
delx	=	np.array([0,0,0,0,0,0,0,0,0,0,0,0])  ##.reshape((4,3)) ###### ????? SHAPE ????   # delta x   in x = x0 + delta x

print('Initial Perturbations to Trim for Step Response')
print('===============================================')

print('Control Vector')
print('--------------')
print(f'Elevator   = {str(delu[0])} rad, Aileron = {str(delu[1])} rad, Rudder = {str(delu[2])} rad')
print(f'Throttle   = {str(delu[3])} x 100%, Asymm Spoiler = {str(delu[4])} rad, Flap = {str(delu[5])} rad')
print(f'Stabilator = {str(delu[6])} rad\n')

print('State Vector')
print('------------')
print(f'u   = {str(delx[0])} m/s, v = {str(delx[1])} m/s, w = {str(delx[2])} m/s')
print(f'x   = {str(delx[3])} m, y = {str(delx[4])} m, z = {str(delx[5])} m')
print(f'p   = {str(delx[6])} rad/s, q = {str(delx[7])} rad/s, r = {str(delx[8])} rad/s')
print(f'Phi = {str(delx[9])} rad, Theta = {str(delx[10])} rad, Psi = {str(delx[11])} rad')

