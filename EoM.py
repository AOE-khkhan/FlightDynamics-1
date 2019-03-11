from numpy import matmul, array, sqrt, arctan, arcsin, pi, iterp, sin, cos
#from event import event   ## write event fucnction
from DCM import DCM
from Atoms import Atoms
from WindField import WindField

## global ?????

def EoM(t,x):
    '''
        FLIGHT Equations of Motion
    '''
    if MODEL == 0:
        from AeroModelAlpha import AeroModelAlpha as AeroModel
    elif MODEL == 1:
        from AeroModelMach import AeroModelMach as AeroModel
    else:
        from AeroModelUser import AeroModelUser as AeroModel

    ##### Event Function ####### ????????????????
    ############################

    # Earth-to-Body-Axis Transformation Matrix
    HEB = DCM(x[9],x[10],x[11])
    # Atmospheric States
    x[5] = min(x[5], 0)     #Limit x[5] to <=0
    (airDens,airPres,temp,soundSpeed) = Atoms(-x[5])
    # Body-Axis Wind Field
	windb	=	WindField(-x[5],x[9],x[10],x[11])
    # Body-Axis Gravity Components
	gb		=	matmul(HEB * [0;0;9.80665])

    # Air-Relative Velocity Vector
    x[0]    =   max(x[0],0)     # Limit axial velocity to >= 0 m/s
    Va		=   array([x[0];x[1];x[2]]) + windb
	V		=	sqrt(matmul(Va', Va))

	alphar	=	arctan(Va[2] / abs(Va[0])
    #alphar  =   min(alphar, (pi/2 - 1e-6))     # Limit angle of attack to <= 90 deg
	alpha 	=	57.2957795 * alphar;
	betar	= 	arcsin(Va[1] / V);
	beta	= 	57.2957795 * betar;
	Mach	= 	V / soundSpeed;
	qbar	=	0.5 * airDens * V**2;

    # Incremental Flight Control Effects
    if CONHIS >=1 and RUNNING ==1:
        ## uInc = array([])
        uInc    =   iterp(t, tuHis,delHis)
        uInc    =   array(uInc).T   # Transpose
        uTotal  =   u + uInc
    else:
        uTotal  =   u
    
    # Force and Moment Coefficients; Thrust
    (CD,CL,CY,Cl,Cm,Cn,Trust) = AeroModel(x,uTotal,Mach,alphar,betar,V)

    qbarS   =   qbar * S

	CX	=	-CD * cos(alphar) + CL * sin(alphar)	# Body-axis X coefficient
	CZ	= 	-CD * sin(alphar) - CL * cos(alphar)	# Body-axis Z coefficient 

    # State Accelerations
	Xb =	(CX * qbarS + Thrust) / m
	Yb =	CY * qbarS / m
	Zb =	CZ * qbarS / m
	Lb =	Cl * qbarS * b
	Mb =	Cm * qbarS * cBar
	Nb =	Cn * qbarS * b
	nz	=	-Zb / 9.80665               # Normal load factor

    # Dynamic Equations
	xd1 = Xb + gb[0] + x[8] * x[1] - x[7] * x[2]
	xd2 = Yb + gb[1] - x[8] * x[0] + x[6] * x[2]
	xd3 = Zb + gb[2] + x[7] * x[0] - x[6] * x[1]
	
	y	=	matmul(HEB', [x[0];x[1];x[2]])
	xd4	=	y[0]
	xd5	=	y[1]
	xd6	=	y[2]
	
	xd7	= 	(Izz * Lb + Ixz * Nb - (Ixz * (Iyy - Ixx - Izz) * x[6] + (Ixz**2 + Izz * (Izz - Iyy)) * x[8]) * x[7]) / (Ixx * Izz - Ixz**2)
	xd8 = 	(Mb - (Ixx - Izz) * x[6] * x[8] - Ixz * (x[6]**2 - x[8]**2)) / Iyy
	xd9 =	(Ixz * Lb + Ixx * Nb + (Ixz * (Iyy - Ixx - Izz) * x[8] + (Ixz**2 + Ixx * (Ixx - Iyy)) * x[6]) * x[7]) / (Ixx * Izz - Ixz**2)

	cosPitch	=	cos(x[10]);
	if abs(cosPitch)	<=	0.00001
		cosPitch	=	0.00001 * (abs(cosPitch)/cosPitch)      # sign(cosPitch) == (abs(cosPitch)/cosPitch) for python
	end
	tanPitch	=	sin(x[10]) / cosPitch
		
	xd10	=	x[6] + (sin(x[9]) * x[7] + cos(x[9]) * x[8]) * tanPitch
	xd11	=	cos(x[9]) * x[7] - sin(x[9]) * x[8]
	xd12	=	(sin(x[9]) * x[7] + cos(x[9]) * x[8]) / cosPitch
	
	xdot	=	array([xd1,xd2,xd3,xd4,xd5,xd6,xd7,xd8,xd9,xd10,xd11,xd12])    


    return xdot
