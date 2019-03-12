from numpy import matmul, array, sqrt, arctan, arcsin, pi, interp, sin, cos, reshape, matrix, zeros, ravel
#from event import event   ## write event fucnction
from DCM import DCM
from Atoms import Atoms
from WindField import WindField
#from globals import GEAR, CONHIS, SPOIL, u, x, V, uInc, tuHis, deluHis, TrimHist, SMI, MODEL, RUNNING


def EoM(t,x):
	'''
		FLIGHT Equations of Motion
	'''

	global m 
	global Ixx 
	global Iyy
	global Izz
	global Ixz
	global S
	global b
	global cBar
	global CONHIS
	global u
	global tuHis
	global deluHis
	global uInc
	global MODEL
	global RUNNING

	MODEL = 2
	CONHIS = 1
	RUNNING = 1
	tuHis	=	array([0, 33, 67, 100])
	deluHis	=	array(zeros(28)).reshape((4,7))
	u = array([0,0,0,0,0,0,0]).reshape((7,1)).ravel()

	print(f'tuHis = {tuHis}')
	print(f'deluHis = {deluHis}')
	print(f'u = {u}')

	print(f'x = {x}')

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
	gb = matmul(HEB, array([0,0,9.80665]).reshape((3,1))).ravel()

	print(f'windb = {windb}')


	# Air-Relative Velocity Vector
	x[0]    =   max(x[0],0)     # Limit axial velocity to >= 0 m/s
	Va		=   array([[x[0],x[1],x[2]]]).reshape(3,1).ravel() + windb

	print(f'Va 1st part = {array([[x[0],x[1],x[2]]]).reshape(3,1).ravel()}')
	print(f'windb.T = {matrix.transpose(windb)}')
	print(f'Va = {Va}')

	V		=	sqrt(matmul(matrix.transpose(Va), Va))

	print(f'V = {V}')

	alphar = arctan(Va[2]/abs(Va[0]))
	#alphar  =   min(alphar, (pi/2 - 1e-6))     # Limit angle of attack to <= 90 deg
	#alpha	=	57.2957795 * float(alphar)
	alpha = 57.2957795 * alphar
	betar	= 	arcsin(Va[1] / V)
	beta	= 	57.2957795 * betar
	Mach	= 	V / soundSpeed
	qbar	=	0.5 * airDens * V**2

	print(f'Mach = {Mach}')


	# Incremental Flight Control Effects
	if CONHIS >=1 and RUNNING ==1:
		## uInc = array([])
		uInc    =   interp(t, tuHis,deluHis[:, 0])
		uInc    =   matrix.transpose(array(uInc))   # Transpose
		uTotal  =   u + uInc
	else:
		uTotal  =   u
	
	# Force and Moment Coefficients; Thrust
	(CD,CL,CY,Cl,Cm,Cn,Thrust) = AeroModel(x,uTotal,Mach,alphar,betar,V)
	print(f'CD = {CD}')

	m               =   4800
	Ixx = 20950
	Iyy=49675
	Izz = 62525
	Ixz = -1710
	cBar            =   3.03
	b               =   10
	S               =   27.77
	lHT             =   5.2
	lVT             =   3.9
	StaticThrust    =   49000

	qbarS   =   qbar * S

	CX	=	-CD * cos(alphar) + CL * sin(alphar)	# Body-axis X coefficient
	CZ	= 	-CD * sin(alphar) - CL * cos(alphar)	# Body-axis Z coefficient 

	# State Accelerations
	Xb =	(CX * qbarS + Thrust) / m
	print(f'CX = {CX}')
	print(f'qbarS = {qbarS}')
	print(f'Thrust = {Thrust}')
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

	print(f'Xb = {Xb}')
	print(f'gb[0] = {gb[0]}')

	print(f'xd1 = {xd1}')

	# xd1 = xd1[0][0]
	# xd2 = xd2[0][0]
	# xd3 = xd3[0][0]

	HEB_T = matrix.transpose(HEB)
	y = matmul(HEB_T, array([x[0], x[1], x[2]]))
	#HEB_T = matrix.transpose(HEB)
	#y = matmul(HEB_T, (array(x[0],x[1],x[2]))


	xd4 = y[0]
	xd5 = y[1]
	xd6 = y[2]

	xd7	= 	((Izz * Lb + Ixz * Nb - (Ixz * (Iyy - Ixx - Izz) * x[6] + (Ixz**2 + Izz * (Izz - Iyy)) * x[8]) * x[7]) / (Ixx * Izz - Ixz**2))
	xd8 = 	((Mb - (Ixx - Izz) * x[6] * x[8] - Ixz * (x[6]**2 - x[8]**2)) / Iyy)
	xd9 =	((Ixz * Lb + Ixx * Nb + (Ixz * (Iyy - Ixx - Izz) * x[8] + (Ixz**2 + Ixx * (Ixx - Iyy)) * x[6]) * x[7]) / (Ixx * Izz - Ixz**2))

	# xd7 = xd7[0][0]
	# xd8 = xd8[0][0]
	# xd9 = xd9[0][0]

	cosPitch	=	cos(x[10])
	if abs(cosPitch)	<=	0.00001:
		cosPitch	=	0.00001 * (abs(cosPitch)/cosPitch)      # sign(cosPitch) == (abs(cosPitch)/cosPitch) for python
	tanPitch	=	sin(x[10]) / cosPitch
		
	xd10	=	x[6] + (sin(x[9]) * x[7] + cos(x[9]) * x[8]) * tanPitch
	xd11	=	cos(x[9]) * x[7] - sin(x[9]) * x[8]
	xd12	=	(sin(x[9]) * x[7] + cos(x[9]) * x[8]) / cosPitch

	
	xdot	=	array([xd1,xd2,xd3,xd4,xd5,xd6,xd7,xd8,xd9,xd10,xd11,xd12])

	for i in range(1,13):
		print(f'xd{str(i)} = {xdot[i-1]} ') 

	return xdot

if __name__ == "__main__":
	print(EoM(1,[1,2,3,1,2,3,1,2,3,1,2,3]))   # Test

