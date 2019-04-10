### B747 Longitudinal Dynamics
## AA271a
## B747 at Mach 0.8, 40000ft, level-flight
## From Etkin and Reid page 166

from math import cos, sin, sqrt
from numpy import matrix, diag
from numpy.linalg import eig
units = "Metric" #or "English"
if units == "Metric":
	Xu = -1.982e3
	Xw = 4.025e3
	Zu = -2.595e4
	Zw = -9.030e4
	Zq = -4.524e5
	Zwd = 1.909e3
	Mu = 1.593e4
	Mw = -1.563e5
	Mq = -1.521e7
	Mwd = -1.702e4

	g, theta0, S, cbar = 9.81, 0, 511, 8.324
	U0, Iyy, m, cbar, rho = 235.9, 0.449e8, 2.83176e6/g, 8.324, 0.3045
	Xdp, Zdp, Mdp = 0.3*m*g, 0, 0
	Xde = -3.818e-6*(0.5 * rho * S * U0**2)
	Zde = -0.3648*(0.5 * rho * S * U0**2)
	Mde = -1.444*(0.5 * rho * S * cbar * U0**2)

else:
    ## English units
	pass

A = [
	[                 Xu/m,                  Xw/m,                            0,                     -g*cos(theta0)],
	[           Zu/(m-Zwd),            Zw/(m-Zwd),              (Zq+m*U0)/(m-Zwd),       (-m*g*sin(theta0))/(m-Zwd)],
	[(Mu + Zu*Mwd/(m-Zwd))/Iyy, (Mw + Zw*Mwd/(m-Zwd))/Iyy, (Mq + Zq+m*U0*Mwd/(m-Zwd))/Iyy, (-m*g*sin(theta0)*Mwd) / (m-Zwd)*Iyy],
	[                    0,                     0,                            1,                               0,]
	]

print(A)

B = [
	[                    Xde/m,                     Xdp/m],
	[              Zde/(m-Zwd),               Zdp/(m-Zwd)],
	[(Mde+Zde*Mwd/(m-Zwd))/Iyy, (Mdp+Zdp*Mwd/(m-Zwd))/Iyy],
	[                        0,                         0]
	]

C = [
	[1, 0, 0],
	[0, 1, 0],
	[0, 0, 1]
	]

# sp - short period
wsp		= sqrt(-U0*Mw/Iyy)
zetasp	= -Mq/2*sqrt(-1/U0/Mw/Iyy)

# p - phugoid
wp 	= sqrt(Zw*Mq/m/Iyy-U0*Mw/Iyy)
zetap 	= -(Zw/m+Mq/Iyy+Mwd/Iyy*U0)/2/wp

# EigenValues - lmda
lmda, lmdaVector = eig(matrix(A))
V = lmda
D = diag(lmdaVector)

# nondimensionalization of the eignvectors
lmda[0] = lmda[0] / U0
lmda[1] = lmda[1] / U0
lmda[2] = lmda[2] / (2*U0/cbar)

