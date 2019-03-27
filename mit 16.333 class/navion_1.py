### Calculate the Static Stability Parameters for the Navion General Aviation Aircraft
### Adapted from Nelson, page 57

import matplotlib.pyplot as plt

pi = 3.14
## Flight Condition
W = 2750        # lb
V = 176         # ft/sec
x_cg = 0.295    # cbar
h_cg = 0.295

## Reference Geometry
S = 184           # ft^2
b = 33.4          # ft
cbar = 5.7        # ft
S_t = 43          #ft^2
l_t = 16          #ft

## Wing Airfoil characteristics
C_m_ac_w = -0.116
C_l_alp_w = 0.097 #/deg
alp_0_L = -5      # deg
X_ac = 0.25       #cbar
h_ac = 0.25
i_w = 1           #deg

## Tail Airfoil section
C_l_alp_t = 0.088 #/deg
C_m_ac_t = 0
i_t = -1         # deg
eta = 1;
b_t = (2.2/5.5)*b

##
V_H = l_t*S_t/cbar/S
AR_w = b**2/S
AR_t = b_t**2/S_t
d2r = 180/pi      # degree to radians

## convert wing section lifts to 3D wing
C_L_alp_w = C_l_alp_w*d2r / (1 + (C_l_alp_w*d2r/(pi*AR_w)))
C_L_alp_t = C_l_alp_t*d2r / (1 + (C_l_alp_t*d2r/(pi*AR_t)))

## Wing Contribution
# C_m_cg_w = C_L_w(h - h_nbar) + C_m_ac_w
#
# with C_L_w = C_L_0_w + C_L_alp_w*alp_w
# where C_L_0_w = lift coeff at zero angle of attack
C_L_0_w = C_L_alp_w * abs(alp_0_L / d2r)
#
C_m_cg_0_w = C_m_ac_w +  C_L_0_w*(h_cg - h_ac)
C_m_cg_alp_w = C_L_alp_w*(h_cg - h_ac)

if C_m_cg_0_w > 0 and C_m_cg_alp_w < 0 :
    print("Wing's contribution to C_m_alpha (Pitch Static Stability), for this particular airplane, is 'destabilizing'.")
else:
    print("Wing's contribution to C_m_alpha (Pitch Static Stability), for this particular airplane, is 'destabilizing'.")


## For Downwash
# eps = eps_0 + eps_alp*alp_w
eps_0 = 2*C_L_0_w / (pi*AR_w) * d2r
eps_alp = 2*C_L_alp_w / (pi*AR_w)

## Tail Contribution
C_m_cg_0_t = (eta*V_H*C_L_alp_t*(eps_0 + i_w - i_t))/d2r
C_m_cg_alp_t = -eta*V_H*C_L_alp_t*(1 - eps_alp)

if C_m_cg_0_t > 0 and C_m_cg_alp_t < 0 :
    print("Tail's contribution to C_m_alpha (Pitch Static Stability), for this particular airplane, is 'destabilizing'.")
else:
    print("Tail's contribution to C_m_alpha (Pitch Static Stability), for this particular airplane, is 'destabilizing'.")


## Total Contribution to Pitch Stability
C_m_cg_0 = C_m_cg_0_w + C_m_cg_0_t
C_m_cg_alp = C_m_cg_alp_w + C_m_cg_alp_t

if C_m_cg_0_w > 0 and C_m_cg_alp_w < 0 :
    print("This Aircraft has Positive Longitudinal Static Stability")
else:
    print("This Aircraft has Negative Longitudinal Static Stability")


## Figure
plt.figure('Contribution to Pitching Moment', figsize=(16,9))
alpha = list(range(0,13))
C_m_cg_wing_total = [C_m_cg_0_w + C_m_cg_alp_w*alpha[i]/d2r for i in alpha]
C_m_cg_tail_total = [C_m_cg_0_t + C_m_cg_alp_t*alpha[i]/d2r for i in alpha]
C_m_cg_total = [C_m_cg_0 + C_m_cg_alp*alpha[i]/d2r for i in alpha]
plt.plot(alpha,C_m_cg_wing_total,'b-', alpha,C_m_cg_tail_total,'g-', alpha,C_m_cg_total,'r-')
plt.xlabel('alpha, deg')
plt.ylabel('C_m_cg')
plt.legend(('Wing', 'Tail', 'Total Airplane'))
plt.title('C_m_cg  Vs  alpha')
plt.grid(True)
plt.show()

## Neutral point
h_np = h_ac + (eta*V_H*C_L_alp_t / C_L_alp_w)*(1 - eps_alp)
