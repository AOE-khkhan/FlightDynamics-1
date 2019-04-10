import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

##################################################################
################  Simulation Settings   ##########################
##################################################################

## Physical Constant
m   =   0.1             # Kg
Ixx =   0.00062         # Kg*m^2
Iyy =   0.00113         # Kg*m^2
Izz =   0.9*(Ixx+Iyy)   # Kg*m^2  (Assume neaerly flat object, z=0)
dx  =   0.114           # m   # distance of rotor from c.m. in x -direction
dy  =   0.0825          # m   # distance of rotor from c.m. in y -direction
g   =   9.81            # m/s/s
DTR =   1/57.3          # deg to rad
RTD =   57.3            # rad to deg


## Simulation time and model parameters
tstep = 0.02    # Sampling time, sec
Simulation_time = 30    # Length of time to run simulation, sec
t = np.arange(0,Simulation_time,tstep)  # time array
#print(t)
#print(np.shape(t))

## Model size
n_states = 12   # number of states
n_inputs = 4    # number of inputs(controls)

## States and Control Inputs
# x      = [u; v; w; p; q; r; phi; theta; psi; xE; yE; hE]
# xIndex = [0; 1; 2; 3; 4; 5;   6;     7;   8;  9; 10; 11]
# u      = [pitch; roll; yaw; climb]
# uIndex = [    0;    1;   2;     3]

## Initialize State Conditions
x = np.zeros((n_states,np.size(t))) # time history of state vectors
x[11,0] = 0.0    # Initial Height
#print(x)
#print(np.shape(x))

## Initialize inputs
u = np.zeros((n_inputs,np.size(t))) # tiem history of input vectors
u[:,0] = np.zeros(4)    # Initialize control inputs
#print(u)
#print(np.shape(u))


##################################################################
################  Aerodynamic Model   ############################
##################################################################

from scipy.optimize import fsolve

## Propeller Thrust equations as a function of propeller induced velocity, vi
def thrustEqn(vi, *prop_params):
    '''
        Input:-    Induced Velocity, vi
              -    propeller parameters = (R,A,rho,a,b,c,eta,theta0,theta1,U,V,W,Omega)

        Output:-   residual Thrust
    '''

    # Unpack parameters
    R,A,rho,a,b,c,eta,theta0,theta1,U,V,W,Omega = prop_params

    # Calculate local airflow velocity at propeller with vi, V'
    Vprime = np.sqrt(U**2 + V**2 + (W-vi)**2)

    # Calculate Total Thrust averaged over one revolution of propeller using vi
    Thrust = 1/4 * rho * a * b * c * R * \
        ( (W - vi) * Omega * R + 2/3 * (Omega * R)**2 * (theta0 + 3/4 * theta1) + \
            (U**2 + V**2) * (theta0 + 1/2 * theta1) )

    # Calculate residual for equation: Thrust = mass flow rate * delta Velocity
    residual = eta * 2 * vi * rho * A * Vprime - Thrust

    return residual

def Fthrust(x, u, dx, dy):
    '''
        Inputs:-    Current state x[k], Commanded Propeller RPM inputs u[k],
                    Propeller location distances dx, dy (m)
        Output:-    Thrust vector for 4 propellers (Newtons)
    '''

    # Propeller Configuration parameters
    R = 0.0762   # propeller length/ disk radius (m) 
    A = np.pi * R ** 2
    rho = 1.225  #kg/m^3  at MSL
    a = 5.7      # Lift curve slope
    b = 2        # number of blades
    c = 0.0274   # mean chord length (m)
    eta = 1      # propeller efficiency
    
    # Manufacturer propeller length x pitch specification:
    p_diameter = 6  #inches
    p_pitch = 3   #inches
    
    theta0 = 2*np.arctan2(p_pitch, (2 * np.pi * 3/4 * p_diameter/2))
    theta1 = -4 / 3 * np.arctan2(p_pitch, 2 * np.pi * 3/4 * p_diameter/2)
    
    
    # Local velocity at propeller from vehicle state information
    ub, vb, wb = x[0], x[1], x[2]
    p, q, r = x[3], x[4], x[5]
    # Transofrm velocity to local propeller location:
    #     [U,V,W] = [ub,vb,wb] + [p,q,r] x [dx,dy,0]
    U = ub - r * dy
    V = vb + r * dx
    W = wb - q * dx + p * dy
    
    # Convert commanded RPM to rad/s
    Omega = 2 * np.pi / 60 * u

    # Collect propeller config, state, and input parameters
    prop_params = (R,A,rho,a,b,c,eta,theta0,theta1,U,V,W,Omega)

    # Numerically solve for propeller induced velocity, vi
    # using nonlinear root finder, fsolve, and prop_params
    vi0 = 0.1    # initial guess for vi
    vi = fsolve(thrustEqn, vi0, args=prop_params)     # help(fsolve) for more info
    
    # Plug vi back into Thrust equation to solve for T
    Vprime = np.sqrt(U**2 + V**2 + (W-vi)**2)
    Thrust = eta * 2 * vi * rho * A * Vprime
    
    return Thrust

## Torque function
def T(F,dx,dy):
    '''
        Output:-    torque about cg given thrust force and 
                    dx, dy distance from cg
    '''
    ## Placeholder ##
    return 0

## Plot Thrust as a fuction of RPM for various vertical velocity conditions
RPM = np.linspace(1000,6000,200)
#print(RPM)
#print(np.shape(RPM))
vertvel = np.array([0,0,1] + 9*[0])   # w(vertical velocity) = 1, all others are zero
#print(vertvel)
#print(np.shape(vertvel))
Thrust_m2vel = np.array([Fthrust(2*vertvel,rpmIn,dx,dy) for rpmIn in RPM])
#print(Thrust_m2vel)
#print(np.shape(Thrust_m2vel))
Thrust_m1vel = np.array([Fthrust(1*vertvel,rpmIn,dx,dy) for rpmIn in RPM])
Thrust_0vel  = np.array([Fthrust(0*vertvel,rpmIn,dx,dy) for rpmIn in RPM])
Thrust_p1vel = np.array([Fthrust(-1*vertvel,rpmIn,dx,dy) for rpmIn in RPM])
Thrust_p2vel = np.array([Fthrust(-2*vertvel,rpmIn,dx,dy) for rpmIn in RPM])
fig = plt.figure(figsize=(8,8))
plt.plot(RPM, 4 * Thrust_m2vel / (m*g) )
plt.plot(RPM, 4 * Thrust_m1vel / (m*g) )
plt.plot(RPM, 4 * Thrust_0vel / (m*g) )
plt.plot(RPM, 4 * Thrust_p1vel / (m*g) )
plt.plot(RPM, 4 * Thrust_p2vel / (m*g) )
plt.plot(RPM, np.ones(np.size(RPM)), 'k--')
plt.legend(('Airspeed = -2 m/s','Airpseed = -1 m/s','Airspeed =  0 m/s', \
            'Airpseed =  1 m/s','Airspeed =  2 m/s'), loc='upper left')
plt.xlabel('Propeller RPM (x4)')
plt.ylabel('Thrust (g)')
plt.title('Quadcopter Thrust for different Vertical Airspeeds')
plt.show()


##################################################################
################  Equation of Motion   ###########################
##################################################################

## Nonlinear Dynamics Equations of Motion
def stateDerivative(x,u):
    '''
    Inputs:-    state vector (x), input vector (u)
    Outputs:-   time derivative of state vector (xdot)
    
    #  State Vector Reference:
    #idx  0, 1, 2, 3, 4, 5,  6,   7,   8,   9, 10, 11
    #x = [u, v, w, p, q, r, phi, the, psi, xE, yE, hE]
    '''

    # Store state variables in a readable format
    ub    = x[0]
    vb    = x[1]
    wb    = x[2]
    p     = x[3]
    q     = x[4]
    r     = x[5]
    phi   = x[6]
    theta = x[7]
    psi   = x[8]
    xE    = x[9]
    yE    = x[10]
    hE    = x[11]
    
    # Calculate forces from propeller inputs (u)
    F1 = Fthrust(x, u[0],  dx,  dy)
    F2 = Fthrust(x, u[1], -dx, -dy)
    F3 = Fthrust(x, u[2],  dx, -dy)
    F4 = Fthrust(x, u[3], -dx,  dy)
    Fz = F1 + F2 + F3 + F4
    L = (F2 + F3) * dy - (F1 + F4) * dy
    M = (F1 + F3) * dx - (F2 + F4) * dx
    N = -T(F1,dx,dy) - T(F2,dx,dy) + T(F3,dx,dy) + T(F4,dx,dy)
    
    # Pre-calculate trig values
    cphi = np.cos(phi);   sphi = np.sin(phi)
    cthe = np.cos(theta); sthe = np.sin(theta)
    cpsi = np.cos(psi);   spsi = np.sin(psi)
    
    # Calculate the derivative of the state matrix using EOM
    xdot = np.zeros(12)
    
    xdot[0] = -g * sthe + r * vb - q * wb  # = udot
    xdot[1] = g * sphi*cthe - r * ub + p * wb # = vdot
    xdot[2] = 1/m * (-Fz) + g*cphi*cthe + q * ub - p * vb # = wdot
    xdot[3] = 1/Ixx * (L + (Iyy - Izz) * q * r)  # = pdot
    xdot[4] = 1/Iyy * (M + (Izz - Ixx) * p * r)  # = qdot
    xdot[5] = 1/Izz * (N + (Ixx - Iyy) * p * q)  # = rdot
    xdot[6] = p + (q*sphi + r*cphi) * sthe / cthe  # = phidot
    xdot[7] = q * cphi - r * sphi  # = thetadot
    xdot[8] = (q * sphi + r * cphi) / cthe  # = psidot
    xdot[9] = cthe*cpsi*ub + (-cphi*spsi + sphi*sthe*cpsi) * vb + \
        (sphi*spsi+cphi*sthe*cpsi) * wb  # = xEdot
    xdot[10] = cthe*spsi * ub + (cphi*cpsi+sphi*sthe*spsi) * vb + \
        (-sphi*cpsi+cphi*sthe*spsi) * wb # = yEdot
    xdot[11] = -1*(-sthe * ub + sphi*cthe * vb + cphi*cthe * wb) # = hEdot
    
    return xdot


##################################################################
################  Control Law   ##################################
##################################################################

def controlInputs(x, t):
    '''
    Inputs:-    Current state x[k], time t
    Output:-    Control inputs u[k]
    
    #### Placeholder Function ####
    '''

    # Trim RPM for all 4 propellers to provide thrust for a level hover
    trim = 3200
    
    pitch_cmd = 0
    roll_cmd = 0
    climb_cmd = 0
    yaw_cmd = 0
    
    # Example open loop control inputs to test dynamics:
    #  Climb
    if t < 11.0:
        climb_cmd = 500
    
    #  Pitch Forward
    if t > 8.0:
        pitch_cmd = -10
    if t > 9.0:
        pitch_cmd = 10
    if t > 10.0:
        pitch_cmd = 0
    
    #  Pitch Backward
    if t > 12.0:
        pitch_cmd = 15
    if t > 13.0:
        pitch_cmd = -15
    if t > 14.0:
        pitch_cmd = 0
    
    #  Increase lift
    if t > 16.0:
        climb_cmd = 150
        
    
    # RPM command based on pitch, roll, climb, yaw commands
    u = np.zeros(4)
    #print(u)
    #print(np.shape(u))
    u[0] = trim + ( pitch_cmd + roll_cmd + climb_cmd - yaw_cmd) / 4
    u[1] = trim + (-pitch_cmd - roll_cmd + climb_cmd - yaw_cmd) / 4
    u[2] = trim + ( pitch_cmd - roll_cmd + climb_cmd + yaw_cmd) / 4
    u[3] = trim + (-pitch_cmd + roll_cmd + climb_cmd + yaw_cmd) / 4
    #print(u)
    
    return u


##################################################################
################  Linearized Equations   #########################
##################################################################

def jacobian(x, PHI_DOT_0, PSI_DOT_0, Ixx, Iyy, Izz):
    U0 = x[0]
    V0 = x[1]
    W0 = x[2]
    P0 = x[3]
    Q0 = x[4]
    R0 = x[5]
    PHI0 = x[6]
    THETA0 = x[7]
    PSI0 = x[8]
    #xE = x[9]   #not needed
    #yE = x[10]
    #hE = x[11]
    
    SPH = np.sin(PHI0)
    CPH = np.cos(PHI0)
    ST = np.sin(THETA0)
    CT = np.cos(THETA0)
    SPS = np.sin(PSI0)
    CPS = np.cos(PSI0)
    tthe = np.tan(THETA0)
    secthe = 1/CT
    
    
    return np.matrix([
        [0, R0, -Q0, 0, -W0, V0, 0, -g * np.cos(THETA0), 0, 0, 0, 0],
        [-R0, 0, P0, W0, 0, -U0, g * np.cos(PHI0) * np.cos(THETA0), -g * np.sin(PHI0) * np.sin(THETA0), 0, 0, 0, 0],
        [Q0, -P0, 0, -V0, U0, 0, g * np.sin(PHI0) * np.cos(THETA0), g * np.cos(PHI0) * np.sin(THETA0), 0, 0, 0, 0],
        [0, 0, 0, 0, (Iyy-Izz)/Ixx * R0, (Iyy-Izz)/Ixx * Q0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, (Izz-Ixx)/Iyy * R0, 0, (Izz-Ixx)/Iyy * P0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, (Ixx-Iyy)/Izz * Q0, (Ixx-Iyy)/Izz * P0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, tthe* np.sin(PHI0), tthe* np.cos(PHI0), tthe*(Q0 * np.cos(PHI0) - R0 * np.sin(PHI0)), (Q0 * np.sin(PHI0) + R0 * np.cos(PHI0) + (PHI_DOT_0 - P0) * tthe), 0, 0, 0, 0],
        [0, 0, 0, 0, np.cos(PHI0), -np.sin(PHI0), -(Q0 * np.sin(PHI0) + R0 * np.cos(PHI0)), 0, 0, 0, 0, 0],
        [0, 0, 0, 0, secthe* np.sin(PHI0), secthe* np.cos(PHI0), -secthe*(R0 * np.sin(PHI0) - Q0 * np.cos(PHI0)), PSI_DOT_0 *tthe, 0, 0, 0, 0],
        [CT*CPS,(SPH*ST*CPS - CPH*SPS),(SPH*SPS + CPH*ST*CPS), 0, 0, 0, (CPH*ST*CPS*V0 - SPH*SPS*V0 + CPH*SPS*W0 - SPH*ST*CPS*W0), (-ST*CPS*U0 + SPH*CT*CPS*V0 + CPH*CT*CPS*W0), (-CT*SPS*U0 + CPH*CPS*V0 - SPH*ST*SPS*V0 + SPH*CPS*W0 - CPH*ST*SPS*W0), 0, 0, 0],
        [CT*SPS, (CPH*CPS + SPH*ST*SPS), (CPH*ST*SPS - SPH*CPS), 0, 0, 0, (CPH*ST*SPS*V0 - SPH*CPS*V0 - CPH*CPS*W0 - SPH*ST*SPS*W0), (-ST*SPS*U0 + SPH*CT*SPS*V0 + CPH*CT*SPS*W0), (CT*CPS*U0 - CPH*SPS*V0 + SPH*ST*CPS*V0 + SPH*SPS*W0 + CPH*ST*CPS*W0), 0, 0, 0],
        [ST, SPH*CT, -CPH*CT, 0, 0, 0, (-CPH*CT*V0 + SPH*CT*W0), (-CT*U0 + SPH*ST*V0 + CPH*ST*W0), 0, 0, 0, 0]
    ])


Bmat = np.concatenate([np.zeros((2,4)), np.diag([-1/m,1/Ixx,1/Iyy,1/Izz]), np.zeros((6,4))])

print(Bmat)
print(f'Shape of Bmat ={np.shape(Bmat)} ')



##################################################################
################  Numerical Integration   ########################
##################################################################

## 4th Order Runge Kutta Calculation
def RK4(x,u,dt):
    '''
    Inputs:-    x[k], u[k], dt (time step, seconds)
    Output:-    x[k+1]
    '''

    # Calculate slope estimates
    K1 = stateDerivative(x, u)
    K2 = stateDerivative(x + K1 * dt / 2, u)
    K3 = stateDerivative(x + K2 * dt / 2, u)
    K4 = stateDerivative(x + K3 * dt, u)
    
    # Calculate x[k+1] estimate using combination of slope estimates
    x_next = x + 1/6 * (K1 + 2*K2 + 2*K3 + K4) * dt
    
    #print(f'xnext = {x_next}')
    #print(f'Shape of xnext = {np.shape(x_next)}')

    return x_next

## March through time array and numerically solve for vehicle states
for k in range(0, np.size(t) - 1): 
    # Determine control inputs based on current state
    u[:,k] = controlInputs(x[:,k], t[k])
    # Predict state after one time step
    x[:,k+1] = RK4(x[:,k], u[:,k], tstep)


## 4th Order Runge Kutta Calculation
def RK4_Linear(x,umat,Amat,Bmat,dt):
    '''
    Inputs:-    x[k], u[k], dt (time step, seconds)
    Output:-    x[k+1]
    '''
    # Calculate slope estimates
    K1 = np.dot(Amat, x) + np.dot(Bmat, umat)
    K2 = np.dot(Amat, x + K1 * dt / 2) + np.dot(Bmat, umat)
    K3 = np.dot(Amat, x + K2 * dt / 2) + np.dot(Bmat, umat)
    K4 = np.dot(Amat, x + K3 * dt) + np.dot(Bmat, umat)
    
    # Calculate x[k+1] estimate using combination of slope estimates
    x_next = x + 1/6 * (K1 + 2*K2 + 2*K3 + K4) * dt
    print(f'xnext = {xnext}')
    print(f'Shape of xnext = {np.shape(xnext)}')
    
    return x_next

## Create linear model comparison for short time
def get_linear_projection(start_time, duration):
    t_offset = round(start_time/tstep+1)
    print(f't_offset = {t_offset}')
    print(f't_offset Shape = {np.shape(t_offset)}')
    t_steps = round(duration/tstep+1)
    print(f't_steps = {t_steps}')
    print(f't_steps Shape = {np.shape(t_steps)}')
    x_0 = x[:, t_offset]
    print(f'x_0 = {x_0}')
    print(f'shape of x_0 = {np.shape(x_0)}')

    # Linear version of simulation, start with same initial conditions
    x_lin = np.zeros((n_states,t_steps))
    print(f'x_lin = {x_lin}')
    print(f'x_lin shape = {np.shape(x_lin)}')
    x_lin[:,0] = x_0
    print(f'x_lin = {x_lin}')
    print(f'x_lin shape = {np.shape(x_lin)}')

    A_save = np.zeros(t_steps)
    print(f'A_save = {A_save}')
    print(f'Shape of A_save = {np.shape(A_save)}')

    def get_Amat(x0):
        PHI_DOT_0 = x0[3] + (x0[4]*np.sin(x0[6]) + x0[5]*np.cos(x0[6])) * np.sin(x0[7]) / np.cos(x0[7])  # = phidot
        PSI_DOT_0 = (x0[4] * np.sin(x0[6]) + x0[5] * np.cos(x0[6])) / np.cos(x0[7])  # = psidot
        Amat =  jacobian(x0,PHI_DOT_0, PSI_DOT_0, Ixx, Iyy, Izz)
        print(f'Amat = {Amat}')
        print(f'Amat Shape = {np.shape(Amat)}')

        #Solve for nominal Fz
        F1_0 = Fthrust(x0, u[0,t_offset],  dx,  dy)
        F2_0 = Fthrust(x0, u[1,t_offset], -dx, -dy)
        F3_0 = Fthrust(x0, u[2,t_offset],  dx, -dy)
        F4_0 = Fthrust(x0, u[3,t_offset], -dx,  dy)

        return (Amat, F1_0, F2_0, F3_0, F4_0)

    Amat, F1_0, F2_0, F3_0, F4_0 = get_Amat(x_0)

    for k in range(0, np.shape(x_lin)[1] - 1):
        #Amat, F1_0, F2_0, F3_0, F4_0 = get_Amat(x[:,k + t_offset])

        A_save[k] = Amat[0,7]
        print(f'A_save[k] == Amat[0,7] = {A_save[k]}')
        print(f'Shape of A_save[k] = {np.shape(A_save[k])}')
        # Calculate forces from propeller inputs (u)
        F1 = Fthrust(x_lin[:,k], u[0,k + t_offset],  dx,  dy) - F1_0
        F2 = Fthrust(x_lin[:,k], u[1,k + t_offset], -dx, -dy) - F2_0
        F3 = Fthrust(x_lin[:,k], u[2,k + t_offset],  dx, -dy) - F3_0
        F4 = Fthrust(x_lin[:,k], u[3,k + t_offset], -dx,  dy) - F4_0
        Fz = F1 + F2 + F3 + F4
        L = (F1 + F4) * dy - (F2 + F3) * dy
        M = (F1 + F3) * dx - (F2 + F4) * dx
        N = -T(F1,dx,dy) - T(F2,dx,dy) + T(F3,dx,dy) + T(F4,dx,dy)

        umat = np.array([Fz, L, M, N])
        print(f'umat = {umat}')
        print(f'shape of umat = {np.shape(umat)}')

        x_lin_dot = np.dot(Amat, x_lin[:,k]) + np.dot(Bmat, umat)
        print(f'x_lin_dot = {x_lin_dot}')
        print(f'x_lin_dot Shape = {np.shape(x_lin_dot)}')
        print(f'x_lin[:,k] = {x_lin[:,k]}')
        x_lin[:,k+1] = tstep * x_lin_dot + x_lin[:,k]
        print(f'x_lin[:,k+1] = {x_lin[:,k+1]}')
        print(f'x_lin[:,k+1] Shape = {np.shape(x_lin[:,k+1])}')
        #x_lin[:,k+1] = RK4_Linear(x_lin[:,k], umat, Amat, Bmat, tstep)
        
    time_segment = t[t_offset:t_offset + t_steps]
    print(f'time_segment = {time_segment}')
    x_nonlinear = x[:,t_offset:t_offset + t_steps]
    print(f'x_nonlinear = {x_nonlinear}')
    print(f'x_nonlinear Shape = {np.shape(x_nonlinear)}')
    u_segment = u[:,t_offset:t_offset + t_steps] - u[:,t_offset,None]
    print(f'u_segment = {u_segment}')
    print(f'u_segment Shape = {np.shape(u_segment)}')
    return time_segment, x_lin, x_nonlinear, u_segment


##################################################################
################  Plot Results   #################################
##################################################################

plt.figure(1, figsize=(8,8))
plt.subplot(311)
plt.plot(t,x[11,:],'b',label='h')
plt.ylabel('h (m)')
#plt.xlabel('Time (sec)')
#plt.legend(loc='best')
plt.title('Time History of Height, X Position, and Pitch')

plt.subplot(312)
plt.plot(t,x[9,:],'b',label='x')
plt.ylabel('x (m)')
#plt.xlabel('Time (sec)')

plt.subplot(313)
plt.plot(t,x[7,:]*RTD,'b',label='theta')
plt.ylabel('Theta (deg)')
plt.xlabel('Time (sec)')

plt.figure(2, figsize=(8,8))
ax = plt.subplot(1,1,1)
plt.plot(x[9,0:-1:20],x[11,0:-1:20],'bo-',label='y')
plt.text(x[9,0] + 0.1, x[11,0],'START')
plt.text(x[9,-1], x[11,-1],'END')
plt.ylabel('h [m]'); plt.xlabel('x [m]')
ax.axis('equal')
#plt.legend(loc='best')
plt.title('Vertical Profile')

plt.figure(3, figsize=(8,4))
plt.plot(t[0:-1],u[0,0:-1],'b',label='T1')
plt.plot(t[0:-1],u[1,0:-1],'g',label='T2')
plt.plot(t[0:-1],u[2,0:-1],'r',label='T3')
plt.plot(t[0:-1],u[3,0:-1],'y',label='T4')
plt.xlabel('Time (sec)')
plt.ylabel('Propeller RPM')
plt.legend(loc='best')
plt.title('Time History of Control Inputs')

plt.show()


def plot_results(t_array, x_linear_array, x_nonlinear_array, u_array, t, x, u):
    
    plt.figure(figsize=(8,8))
    ax = plt.subplot(1,1,1)
    plt.plot(x[9,0:-1:20],x[11,0:-1:20],'ko-',label='nonlinear')
    plt.text(x[9,0] + 0.1, x[11,0],'START')
    plt.text(x[9,-1], x[11,-1],'END')
    plt.ylabel('h [m]'); plt.xlabel('x [m]')
    ax.axis('equal')
    #plt.legend(loc='best')
    plt.title('Vertical Profile')
    
    color_wheel = 'bgrmcbgrmc'
    for i in range(0,len(x_linear_array)):
        plt.plot(x_linear_array[i][9,0:-1:10],x_linear_array[i][11,0:-1:10],color_wheel[i]+'o-',label='linear')

    plt.figure(figsize=(8,4))
    plt.plot(t[0:-1],u[0,0:-1],'b.-',label='T1')
    plt.plot(t[0:-1],u[1,0:-1],'r.-',label='T2')
    plt.plot(t[0:-1],u[2,0:-1],'g',label='T3')
    plt.plot(t[0:-1],u[3,0:-1],'y',label='T4')
    plt.xlabel('Time (sec)')
    plt.ylabel('Propeller RPM')
    plt.legend(loc='best')
    plt.title('Time History of Control Inputs')

    for i in range(0,len(x_linear_array)):
        plt.figure(figsize=(10,20))
        plt.rcParams["axes.edgecolor"] = "black"
        
        plt.subplot(711)
        plt.plot(t_array[i],u_array[i][0,:],'b.-',label='T1')
        plt.plot(t_array[i],u_array[i][1,:],'g.-',label='T2')
        plt.plot(t_array[i],u_array[i][2,:],'r',label='T3')
        plt.plot(t_array[i],u_array[i][3,:],'y',label='T4')
        plt.xlabel('Time (sec)')
        plt.ylabel('delta Propeller RPM')
        plt.legend(loc='best')
        plt.title('Time History Comparison of Nonlinear and Linear Simulations')
        
        plt.subplot(712)
        plt.plot(t_array[i],x_nonlinear_array[i][0,:],'k',label='nonlinear')
        plt.plot(t_array[i],x_linear_array[i][0,:],color_wheel[i],label='linear')
        plt.legend(loc='best')
        plt.ylabel('u (m/s)')
        plt.xlabel('Time (sec)')
        

        plt.subplot(713)
        plt.plot(t_array[i],x_nonlinear_array[i][4,:],'k',label='nonlinear')
        plt.plot(t_array[i],x_linear_array[i][4,:],color_wheel[i],label='linear')
        plt.legend(loc='best')
        plt.ylabel('q (rad/s)')
        plt.xlabel('Time (sec)')

        plt.subplot(714)
        plt.plot(t_array[i],x_nonlinear_array[i][7,:]*RTD,'k',label='nonlinear')
        plt.plot(t_array[i],x_linear_array[i][7,:]*RTD,color_wheel[i],label='linear')
        plt.legend(loc='best')
        plt.ylabel('Theta (deg)')
        plt.xlabel('Time (sec)')

        plt.subplot(715)
        plt.plot(t_array[i],x_nonlinear_array[i][9,:],'k',label='nonlinear')
        plt.plot(t_array[i],x_linear_array[i][9,:],color_wheel[i],label='linear')
        plt.legend(loc='best')
        plt.ylabel('x (m)')
        plt.xlabel('Time (sec)')

        plt.subplot(716)
        plt.plot(t_array[i],x_nonlinear_array[i][2,:],'k',label='nonlinear')
        plt.plot(t_array[i],x_linear_array[i][2,:],color_wheel[i],label='linear')
        plt.legend(loc='best')
        plt.ylabel('w (m/s)')
        plt.xlabel('Time (sec)')

        plt.subplot(717)
        plt.plot(t_array[i],x_nonlinear_array[i][11,:],'k',label='nonlinear')
        plt.plot(t_array[i],x_linear_array[i][11,:],color_wheel[i],label='linear')
        plt.legend(loc='best')
        plt.ylabel('h(m)')
        plt.xlabel('Time (sec)')
        
    plt.show()


## Compare a selection of starting points projecting linear model 2 seconds
num_conditions = 5
t_array = [None] * num_conditions
x_linear_array = [None] * num_conditions
x_nonlinear_array = [None] * num_conditions
u_array = [None] * num_conditions

duration = 2 #sec

start_time = 9 #sec
t_array[0], x_linear_array[0], x_nonlinear_array[0], u_array[0] = get_linear_projection(start_time, duration)

start_time = 11 #sec
t_array[1], x_linear_array[1], x_nonlinear_array[1], u_array[1] = get_linear_projection(start_time, duration)

start_time = 13 #sec
t_array[2], x_linear_array[2], x_nonlinear_array[2], u_array[2] = get_linear_projection(start_time, duration)

duration = 5 #sec

start_time = 15 #sec
t_array[3], x_linear_array[3], x_nonlinear_array[3], u_array[3] = get_linear_projection(start_time, duration)

duration = 9 #sec

start_time = 20 #sec
t_array[4], x_linear_array[4], x_nonlinear_array[4], u_array[4] = get_linear_projection(start_time, duration)

plot_results(t_array, x_linear_array, x_nonlinear_array, u_array, t, x, u)


## Compare a selection of starting points projecting linear model 2 seconds
duration = 10 #sec

num_conditions = 4
t_array = [None] * num_conditions
x_linear_array = [None] * num_conditions
x_nonlinear_array = [None] * num_conditions
u_array = [None] * num_conditions

start_time = 0 #sec
t_array[0], x_linear_array[0], x_nonlinear_array[0], u_array[0] = get_linear_projection(start_time, duration)

start_time = 0.5 #sec
t_array[1], x_linear_array[1], x_nonlinear_array[1], u_array[1] = get_linear_projection(start_time, duration)

start_time = 1 #sec
t_array[2], x_linear_array[2], x_nonlinear_array[2], u_array[2] = get_linear_projection(start_time, duration)

start_time = 2 #sec
t_array[3], x_linear_array[3], x_nonlinear_array[3], u_array[3] = get_linear_projection(start_time, duration)

plot_results(t_array, x_linear_array, x_nonlinear_array, u_array, t, x, u)


### USING A DIFFERENT SET OF CONTROL INPUTS (starting at 0 RPM) ###
## Compare a selection of starting points projecting linear model 2 seconds
duration = 10 #sec

num_conditions = 4
t_array = [None] * num_conditions
x_linear_array = [None] * num_conditions
x_nonlinear_array = [None] * num_conditions
u_array = [None] * num_conditions

start_time = 0 #sec
t_array[0], x_linear_array[0], x_nonlinear_array[0], u_array[0] = get_linear_projection(start_time, duration)

start_time = 0.5 #sec
t_array[1], x_linear_array[1], x_nonlinear_array[1], u_array[1] = get_linear_projection(start_time, duration)

start_time = 1 #sec
t_array[2], x_linear_array[2], x_nonlinear_array[2], u_array[2] = get_linear_projection(start_time, duration)

start_time = 2 #sec
t_array[3], x_linear_array[3], x_nonlinear_array[3], u_array[3] = get_linear_projection(start_time, duration)

plot_results(t_array, x_linear_array, x_nonlinear_array, u_array, t, x, u)


#Store time, input, and state data to a json file for visualization
import json, codecs
log_filename = "visualizer/sim_log.json"
with open(log_filename, "w") as logfile:
    logfile.write("var sim_data = ")
json.dump([t.tolist(), u.tolist(),x.tolist()], \
          codecs.open(log_filename, 'a', encoding='utf-8'), \
          separators=(',', ':'), indent=4)


