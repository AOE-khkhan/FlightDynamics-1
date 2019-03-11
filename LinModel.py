from numpy import array, concatenate

# globals ?????

def LinModel(tj,xj):
    '''
        Flight Equations of Motion for Linear Model (Jacobian) Evaluation,
        with dummy state elements added for controls
    '''

    x   =   xj[0:11+1]
    u   =   xj[12:18+1]

    xdot    =   EoM(tj,x)
    xdotj   =   concatenate(([xdot], np.array([0;0;0;0;0;0;0])))

    return xdotj
    