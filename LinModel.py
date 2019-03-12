from numpy import array, concatenate
#from globals import GEAR, CONHIS, SPOIL, u, x, V, uInc, tuHis, deluHis, TrimHist, SMI, MODEL, RUNNING


# globals ?????

def LinModel(tj,xj):
    '''
        Flight Equations of Motion for Linear Model (Jacobian) Evaluation,
        with dummy state elements added for controls
    '''

    global u

    x   =   xj[0:11+1]
    u   =   xj[12:18+1]

    xdot    =   EoM(tj,x)
    xdotj   =   concatenate(([xdot], np.array([0;0;0;0;0;0;0])))

    return xdotj

if __name__ == "__main__":
    LinModel(tj,xj)