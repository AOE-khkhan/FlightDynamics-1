from numpy import array, reshape, cos, sin, matmul, concatenate, hstack

# globals ?????

def TrimCost(OptParam):
    '''
        FLIGHT Cost Function for Longitudinal Trim in Steady Level Flight
    '''
    R = array([[1,0,0],[0,1,0],[0,0,1]])
    R = R.reshape((3,3))

    # Optimizing Vector:
    #   1 = Stabilator, rad
    #   2 = Throttle, #
    #   3 = Pitch Angle, rad
    OptParam

    u = [u[0]
        u[1]
        u[2]
        OptParam[1]
        u[4]
        u[5]
        OptParam[0]]

    x	=	[V * cos(OptParam[2])
            x[1]
            V * sin(OptParam[2])
            x[3]
            x[4]
            x[5]
            x[6]
            x[7]
            x[8]
            x[9]
            OptParam[2]
            x[11]]
    
    xdot    = EoM(1,x)
    xCost   = array.([xdot[0], xdot[2], xdot[7]]).reshape((3,1))
    # J = quadretic cost function    # otehr choices avilable to compute cost function(J)
    J       =   matmul(matmul(xCost.T, R), xCost)    ## xCost' * R * xCost
    ParamCost   =   concatenate((OptParam, J))   ##[OptParam;J];   # column vector
    TrimHist    =   hstack((TrimHist, ParamCost))   # 4 rows showing history of trim computation over "INDEX" times 
    ##set TrimHist = [0,0,0,0]

    
    return J

