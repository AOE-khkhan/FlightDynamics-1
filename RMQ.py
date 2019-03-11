from numpy import zeros

def RMQ(q1,q2,q3,q4):
    '''
    FLIGHT Earth-to-Body-Axis Rotation (Direction-Cosine) Matrix from Quaternion
    #   Quaternion:
            q1, x Component of quaternion
            q2, y Component of quaternion
            q3, z Component of quaternion
            q4, cos(Euler) Component of quaternion   
    '''

    H = zeros((3,3))

    H[0,0] = q1**2 - q2**2 - q3**2 + q4**2
    H[0,1] = 2*(q1*q2 + q3*q4)
    H[0,2] = 2*(q1*q3 - q2*q4)
    H[1,0] = 2*(q1*q2 - q3*q4)
    H[1,1] = -q1**2 + q2**2 - q3**2 + q4**2
    H[1,2] = 2*(q2*q3 + q1*q4)
    H[2,0] = 2*(q1*q3 + q2*q4)
    H[2,1] = 2*(q2*q3 - q1*q4)
    H[2,2] = -q1**2 - q2**2 + q3**2 + q4**2    

    return H
