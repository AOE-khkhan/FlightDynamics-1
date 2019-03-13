from numpy import zeros, sin, cos


def DCM(Phi,Theta,Psi):
        '''
                FLIGHT Earth-to-Body-Axis Direction-Cosine Matrix

                # Euler Angles:
                        # phi,    Roll Angle,  rad
                        # Theta,  Pitch Angle, rad
                        # Psi,    Yaw Angle,   rad
        '''
        sinR = sin(Phi)
        cosR = cos(Phi)
        sinP = sin(Theta)
        cosP = cos(Theta)
        sinY = sin(Psi)
        cosY = cos(Psi)

        H = zeros((3,3))

        H[0,0] = cosP*cosY
        H[0,1] = cosP*sinY
        H[0,2] = -sinP
        H[1,0] = sinR*sinP*cosY - cosR*sinY
        H[1,1] = sinR*sinP*sinY + cosR*cosY
        H[1,2] = sinR*cosP
        H[2,0] = cosR*sinP*cosY + sinR*sinY
        H[2,1] = cosR*sinP*sinY - sinR*cosY
        H[2,2] = cosR*cosP

        return H

if __name__ == "__main__":
        DCM(Phi,Theta,Psi)
        