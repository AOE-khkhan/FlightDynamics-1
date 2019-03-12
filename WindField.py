from numpy import array, interp, matmul, reshape, ravel
from DCM import DCM


def WindField(height,phir,thetar,psir):
    '''
        Flight Wind Field Interpolation for 3-D wind as a Function of Altitude
    '''

    windh = array([-10, 0, 100, 200, 500, 1000, 2000, 4000, 8000, 16000])	# Height, m
    windx = array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])	# Northerly wind, m/s
    windy =	array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])	# Easterly wind, m/s
    windz =	array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])	# Vertical wind. m/s
    
    winde = [interp(height, windh,windx), interp(height, windh,windy), interp(height, windh,windz)];   # Earth-relative frame
    winde = reshape(winde, (3,1))
    HEB     =   DCM(phir,thetar,psir)
    windb   =   matmul(HEB, winde).ravel()                # Body-axis frame  (3x1) matrix

    return windb

if __name__ == "__main__":
    WindField(height,phir,thetar,psir)
