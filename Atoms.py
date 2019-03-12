# import
import numpy as np

def Atoms(geomAlt):
    '''
    Note:   Function does not extrapolate outside altitude range
    Input:  Geometric Altitude, m (positive up)
    Output: Air Density, kg/m^3
            Air Pressure, N/m^2
            Air Temprature, K
            Speed of Sound, m/s
    '''

    # Values Tabualated by Geometric Altitude
    Z       = np.array([-1000,0,2500,5000,10000,11100,15000,20000,47400,51000])  # geometric altitude, m
    H       = np.array([-1000,0,2499,4996,9984,11081,14965,19937,47049,50594])   # geopotential altitude, m
    ppo     = np.array([1,1,0.737,0.533,0.262,0.221,0.12,0.055,0.0011,0.0007])   # p/p0*  #p=standard pressure  #p0*=1.013*10^5 N/m^2
    rro     = np.array([1,1,0.781,0.601,0.338,0.293,0.159,0.073,0.0011,0.0007])  # rho/rho0*  #r=standard density  #r0*=1.225 kg/m^3
    T	    = np.array([288.15,288.15,271.906,255.676,223.252,216.65,216.65,216.65,270.65,270.65])  # Tempraure, K
    a	    = np.array([340.294,340.294,330.563,320.545,299.532,295.069,295.069,295.069,329.799,329.799])  # speed of sound, m/s
    R		= 6367435	# Mean radius of the earth, m
    Dens	= 1.225	    # Air density at sea level, Kg/m^3
    Pres	= 101300	# Air pressure at sea level, N/m^2    

    # Geopotential Altitude, m
    geopAlt = R * geomAlt / (R + geomAlt)

    # Linear Interpolation in Geopotential Altitude for Temprature and Speed of Sound
    temp        = np.interp(geopAlt, Z,T)
    soundSpeed  = np.interp(geopAlt, Z,a)

    # Exponential Interpolation in Geometric Altitude for Air Density and Pressure
    for k in range(2, 10+1):
        if geomAlt <= Z[k]:
            betap	= np.log(ppo[k] / ppo[k-1]) / (Z[k] - Z[k-1])
            betar   = np.log(rro[k] / rro[k-1]) / (Z[k] - Z[k-1])
            airPres = Pres * ppo[k-1] * np.exp(betap * (geomAlt - Z[k-1]))
            airDens = Dens * rro[k-1] * np.exp(betar * (geomAlt - Z[k-1]))
            break
    
    return (airDens, airPres, temp, soundSpeed)

if __name__ == "__main__":
    Atoms(geomAlt)