import numpy as np
#
mud = 0.1
mus = 0.677
stsxx_max = -135e6
stsyy_max = -120e6
stsxy_max = 70e6
nu = 0.25
G = 32.04e9
Dc = 0.4
#
stszz_max = 0
mu = np.abs( stsxy_max / stsyy_max )
S = ( mus - mu ) / ( mu - mud )
#
shear_stress = np.sqrt( (stsxy_max/2)**2 )
shear_strength = abs( mus * (stsyy_max/2) )
#
L = G * Dc / ( ( mus - mud ) * abs(stsyy_max/2) )
#
print("stszz_max: ",stszz_max)
print("mu: ",mu)
print("S: ",S)
print("shear_stress: ",shear_stress)
print("shear_strength: ",shear_strength)
print("L: ",L)