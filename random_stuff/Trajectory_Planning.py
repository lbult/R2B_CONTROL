from Dubin_Path import _All_Dubin_Paths
import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np

#boundary conditions
init_height = 500
initial_conditions = np.array([0,0,0])
final_conditions = np.array([100,100,pi/2])
tau_f = init_height

#necessary system parameters
V_g = 9.27 #m/s
gamma_g = -25.2 * pi /180 #deg
sigma = pi/6

minimum_conditions = _All_Dubin_Paths(pos_init=initial_conditions, 
        pos_final=final_conditions, gamma_g_traj = gamma_g, altitude=init_height, v_g=V_g, sigma_max=sigma)

minimum_conditions._Minimum_Tau()

plt.plot(minimum_conditions.pos_xs, minimum_conditions.pos_ys, 'b')
plt.plot(minimum_conditions.pos_x, minimum_conditions.pos_y, 'r')

#print(minimum_conditions.pos_x)
#print(minimum_conditions.heading)
print(minimum_conditions.alt)

plt.gca().set_aspect('equal', adjustable='box')

plt.show()