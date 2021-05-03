from Dubin_Path import _All_Dubin_Paths
import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos, sqrt
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
        pos_final=final_conditions, gamma_g_traj = gamma_g, altitude=init_height, v_g=V_g, sigma_max=sigma,
        wind_direction=45, wind_magnitude=0.4, monte_carlo=1000)

minimum_conditions._Minimum_Tau()

plt.plot(minimum_conditions.pos_xs, minimum_conditions.pos_ys, 'b')
plt.plot(minimum_conditions.pos_x, minimum_conditions.pos_y, 'r')

## plotting the first path that monte carlo feeds back, i made it so that this first path would just be the trajectory with no noise in the wind
plt.plot(minimum_conditions.pos_x_w[0], minimum_conditions.pos_y_w[0], 'g', alpha=1)

for i in range(len(minimum_conditions.pos_x_w)):
        ## plotting the distribution of landing coordinates
        # plt.plot(minimum_conditions.pos_x_w[i], minimum_conditions.pos_y_w[i], 'g', alpha=0.01)
        plt.scatter(minimum_conditions.pos_x_w[i][-1], minimum_conditions.pos_y_w[i][-1], c='r', alpha=0.05)

# print(minimum_conditions.pos_x_w, minimum_conditions.pos_y_w)
# print(minimum_conditions.pos_y, "<--pos_x")
# print(minimum_conditions.pos_y_w, "<--pos_x_w")
# print(minimum_conditions.heading)
# print(minimum_conditions.alt)

plt.gca().set_aspect('equal', adjustable='box')

plt.show()