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
        pos_final=final_conditions, gamma_g_traj = gamma_g, altitude=init_height, v_g=V_g, sigma_max=sigma,
        wind_direction=12, wind_magnitude=0.4, monte_carlo=100)

minimum_conditions._Minimum_Tau()

plt.plot(minimum_conditions.pos_xs, minimum_conditions.pos_ys, 'b')
plt.plot(minimum_conditions.pos_x, minimum_conditions.pos_y, 'r')
plt.plot(minimum_conditions.pos_x_w[0], minimum_conditions.pos_y_w[0], 'g', alpha=1)
# plt.plot(minimum_conditions.pos_x_ws[0], minimum_conditions.pos_y_ws[0], 'g', alpha=1)

for i in range(len(minimum_conditions.pos_x_w)):
        ## plotting the land_fall coordinates of the parafoil after monte carlo
        plt.scatter(minimum_conditions.pos_x_w[i][-1], minimum_conditions.pos_y_w[i][-1], c='r', alpha=0.01)
        # plt.scatter(minimum_conditions.pos_x_ws[i][-1], minimum_conditions.pos_y_ws[i][-1], c='r', alpha=0.3)

# print(minimum_conditions.pos_x_w, minimum_conditions.pos_y_w)
# print(minimum_conditions.pos_y, "<--pos_x")
# print(minimum_conditions.pos_y_w, "<--pos_x_w")
# print(minimum_conditions.heading)
# print(minimum_conditions.alt)

plt.gca().set_aspect('equal', adjustable='box')

plt.show()