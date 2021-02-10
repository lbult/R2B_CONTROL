from Wind_Reference_System import _All_Dubin_Paths
import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np

#boundary conditions
init_height = 10000
initial_conditions = np.array([50,50,pi/4])
final_conditions = np.array([0,0,0])
tau_f = init_height

#necessary system parameters
V_g = 9.27 #m/s
gamma_g = -25.2 #deg
sigma = pi/6

gamma_min = atan(tan(gamma_g)/cos(sigma))
V_min = sqrt(V_g**2 * cos(gamma_min) / (cos(gamma_g)*cos(sigma)))
r_min = V_min**2 * cos(gamma_min) / (9.81* tan(sigma))
tau_full = abs(2*pi*r_min/(V_min*cos(gamma_min)) * V_min * sin(gamma_min))

print(V_min * sin(gamma_min))

minimum_conditions = _All_Dubin_Paths(pos_init=initial_conditions, 
        pos_final=final_conditions,
        r_traj = r_min, gamma_g_traj = gamma_g,  gamma_traj = gamma_min)

minimum_conditions._Minimum_Tau()
tau_min = minimum_conditions.tau_min

r=r_min

arclengths = minimum_conditions.rsl_traj

print(arclengths)


x_1l = initial_conditions[0] + r*sin(initial_conditions[2]+arclengths[0]) - r*sin(initial_conditions[2])
y_1l = initial_conditions[1] - r*cos(initial_conditions[2]+arclengths[0]) + r*cos(initial_conditions[2])
phi_1l = initial_conditions[2]+arclengths[0]

x_2l = x_1l-arclengths[1]*cos(phi_1l)
y_2l = y_1l+arclengths[1]*sin(phi_1l)
phi_2l = phi_1l

x_3l = x_2l - r*sin(phi_1l-arclengths[2]) + r*sin(phi_1l)
y_3l = y_2l + r*cos(phi_1l-arclengths[2]) - r*cos(phi_1l)
phi_3l = phi_1l - arclengths[2]

update = [initial_conditions[0],x_1l, x_2l, x_3l]

update_2 = [initial_conditions[1],y_1l, y_2l, y_3l]
plt.plot(update,update_2, 'r')

plt.show()