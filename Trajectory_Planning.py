from Wind_Reference_System import _All_Dubin_Paths
import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np

#boundary conditions
init_height = 10000
initial_conditions = np.array([0,0,0])
final_conditions = np.array([1,40,3*pi/2])
tau_f = init_height

#necessary system parameters
V_g = 9.27 #m/s
gamma_g = -25.2 * pi /180 #deg
sigma = pi/6

gamma_min = np.arctan2(tan(gamma_g), cos(sigma))
V_min = sqrt(V_g**2 * cos(gamma_min) / (cos(gamma_g)*cos(sigma)))
r_min = (V_min)**2 * cos(gamma_min) / (9.81* tan(sigma))
tau_full = abs(2*pi*r_min/(V_min*cos(gamma_min)) * V_min * sin(gamma_min))

print(r_min)

minimum_conditions = _All_Dubin_Paths(pos_init=initial_conditions, 
        pos_final=final_conditions,
        r_traj =r_min, gamma_g_traj = gamma_g,  gamma_traj = gamma_min)

minimum_conditions._LSL()
minimum_conditions._RSR()

r=r_min

arclengths = minimum_conditions.lsl_traj
arclengths_1 = minimum_conditions.rsr_traj

print(arclengths_1)


if arclengths[0] < 0:
        arclengths += np.array([2*pi,0,0])

minimum_conditions._Go_Left(arclengths[0])
minimum_conditions._Straight(arclengths[1])
minimum_conditions._Go_Left(arclengths[2])
plt.plot(minimum_conditions.pos_x,minimum_conditions.pos_y,'r')
minimum_conditions._Remove_Path()

minimum_conditions._Go_Right(arclengths_1[0])
minimum_conditions._Straight(arclengths_1[1])
minimum_conditions._Go_Right(arclengths_1[2])
plt.plot(minimum_conditions.pos_x,minimum_conditions.pos_y,'y')
minimum_conditions._Remove_Path()

final_conditions = np.array([-40,1,pi])

minimu_conditions = _All_Dubin_Paths(pos_init=initial_conditions, 
        pos_final=final_conditions,
        r_traj =r_min, gamma_g_traj = gamma_g,  gamma_traj = gamma_min)

minimu_conditions._LSR()
minimu_conditions._RSL()

arclengths = minimu_conditions.lsr_traj
arclengths_1 = minimu_conditions.rsl_traj

print(arclengths_1)

minimu_conditions._Go_Right(arclengths_1[0])
minimu_conditions._Straight(arclengths_1[1])
minimu_conditions._Go_Left(arclengths_1[2])
plt.plot(minimu_conditions.pos_x,minimu_conditions.pos_y,'b')
minimu_conditions._Remove_Path()

minimu_conditions._Go_Left(arclengths[0])
minimu_conditions._Straight(arclengths[1])
minimu_conditions._Go_Right(arclengths[2])
plt.plot(minimu_conditions.pos_x,minimu_conditions.pos_y,'g')
minimu_conditions._Remove_Path()

plt.gca().set_aspect('equal', adjustable='box')

plt.show()
