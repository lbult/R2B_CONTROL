from Wind_Reference_System import _All_Dubin_Paths
import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np

#boundary conditions
init_height = 10000
initial_conditions = np.array([0,0,0])
final_conditions = np.array([-30,-40,pi])
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
minimum_conditions._RSL()

print(minimum_conditions.tau_lsl, minimum_conditions.lsl_traj)

#minimum_conditions._Minimum_Tau()
#tau_min = minimum_conditions.tau_min

r=r_min

arclengths = minimum_conditions.lsl_traj
arclengths_1 = minimum_conditions.rsl_traj

print(arclengths)

if arclengths[0] < 0:
        arclengths += np.array([2*pi,0,0])

print(arclengths_1)

#initial_conditions = initial_conditions + np.array([0,0,-pi/2])
'''
#x_1l = initial_conditions[0] + r*sin(initial_conditions[2]+arclengths[0]) - r*sin(initial_conditions[2])
x_1l = -50 + r*sin(arclengths[0]) - r*sin(initial_conditions[2])
#y_1l = initial_conditions[1] + r*cos(initial_conditions[2]+arclengths[0]) - r*cos(initial_conditions[2])
y_1l = 40 - r*cos(arclengths[0]) + r*cos(initial_conditions[2])
phi_1l = arclengths[0]

x_2l = x_1l+arclengths[1]*cos(phi_1l)
y_2l = y_1l+arclengths[1]*sin(phi_1l)
phi_2l = phi_1l'''

x=[]
y=[]

f=[]
g=[]

def draw_circle(theta, theta_offset, x_1, y_1):
        dtheta = theta / 100
        f = 0
        theta_m = theta_offset
        while f < 100:
                f += 1
                x.append(x_1+r*cos(theta_m))
                y.append(y_1+r*sin(theta_m))
                theta_m += dtheta
        

draw_circle(-arclengths_1[0], pi, r, 0)
phi_1l = pi/2-arclengths_1[0]

f=x
g=y

x_2l = f[99]+arclengths[1]*cos(phi_1l)
y_2l = g[99]+arclengths[1]*sin(phi_1l)

draw_circle(-arclengths_1[2]-pi/2, phi_1l+pi/2, x_2l-r*sin(phi_1l), y_2l-r*cos(phi_1l))

plt.gca().set_aspect('equal', adjustable='box')
plt.plot(x,y,'r')
plt.plot(f,g,'r')

plt.show()
'''
For plotting LSL
draw_circle(arclengths[0], -pi/2, 0, r)
f=x
g=y

x_2l = f[99]+arclengths[1]*cos(phi_1l)
y_2l = g[99]+arclengths[1]*sin(phi_1l)

draw_circle(arclengths[2], -pi/2+arclengths[0], x_2l-r*sin(arclengths[0]), y_2l+r*cos(arclengths[0]))

plt.gca().set_aspect('equal', adjustable='box')
plt.plot(x,y,'r')
plt.plot(f,g,'r')

plt.show()
'''