import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np
import scipy

def _Gen_Wind_Field(x_speed, y_speed):
    #try later with periodical disturbances at time tau
    return np.array([x_speed, y_speed, 0])


def _Return_Minimum_Value(values):
    return min(values)

def _RSL(lambdas, d, r_traj, gamma_traj, gamma_g_traj):
    p_rsl = (-2*r_traj**2 + d**2 + 2*r_traj**2*cos(lambdas) + 2*r_traj*d*(sin(lambdas)))**0.5
    q_rsl = -atan(2*r_traj/p_rsl) + atan((-r_traj*cos(lambdas)-r_traj)/(d+r_traj*sin(lambdas)))
    t_rsl = - lambdas - atan(2*r_traj/p_rsl) + atan((-r_traj*cos(lambdas)-r_traj)/(d+r_traj*sin(lambdas)))
    tau_rsl = r_traj*tan(gamma_traj*(lambdas+2*t_rsl)) + p_rsl*tan(gamma_g_traj)
    return tau_rsl, np.array([t_rsl/r_traj, p_rsl/r_traj, q_rsl/r_traj])

#necessary system parameters
V_g = 9.27 #m/s
gamma_g = -25.2 #deg
psi_dot_max = 30 #deg/s

#boundary conditions
init_height = 10
final_conditions = np.array([0,0,0])
wind = _Gen_Wind_Field(2,2)
tau_f = init_height

#initial values of time changing variables
gamma = 0
sigma = 1
V = 0
#tau = init_height
x_dot = 0
y_dot = 0
phi_dot = 0
init_conditions = np.array([
    -20, #-abs(tau_f*wind[0]/(V_g*sin(gamma_g))) - 10,
    -20, #abs(tau_f*wind[1]/(V_g*sin(gamma_g))) - 10,
    0
    ])

#assume continuous density over h
density = 1.225 #kg/m^3

dt = 0.01
tau_min = 0.01
sim = True

while sim:
    
    gamma = atan(tan(gamma_g)/cos(sigma))
    V = sqrt(V_g**2 * cos(gamma) / (cos(gamma_g)*cos(sigma)))

    r = V**2 * cos(gamma) / (9.81* tan(sigma))
    
    tau_min = _RSL(init_conditions[2], 
    sqrt(abs(init_conditions[1])**2 + abs(init_conditions[0])**2),
    r, gamma, gamma_g
    )[0]

    arclengths = _RSL(init_conditions[2], 
    sqrt(abs(init_conditions[1])**2 + abs(init_conditions[0])**2),
    r, gamma, gamma_g
    )[1]

    #parameters that change with time
    if abs(tau_min) > 1.01*tau_f:
        sigma += 0.001
    elif abs(tau_min) < 0.99*tau_f:
        sigma -= 0.001
    else:
        sim = False 
        print(sigma)
        print(r)

    
    #tau += -V*sin(gamma) * dt
    
    '''a_tilda = -1/(tan(gamma))
    u = 1/(r*cos(gamma))

    init_conditions = init_conditions + np.array([
        a_tilda * cos(gamma),
        a_tilda * sin(gamma),
        a_tilda * u
    ])'''

    '''
    init_conditions = init_conditions + np.array([
        -(cos(init_conditions[2])/tan(gamma) + wind[0]/V*sin(gamma)),
        -(sin(init_conditions[2])/tan(gamma) + wind[1]/V*sin(gamma)),
        -9.81*tan(sigma)/(V**2*sin(gamma))
    ])

    phi_dot = 9.81 * tan(sigma) / V'''

total_length = V * cos(gamma) * init_height #/(V*sin(gamma))

arclengths = arclengths*init_height

print(arclengths)

x_1l = []
y_1l = []

x_2l = []
y_2l = []

x_3l = []
y_3l = []

dtheta_r = -arclengths[0]
theta_r = pi + dtheta_r
dtheta_l = -arclengths[2]
theta_l= -pi/2 + dtheta_l

print(dtheta_l)
print(dtheta_r)

def func_1():
    my_one = 0
    interval_one = -pi/500#-arclengths[0]/(1000*r)
    my_theta_r = pi
    while my_one < 1000:
        x_1l.append(init_conditions[0] + r + r*cos(my_theta_r))
        y_1l.append(init_conditions[1] + r*sin(my_theta_r))
        my_theta_r += interval_one
        my_one += 1

def func_2():
    my_two = 0
    interval_two = -pi/500#-arclengths[2]/(1000*r)
    my_theta_l = -pi/2
    while my_two < 1000:
        x_2l.append(r*cos(my_theta_l ))
        y_2l.append(r + r*sin(my_theta_l))
        my_theta_l += interval_two
        my_two += 1

def func_3():
    dy = y_2l[-1] - y_1l[-1]
    dx = x_2l[-1] - x_1l[-1]
    a = dy/dx
    b = y_1l[-1] - dy/dx * x_1l[-1]

    x_33 = x_1l[-1]

    my_three = 0
    interval_three = dx/1000
    while my_three < 1000:
        y_3l.append(+dy/dx*x_33+b)
        x_3l.append(x_33)
        x_33 += interval_three
        my_three += 1

func_1()
func_2()
func_3()

plt.plot(x_1l,y_1l, 'r')
plt.plot(x_2l,y_2l, 'r')
plt.plot(x_3l,y_3l, 'r')

# show the plot
plt.show()