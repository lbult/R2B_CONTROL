import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np
import scipy

def _Gen_Wind_Field(x_speed, y_speed):
    #try later with periodical disturbances at time tau
    return np.array([x_speed, y_speed, 0])

#necessary system parameters
V_g = 9.27 #m/s
gamma_g = -25.2 #deg
psi_dot_max = 30 #deg/s

#boundary conditions
init_height = 50
final_conditions = np.array([0,0,0])
wind = _Gen_Wind_Field(2,2)
tau_f = init_height

#initial values of time changing variables
gamma = 0
V = 0
tau = init_height
x_dot = 0
y_dot = 0
phi_dot = 0
init_conditions = np.array([
    -tau_f*wind[0]/(V_g*sin(gamma_g)),
    -tau_f*wind[1]/(V_g*sin(gamma_g)),
    0
    ])

#assume continuous density over h
density = 1.225 #kg/m^3

sim_iterations = 0
dt = 0.01

while sim_iterations < 0:
    #parameters that change with time
    sigma = 0 #command changing turn rate
    sigma_com = 0 #commanded turn rate
    gamma = atan(tan(gamma_g)/cos(sigma))
    V = sqrt(V_g**2 * cos(gamma) / (cos(gamma_g)*cos(sigma)))

    tau += -V*sin(gamma) * dt

    r = V**2 * cos(gamma) / (9.81* tan(sigma))
    a_tilda = -1/(tan(gamma))
    u = 1/(r*cos(gamma))

    init_conditions = init_conditions + np.array([
        a_tilda * cos(gamma),
        a_tilda * sin(gamma),
        a_tilda * u
    ])

    '''
    init_conditions = init_conditions + np.array([
        -(cos(init_conditions[2])/tan(gamma) + wind[0]/V*sin(gamma)),
        -(sin(init_conditions[2])/tan(gamma) + wind[1]/V*sin(gamma)),
        -9.81*tan(sigma)/(V**2*sin(gamma))
    ])

    phi_dot = 9.81 * tan(sigma) / V'''