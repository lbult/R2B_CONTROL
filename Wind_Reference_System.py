import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np
import scipy

def _Gen_Wind_Field(x_speed, y_speed):
    #try later with periodical disturbances at time tau
    return np.array([x_speed, y_speed, 0])


def _Return_Minimum_Value(values):
    return min(values)

def _RSL(phi_init, d, r_traj, gamma_traj, gamma_g_traj, gammas_input):
    gammas = atan(gammas_input)
    lambdas = gammas-phi_init
    mus = gammas
    
    p_rsl = (abs(-2*r_traj**2 + d**2 + cos(lambdas-mus)*2*r_traj**2 + 2*r_traj*d*(sin(lambdas)+ sin(mus))))**0.5
    q_rsl = - mus - atan(2*r_traj/p_rsl) + atan((-r_traj*cos(lambdas)-r_traj*cos(mus))/(d+r_traj*sin(lambdas)+r_traj*sin(mus)))
    t_rsl = - lambdas - atan(2*r_traj/p_rsl) + atan((-r_traj*cos(lambdas)-r_traj*cos(mus))/(d+r_traj*sin(lambdas)+r_traj*sin(mus)))
    
    tau_rsl = r_traj*tan(gamma_traj)*(lambdas-mus+2*t_rsl) + p_rsl*tan(gamma_g_traj)
    
    return tau_rsl, np.array([t_rsl/r_traj, p_rsl/r_traj, q_rsl/r_traj])

#necessary system parameters
V_g = 9.27 #m/s
gamma_g = -25.2 #deg
psi_dot_max = 30 #deg/s

#boundary conditions
init_height = 70
final_conditions = np.array([0,0,0])
wind = _Gen_Wind_Field(2,2)
tau_f = init_height

#initial values of time changing variables
gamma = 0
sigma = pi/6
V = 0
#tau = init_height
x_dot = 0
y_dot = 0
phi_dot = 0
init_conditions = np.array([
    -10, #-abs(tau_f*wind[0]/(V_g*sin(gamma_g))) - 10,
    -10, #abs(tau_f*wind[1]/(V_g*sin(gamma_g))) - 10,
    0
    ])

#assume continuous density over h
density = 1.225 #kg/m^3

dt = 0.01
tau_min = 0.01
sim = True

gamma_min = atan(tan(gamma_g)/cos(pi/6))
V_min = sqrt(V_g**2 * cos(gamma_min) / (cos(gamma_g)*cos(pi/6)))
r_min = V_min**2 * cos(gamma_min) / (9.81* tan(pi/6))
tau_full = abs(2*pi*r_min/(V_min*cos(gamma_min)) * V_min * sin(gamma_min))

tau_min = _RSL(init_conditions[2], 
    sqrt((init_conditions[1])**2 + (init_conditions[0])**2),
    r_min, gamma_min, gamma_g, init_conditions[1]/init_conditions[0]
    )[0]

print(tau_min)
print(tau_full)
print(tau_f)

eta = (tau_f - tau_min)/tau_full

print(eta)
print(eta)

if eta > 0 and tau_min>0:
    while sim:
        
        gamma = atan(tan(gamma_g)/cos(sigma))
        V = sqrt(V_g**2 * cos(gamma) / (cos(gamma_g)*cos(sigma)))

        r = V**2 * cos(gamma) / (9.81* tan(sigma))
        
        tau = _RSL(init_conditions[2], 
        sqrt(abs(init_conditions[1])**2 + abs(init_conditions[0])**2),
        r, gamma, gamma_g, init_conditions[1]/init_conditions[0]
        )[0]

        arclengths = _RSL(init_conditions[2], 
        sqrt(abs(init_conditions[1])**2 + abs(init_conditions[0])**2),
        r, gamma, gamma_g, init_conditions[1]/init_conditions[0]
        )[1]

        #parameters that change with time
        if tau > 1.01*tau_f:
            sigma += 0.001
        elif tau < 0.99*tau_f:
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

    #total_length = V * cos(gamma) * init_height #/(V*sin(gamma))

    arclengths = arclengths

    print(arclengths)


    x_1l = init_conditions[0] - r*sin(init_conditions[2]+arclengths[0]) + r*sin(init_conditions[2])
    y_1l = init_conditions[1] + r*cos(init_conditions[2]+arclengths[0]) - r*cos(init_conditions[2])
    phi_1l = init_conditions[2]+arclengths[0]

    x_2l = x_1l+arclengths[1]*cos(phi_1l)
    y_2l = y_1l+arclengths[1]*sin(phi_1l)
    phi_2l = phi_1l

    x_3l = x_2l + r*sin(init_conditions[2]-arclengths[2]) - r*sin(init_conditions[2])
    y_3l = y_2l - r*cos(init_conditions[2]-arclengths[2]) + r*cos(init_conditions[2])
    phi_3l = phi_1l - arclengths[2]

    update = [init_conditions[0],x_1l, x_2l, x_3l]

    update_2 = [init_conditions[1],y_1l, y_2l, y_3l]

    plt.plot(init_conditions[0], init_conditions[1], 'bo')
    plt.plot(update,update_2, 'r')
    
    #plt.plot(x_2l,y_2l, 'bo')
    #plt.plot(x_3l,y_3l, 'bo')

    # show the plot
    plt.show()