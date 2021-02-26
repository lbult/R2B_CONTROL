import numpy as np
from math import cos, sin
import matplotlib.pyplot as plt

def PID(kp, kd, reference_trajectory, actual_position, sigma_ref, velocity):
    '''
    actual position: nd_array(3) --> x,y,heading
    kp --> proportian gain
    kd --> derivative gain
    sigma_ref --> control input suggested by the trajectory planning
    reference_trajectory --> trajectory generated
    '''

    #initiate errors in position and heading
    e1 = actual_position[0] - reference_trajectory[0]
    e2 = actual_position[1] - reference_trajectory[1]
    e3 = actual_heading - reference_trajectory[2]

    #define proportional and derivative
    e_perp = e1*cos(reference_trajectory[2]) + e2*sin(reference_trajectory[2])
    #e_dot_perp =

    #define control input
    sigma_com = sigma_ref + kp*e_perp #+ kd*e_dot_perp
    return sigma_com


hover_altitude = 100
sc_altitude = 0
sc_velocity = 0

kp = 6
kd = 0.05

dt = 0.1

altitude_save = []
error_save = [100]
error_plot = []
time = []
t = 0

run = 0

while run < 500:
    e = hover_altitude-sc_altitude

    e_dot = (error_save[-1]-e)/dt

    error_save.append(e)
    error_plot.append(abs(e))

    u_t = kp*e + kd*sc_velocity - 9.81

    altitude_save.append(sc_altitude)

    sc_velocity += u_t*dt
    sc_altitude += sc_velocity*dt + 0.5 * dt**2 * u_t
    t+=dt
    time.append(t)
    run+=1

plt.plot(time, error_plot)
plt.show()