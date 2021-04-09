import math
import numpy as np
from math import pi, sqrt, asin
from matplotlib import pyplot as plt

#import different support files
from SelfFly import ParafoilProperties, Dynamics, Payload
from support import Quaternion, Variable, takeClosest
from Dubin_Path import _All_Dubin_Paths
from Control_Filter import _Controller


ts = 0.02
#define parafoil-payload system properties
parafoil = ParafoilProperties(m=635, alpha_0=(-3.5*pi/180), surface=697, 
    a_0=6.41, Cd0=0.005,rigging=10*pi/180, ts=ts)
mpayload = Payload(M=9979,payload_cd=0.02)

#start dynamics and attitude tracking
parafoil_dynamics = Dynamics(dt=ts, mass=(parafoil.m+mpayload.M))
parafoil_attitude = Quaternion()
parafoil_attitude.psi = 0
parafoil_attitude._to_quaternion()

#intialize variables to be tracked
alpa = Variable("alpha")
pos_x = Variable("Position [x]")
pos_y = Variable("Position [y]")
pos_z = Variable("Position [z]")
vel_x = Variable("Velocity [x]")
vel_y = Variable("Velocity [y]")
vel_z = Variable("Velocity [z]")

#Noisy Data - Gaussian Noise with sensor standard deviation
pos_x_noise = Variable("Noisy Position [x]")
pos_y_noise = Variable("Noisy Position [y]")
pos_z_noise = Variable("Noisy Position [z]")
vel_x_noise = Variable("Noisy Velocity [x]")
vel_y_noise = Variable("Noisy Velocity [y]")
vel_z_noise = Variable("Noisy Velocity [z]")

#Control code variables:
error_time = Variable("Error")
error_time.update_history(0)
Desired_heading = Variable("Desired heading")
Desired_heading.update_history(0)
TE_deflection = Variable(var_name="TE Deflection", limit=math.pi / 2 )
Control_Input = Variable(var_name="Control_Input", limit=0.5)

print(parafoil._Calc_Alpha_Trim(0.02, 7*pi/180))

current_alt = 0

start = True
controls = False
calc_dubin = True
first_turn = True

while start or pos_z.history[-1] > 0:
    
    #after initialisation of onboard navigation
    if controls:
        index = minimum_conditions.alt.index(takeClosest(minimum_conditions.alt, 
            parafoil_dynamics.pos[2]*2))
        control_input = minimum_conditions.control[index]
        
        _Controller(minimum_conditions.heading[index], 
        parafoil_dynamics.pos, minimum_conditions.pos_x[index],  
        minimum_conditions.pos_y[index], parafoil_attitude.psi,
        error_time, Control_Input, Desired_heading,
        parafoil_dynamics)

        Desired_heading.update_history(minimum_conditions.heading[index])

        if control_input != 0:
            TE_deflection.update_history(-control_input*TE + Control_Input.history[-1])
            parafoil.TE_DEF = TE_deflection.history[-1]
            parafoil.bank = minimum_conditions.sigma_max
        else:
            TE_deflection.update_history(Control_Input.history[-1])
            parafoil.TE_DEF = TE_deflection.history[-1]

    #update quaternion and gravity matrix
    parafoil_attitude.omega = parafoil._Parafoil_Control(parafoil_dynamics.turn_vel)
    parafoil_attitude._to_quaternion()
    parafoil_attitude._update_quaternion()
    parafoil_attitude._to_euler()
    
    #update parafoil forces and moments
    Payload_Vector = mpayload._Calc_Forces(parafoil_dynamics.vel)
    Parafoil_Vector = parafoil._Parafoil_Forces_Moments(parafoil_dynamics.vel)
    Gravity_Vector = parafoil_attitude.body_g * 9.80665 * parafoil_dynamics.mass
    #total force vector
    parafoil_dynamics.forces = Parafoil_Vector + Payload_Vector + Gravity_Vector

    # input them into Dynamics
    parafoil_dynamics.update_dynamics(parafoil_attitude._rot_b_v())
    parafoil_dynamics._next_time_step()

    '''
    Continue on this
    if error_time.history[-1] > 30:
        error_time.update_history(0)
        parafoil_dynamics.acc = np.array([0,0,0])
        calc_dubin = True 
        parafoil.TE_DEF = 0
        parafoil.bank = 0
        TE_deflection.update_history(0)'''

    if sqrt(parafoil_dynamics.acc.dot(parafoil_dynamics.acc)) < .5 and calc_dubin:
        #calculate maximum bank angle during turn
        TE_temp = pi/2
        sigma_maxx = np.arctan2(sqrt(parafoil_dynamics.vel_mag * parafoil_dynamics.vel_mag * 0.71 * TE_temp), sqrt(9.81 * parafoil.b))

        #set current position to x,y zero in the reference system
        #parafoil_dynamics.pos = np.array([0,0,parafoil_dynamics.pos[2]])

        #calculate trajectory
        minimum_conditions = _All_Dubin_Paths(pos_init=np.array([parafoil_dynamics.pos[0],parafoil_dynamics.pos[1],parafoil_attitude.psi]), 
        pos_final=np.array([200,200,-pi/2]), 
        altitude=parafoil_dynamics.pos[2],sigma_max=sigma_maxx,
        v_g=parafoil_dynamics.vel_mag,
        gamma_g_traj=parafoil_dynamics.gamma)
        minimum_conditions._Minimum_Tau()
        
        #initiate control
        TE = abs(math.cos(parafoil_dynamics.gamma)*parafoil.b/(0.71*minimum_conditions.r_traj))

        controls = True
        calc_dubin = False


    #update alpha var
    alpa.update_history(sqrt(parafoil_dynamics.acc.dot(parafoil_dynamics.acc)))

    #update position vars
    pos_x.update_history(parafoil_dynamics.pos[0])
    pos_y.update_history(parafoil_dynamics.pos[1])
    pos_z.update_history(parafoil_dynamics.pos[2])

    #update velocity vars
    vel_x.update_history(parafoil_dynamics.vel_r[0])
    vel_y.update_history(parafoil_dynamics.vel_r[1])
    vel_z.update_history(parafoil_dynamics.vel_r[2])

    #update noisy variables
    pos_x_noise.update_history(parafoil_dynamics.pos_noise[0])
    pos_y_noise.update_history(parafoil_dynamics.pos_noise[1])
    pos_z_noise.update_history(parafoil_dynamics.pos_noise[2])

    vel_x_noise.update_history(parafoil_dynamics.vel_noise[0])
    vel_y_noise.update_history(parafoil_dynamics.vel_noise[1])
    vel_z_noise.update_history(parafoil_dynamics.vel_noise[2])

    print(parafoil.Cl/parafoil.Cd)
    #update counter
    ts+= 0.01
    start=False

pos_x.plot(pos_y.history, None, "x", pos_y.var_name, True)
error_time.plot(None, None, "x", None, False)

pos_x_noise.plot(None, None, "x", None, False)
pos_y_noise.plot(None, None, "x", None, False)
pos_z_noise.plot(None, None, "x", None, False)
vel_x_noise.plot(None, None, "x", None, False)
vel_y_noise.plot(None, None, "x", None, False)
vel_z_noise.plot(None, None, "x", None, False)

fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter3D(pos_x.history, pos_y.history, pos_z.history, c=pos_z.history, cmap='Greens');
#ax.scatter3D(minimum_conditions.pos_x, minimum_conditions.pos_y, minimum_conditions.alt, c=minimum_conditions.alt, cmap='Greens');
#plt.savefig("First")
'''
plt.plot(minimum_conditions.pos_x, minimum_conditions.pos_y)
Control_Input.plot(None, None, "x", None, False)
Plots for control:
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.plot(pos_x.history, pos_y.history)
ax1.plot(minimum_conditions.pos_x, minimum_conditions.pos_y)
ax2.plot(error_time.history)
plt.show()

plt.plot(pos_x.history, pos_z.history)
plt.plot(minimum_conditions.pos_x, minimum_conditions.alt)
plt.show()

plt.plot(minimum_conditions.heading)
plt.plot(Psi.history)'''

plt.show()
    