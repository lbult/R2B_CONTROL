import SelfFly
import math
import numpy as np
from math import pi, sqrt, asin
from matplotlib import pyplot as plt
from SelfFly import ParafoilProperties, Quaternion, Dynamics, Variable, Payload
from Wind_Reference_System import _All_Dubin_Paths

ts = 0.05
#define parafoil-payload system properties
parafoil = ParafoilProperties(m=635, alpha_0=(-3.5*pi/180), surface=697, a_0=6.41, Cd0=0.005,rigging=10*pi/180, ts=ts)
mpayload = Payload(M=9979)

#start dynamics and attitude tracking
parafoil_dynamics = Dynamics(dt=ts, mass=(parafoil.m+mpayload.M))
parafoil_attitude = Quaternion()
parafoil_attitude.psi = 3.14
parafoil_attitude._to_quaternion()

#system two
parafoil2 = ParafoilProperties(m=635, alpha_0=(-3.5*pi/180), surface=697, a_0=6.41, Cd0=0.005, rigging=10*pi/180, ts=ts)
parafoil_dynamics_2 = Dynamics(dt=ts, mass=(parafoil2.m+mpayload.M))

#intialize variables to be tracked
alpa = Variable("alpha")
pos_x = Variable("Position [x]")
pos_y = Variable("Position [y]")
pos_z = Variable("Position [z]")
vel_x = Variable("Velocity [x]")
vel_y = Variable("Velocity [y]")
vel_z = Variable("Velocity [z]")

vel_x_2 = Variable("Velocity 2 [x]")

force_x = Variable("Drag [x]")
force_z = Variable("Lift [z]")

change_TE_1  = True
change_TE_2 = True
change_TE_3 = True
current_alt = 0

calc_dubin = True

start=True

while start == True or pos_z.history[-1] > 0:
    '''
    if parafoil_dynamics.time > 2.0 and calc_dubin:
        gamma_traje = abs(np.arctan2(parafoil_dynamics.vel_r[2],sqrt((parafoil_dynamics.vel_r[0])**2+(parafoil_dynamics.vel_r[1])**2)))
        V_g = sqrt((parafoil_dynamics.vel_r[0])**2+(parafoil_dynamics.vel_r[1])**2+(parafoil_dynamics.vel_r[2])**2)

        parafoil.Left_TE = pi/2
        sigma_maxx = np.arcsin(sqrt(V_g * parafoil._Parafoil_Control(V_g)[2]/ 9.81))
        parafoil.Left_TE = 0

        minimum_conditions = _All_Dubin_Paths(pos_init=np.array([0,0,0]), pos_final=np.array([100,100,0]), 
        altitude=parafoil_dynamics.pos[2],sigma_max=sigma_maxx,v_g=V_g,gamma_g_traj=gamma_traje
        )
        minimum_conditions._Minimum_Tau()

        current_alt = parafoil_dynamics.pos[2]
        TE = parafoil.b* 9.81 *np.sin(minimum_conditions.sigma_max)/(0.71*V_g**2*np.cos(gamma_traje))
        
        calc_dubin = False

    try:
        if parafoil_dynamics.pos[2] == current_alt and change_TE_1:
            if minimum_conditions.chosen_traj[0] < 0:
                parafoil.Left_TE = TE
                change_TE_1 = False
            else:
                parafoil.Right_TE = TE
                change_TE_1 = False
    except:
        xyx =1

    try:
        if parafoil_dynamics.pos[2] < current_alt-abs(minimum_conditions.chosen_traj[0]) and parafoil_dynamics.pos[2] > current_alt -abs(minimum_conditions.chosen_traj[0])-abs(minimum_conditions.chosen_traj[1]) and change_TE_3:
            parafoil.Right_TE = 0
            parafoil.Left_TE = 0
            change_TE_3 = False
    except:
        xyx =1

    try:
        if parafoil_dynamics.pos[2] < current_alt-abs(minimum_conditions.chosen_traj[1])-abs(minimum_conditions.chosen_traj[0]) and change_TE_2:
            parafoil.Right_TE = 0
            parafoil.Left_TE = 0
            if minimum_conditions.chosen_traj[2] < 0:
                parafoil.Left_TE = TE
                change_TE_2 = False
            else:
                parafoil.Right_TE = TE
                change_TE_2 = False
    except:
        xyx = 1'''

    #update quaternion and gravity matrix
    parafoil_attitude.omega = parafoil._Parafoil_Control(parafoil_dynamics.turn_vel)
    parafoil_attitude.theta = parafoil.gamma + parafoil.rigging + parafoil.alpa
    parafoil_attitude._to_quaternion()
    parafoil_attitude._update_quaternion()
    parafoil_attitude._to_euler()

    #update parafoil forces and moments
    Payload_Vector = mpayload._Calc_Forces(parafoil_dynamics.vel)
    Parafoil_Vector = parafoil._Parafoil_Forces_Moments(parafoil_dynamics.vel, Payload_Vector, mpayload.M)
    Gravity_Vector = parafoil_attitude.body_g * 9.80665 * parafoil_dynamics.mass
    parafoil_dynamics.forces = Parafoil_Vector + Payload_Vector + Gravity_Vector
    
    #print(parafoil_dynamics.forces)
    #second parafoil
    #Parafoil_Vector_2 = parafoil2._Parafoil_Forces_Moments(parafoil_dynamics_2.vel)
    #Payload_Vector_2 = mpayload._Calc_Forces(parafoil_dynamics_2.vel)
    #Gravity_Vector_2 = parafoil_attitude.body_g * 9.80665 * parafoil_dynamics.mass
    #parafoil_dynamics_2.forces = Parafoil_Vector_2 + Payload_Vector_2 + Gravity_Vector_2

    #parafoil_dynamics_2.update_dynamics(parafoil_attitude._rot_b_v())
    #parafoil_dynamics_2._next_time_step()

    # input them into Dynamics
    parafoil_dynamics.update_dynamics(parafoil_attitude._rot_b_v())
    parafoil_dynamics._next_time_step()

    #update alpha var
    #alpa.update_history(parafoil.alpa)

    #update position vars
    pos_x.update_history(parafoil_dynamics.pos[0])
    pos_y.update_history(parafoil_dynamics.pos[1])
    pos_z.update_history(parafoil_dynamics.pos[2])

    #update velocity vars
    vel_x.update_history(parafoil_dynamics.vel[0])
    vel_y.update_history(parafoil_dynamics.vel[1])
    vel_z.update_history(parafoil_dynamics.vel[2])

    #update force vars
    force_x.update_history(parafoil.Parafoil_Forces[0])
    force_z.update_history(parafoil.Parafoil_Forces[2])
    print(parafoil.gamma)
    #update counter
    ts+= 0.05
    start=False

fig = plt.figure()
plt.plot(vel_x.history)


alpa.plot(None, None, "y", alpa.var_name, False)
pos_x.plot(pos_y.history, None, "x", pos_y.var_name, False)
pos_z.plot(pos_x.history, None, "y", pos_x.var_name, False)
vel_x.plot(None, None, "y", vel_x.var_name, False)
#vel_y.plot(None, None, "y", vel_y.var_name)
vel_z.plot(None, None, "y", vel_z.var_name, False)
force_x.plot(None, None, "y", force_x.var_name, False)
force_z.plot(None, None, "y", force_z.var_name, False)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter3D(pos_x.history, pos_y.history, pos_z.history, c=pos_z.history, cmap='Greens');
#plt.savefig("First")

plt.show()

mm = 1