import SelfFly
import math
from math import pi
from matplotlib import pyplot as plt
from SelfFly import ParafoilProperties, Quaternion, Dynamics, Variable, Payload

ts = 0.025
#define parafoil-payload system properties
parafoil = ParafoilProperties(m=635, alpha_0=(-3.5*pi/180), surface=697, a_0=6.41, Cd0=0.005,rigging=0)
mpayload = Payload(M=9979)

#start dynamics and attitude tracking
parafoil_dynamics = Dynamics(dt=ts, mass=(parafoil.m+mpayload.M))
parafoil_attitude = Quaternion()

#system two
parafoil2 = ParafoilProperties(m=635, alpha_0=(-3.5*pi/180), surface=697, a_0=6.41, Cd0=0.005, rigging=10*pi/180)
parafoil_dynamics_2 = Dynamics(dt=ts, mass=(parafoil2.m+mpayload.M))

#intialize variables to be tracked
pos_x = Variable("Position [x]")
pos_y = Variable("Position [y]")
pos_z = Variable("Position [z]")
vel_x = Variable("Velocity [x]")
vel_y = Variable("Velocity [y]")
vel_z = Variable("Velocity [z]")

vel_x_2 = Variable("Velocity 2 [x]")

force_x = Variable("Drag [x]")
force_z = Variable("Lift [z]")

parafoil.Left_TE = 0#5*pi/180
change_TE  = False

start=True

while start == True or pos_z.history[-1] > 0:
    
    if parafoil_dynamics.pos[2] <= 150 and change_TE:
        parafoil.Left_TE = 0
        parafoil.Right_TE = 5*pi/180
        change_TE = False

    #update quaternion and gravity matrix
    parafoil_attitude.omega = parafoil._Parafoil_Control(parafoil_dynamics.turn_vel)
    parafoil_attitude._update_quaternion()
    parafoil_attitude._to_euler()

    #update parafoil forces and moments
    Parafoil_Vector = parafoil._Parafoil_Forces_Moments(parafoil_dynamics.vel)
    Payload_Vector = mpayload._Calc_Forces(parafoil_dynamics.vel)
    Gravity_Vector = parafoil_attitude.body_g * 9.80665 * parafoil_dynamics.mass
    parafoil_dynamics.forces = Parafoil_Vector + Payload_Vector + Gravity_Vector
    
    #second parafoil
    Parafoil_Vector_2 = parafoil2._Parafoil_Forces_Moments(parafoil_dynamics_2.vel)
    Payload_Vector_2 = mpayload._Calc_Forces(parafoil_dynamics_2.vel)
    Gravity_Vector_2 = parafoil_attitude.body_g * 9.80665 * parafoil_dynamics.mass
    parafoil_dynamics_2.forces = Parafoil_Vector_2 + Payload_Vector_2 + Gravity_Vector_2

    parafoil_dynamics_2.update_dynamics(parafoil_attitude._rot_b_v())
    parafoil_dynamics_2._next_time_step()

    # input them into Dynamics
    parafoil_dynamics.update_dynamics(parafoil_attitude._rot_b_v())
    parafoil_dynamics._next_time_step()
    
    #update position vars
    pos_x.update_history(parafoil_dynamics.pos[0])
    pos_y.update_history(parafoil_dynamics.pos[1])
    pos_z.update_history(parafoil_dynamics.pos[2])

    #update velocity vars
    vel_x.update_history(parafoil_dynamics.vel_r[0])
    vel_y.update_history(parafoil_dynamics.vel_r[1])
    vel_z.update_history(parafoil_dynamics.vel_r[2])

    force_x.update_history(parafoil.Parafoil_Forces[0])
    force_z.update_history(parafoil.Parafoil_Forces[2])

    vel_x_2.update_history(parafoil_dynamics_2.vel_r[0])

    #update counter
    ts+= 0.025
    start=False

fig = plt.figure()
plt.plot(vel_x.history)
plt.plot(vel_x_2.history)
plt.show()

pos_x.plot(pos_y.history, None, "x", pos_y.var_name)
pos_z.plot(pos_x.history, None, "y", pos_x.var_name)
vel_x.plot(None, None, "y", vel_x.var_name)
#vel_y.plot(None, None, "y", vel_y.var_name)
vel_z.plot(None, None, "y", vel_z.var_name)
force_x.plot(None, None, "y", force_x.var_name)
force_z.plot(None, None, "y", force_z.var_name)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter3D(pos_x.history, pos_y.history, pos_z.history, c=pos_z.history, cmap='Greens');
#plt.savefig("First")
plt.show()

mm = 1