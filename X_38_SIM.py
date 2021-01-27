import SelfFly
import math
from math import pi
from matplotlib import pyplot as plt
from SelfFly import ParafoilProperties, Quaternion, Dynamics, Variable, Payload

ts = 0.05
#define parafoil-payload system properties
parafoil = ParafoilProperties(m=635, alpha_0=(-3.5*pi/180), surface=697, a_0=6.41, Cd0=0.005)
mpayload = Payload(M=9979)

#start dynamics and attitude tracking
parafoil_dynamics = Dynamics(dt=ts, mass=(parafoil.m+mpayload.M))
parafoil_attitude = Quaternion()

#intialize variables to be tracked
pos_x = Variable("Position [x]")
pos_y = Variable("Position [y]")
pos_z = Variable("Position [z]")
vel_x = Variable("Velocity [x]")
vel_y = Variable("Velocity [y]")
vel_z = Variable("Velocity [z]")

force_x = Variable("Drag [x]")
force_z = Variable("Lift [z]")

parafoil.Left_TE = 0 #5*pi/180
change_TE  = False

start=True

while start == True or pos_z.history[-1] > 0:
    
    if ts >= 3.25 and change_TE:
        parafoil.Left_TE = 0
        parafoil.Right_TE = 5*pi/180
        change_TE = False
    
    #update parafoil forces and moments
    parafoil._Parafoil_Forces_Moments(3*pi/180, parafoil_dynamics.vel_mag, parafoil_dynamics.gamma)
    parafoil_dynamics.forces = parafoil.Parafoil_Forces
    print(parafoil.Parafoil_Forces)

    # input them into Dynamics
    parafoil_attitude.omega = parafoil._Parafoil_Control(parafoil_dynamics.turn_vel)
    parafoil_attitude._update_quaternion()
    parafoil_dynamics.update_dynamics(parafoil_attitude._rot_b_v())
    parafoil_dynamics._next_time_step()
    
    #update position vars
    pos_x.update_history(parafoil_dynamics.pos[0])
    pos_y.update_history(parafoil_dynamics.pos[1])
    pos_z.update_history(parafoil_dynamics.pos[2])
    
    #update velocity vars
    vel_x.update_history(parafoil_dynamics.vel[0])
    vel_y.update_history(parafoil_dynamics.vel[1])
    vel_z.update_history(parafoil_dynamics.vel[2])

    force_x.update_history(parafoil.Parafoil_Forces[0])
    force_z.update_history(parafoil.Parafoil_Forces[2])

    #update counter
    ts+= 0.05
    start=False

pos_x.plot(pos_y.history, None, "x", pos_y.var_name)
pos_z.plot(pos_x.history, None, "y", pos_x.var_name)
vel_x.plot(None, None, "y", vel_x.var_name)
vel_y.plot(None, None, "y", vel_y.var_name)
vel_z.plot(None, None, "y", vel_z.var_name)
force_x.plot(None, None, "y", force_x.var_name)
force_z.plot(None, None, "y", force_z.var_name)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter3D(pos_x.history, pos_y.history, pos_z.history, c=pos_z.history, cmap='Greens');
#plt.savefig("First")
plt.show()

mm = 1