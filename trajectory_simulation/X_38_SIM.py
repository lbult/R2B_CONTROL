import SelfFly
import math
import numpy as np
from math import pi, sqrt, asin
from matplotlib import pyplot as plt
from SelfFly import ParafoilProperties, Dynamics, Payload
from support import Quaternion, Variable
from Dubin_Path import _All_Dubin_Paths

ts = 0.02
#define parafoil-payload system properties
parafoil = ParafoilProperties(m=635, alpha_0=(-3.5*pi/180), surface=697, 
    a_0=6.41, Cd0=0.005,rigging=10*pi/180, ts=ts)
mpayload = Payload(M=9979,payload_cd=0.02)

#start dynamics and attitude tracking
parafoil_dynamics = Dynamics(dt=ts, mass=(parafoil.m+mpayload.M))
parafoil_attitude = Quaternion()
parafoil_attitude.psi = 3.14
parafoil_attitude._to_quaternion()

#intialize variables to be tracked
alpa = Variable("alpha")
pos_x = Variable("Position [x]")
pos_y = Variable("Position [y]")
pos_z = Variable("Position [z]")
vel_x = Variable("Velocity [x]")
vel_y = Variable("Velocity [y]")
vel_z = Variable("Velocity [z]")

error_time = Variable("Error")
error_time.update_history(0)

force_x = Variable("Drag [x]")
force_z = Variable("Lift [z]")

current_alt = 0

start = True
controls = False
calc_dubin = True

def takeClosest(myList, myNumber):
    closest = myList[0]
    for i in range(1, len(myList)):
        if abs(myList[i] - myNumber) < closest:
            closest = myList[i]
    return closest

while start or pos_z.history[-1] > 0:
    
    #after initialisation of onboard navigation
    if controls:
        index = minimum_conditions.alt.index(takeClosest(minimum_conditions.alt, 
            parafoil_dynamics.pos[2]*2))
        control_input = minimum_conditions.control[index]

        def _check_control(desired_radius, desired_heading, actual_position, arc_center, left_or_right):
            '''
            :param left_or_right: negative if left turn, positive if right turn
            :param actual_position: actual position ([0,1]) and heading ([2]), nparray(3)
            :param arc_center: desired arc center of turn, list(2)
            :param desired_radius: desired turn radius, float
            '''
            dx = arc_center[0] - actual_position[0]
            dy = arc_center[1] - actual_position[1]
            actual_radius = sqrt(dx**2 + dy**2)
            
            #set PID gains
            kp = 0.0005
            kd = 0.00001

            if left_or_right < 0:
                e = (desired_radius - actual_radius)*math.cos(desired_heading)
                de_dt = (error_time.history[-1]-e) / ts
            elif left_or_right > 0:
                e = (actual_radius - desired_radius)*math.cos(desired_heading)
                de_dt = (e-error_time.history[-1]) / ts
            
            delta_control = kp * e + kd*de_dt

            error_time.update_history(e)

            return delta_control

        #account for banking, which balances the centrifugal force
        if control_input < 0:
            parafoil.Left_TE = TE + -abs(_check_control(minimum_conditions.r_traj, minimum_conditions.heading[index], parafoil_dynamics.pos, [minimum_conditions.arc_centers[1][0], -1*minimum_conditions.arc_centers[1][1]], -1))
            parafoil.Right_TE = 0
            parafoil.bank = -minimum_conditions.sigma_max
        elif control_input > 0:
            parafoil.Left_TE = 0
            parafoil.Right_TE = TE + abs(_check_control(minimum_conditions.r_traj,  minimum_conditions.heading[index], parafoil_dynamics.pos, minimum_conditions.arc_centers[0], 1))
            parafoil.bank = minimum_conditions.sigma_max
        else:
            parafoil.Left_TE = 0
            parafoil.Right_TE = 0
            parafoil.bank = 0

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


    if sqrt(parafoil_dynamics.acc.dot(parafoil_dynamics.acc)) < .5 and calc_dubin:
        #calculate maximum bank angle during turn
        parafoil.Left_TE = pi/2
        sigma_maxx = np.arcsin(sqrt(parafoil_dynamics.vel_mag * parafoil._Parafoil_Control(parafoil_dynamics.vel_mag)[2]/ 9.81))
        parafoil.Left_TE = 0

        #set current position to x,y zero in the reference system
        parafoil_dynamics.pos = np.array([0,0,parafoil_dynamics.pos[2]])

        #calculate trajectory
        minimum_conditions = _All_Dubin_Paths(pos_init=np.array([0,0,0]), 
        pos_final=np.array([100,100,pi/2]), 
        altitude=parafoil_dynamics.pos[2],sigma_max=sigma_maxx,
        v_g=parafoil_dynamics.vel_mag,
        gamma_g_traj=parafoil_dynamics.gamma)
        minimum_conditions._Minimum_Tau()
        
        #initiate control
        TE = math.cos(parafoil_dynamics.gamma)*parafoil.b/(0.71*minimum_conditions.r_traj) #parafoil.b * 9.81 *np.sin(minimum_conditions.sigma_max)/(0.71*(minimum_conditions.v_min)**2*np.cos(parafoil_dynamics.gamma))
        print(minimum_conditions.chosen_traj)
        print(minimum_conditions.arc_centers)
        print(minimum_conditions.r_traj)

        print(minimum_conditions.sigma_max)
        parafoil.Left_TE = TE
        print(np.arcsin(parafoil_dynamics.vel_mag**2*math.cos(parafoil_dynamics.gamma)*0.71*TE/(9.81*parafoil.b)))
        parafoil.Left_TE = 0
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

    #update force vars
    force_x.update_history(Parafoil_Vector[2])
    force_z.update_history(parafoil.Parafoil_Forces[2])

    #update counter
    ts+= 0.01
    start=False

#parafoil._Calc_Alpha_Trim(0.02)

fig = plt.figure()

#alpa.plot(None, None, "y", alpa.var_name, False)
pos_x.plot(pos_y.history, None, "x", pos_y.var_name, True)
error_time.plot(None, None, "x", None, False)
#pos_z.plot(pos_x.history, None, "y", pos_x.var_name, True)
#vel_x.plot(None, None, "y", vel_x.var_name, False)
#vel_y.plot(None, None, "y", vel_y.var_name)
#vel_z.plot(None, None, "y", vel_z.var_name, False)
#force_x.plot(None, None, "y", force_x.var_name, False)
#force_z.plot(None, None, "y", force_z.var_name, False)

'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter3D(pos_x.history, pos_y.history, pos_z.history, c=pos_z.history, cmap='Greens');
#plt.savefig("First")
'''
plt.show()

mm = 1
    