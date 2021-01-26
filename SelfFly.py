import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, tan, pi, tanh, atan, asin, acos, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d

density = 1.225 #kg/m^3

class ParafoilProperties():
    def __init__(self, AR=2.5, alpha_0=0, a_0=2*pi, tau=0.25, epsilon=5*pi/180, surface=5, Cd0=0.01, delta=0.02, rigging=0, m=10):
        #parameters set for thin airfoil theory (see DARE-PRG_R2B Report and Anderson)
        self.Parafoil_Forces = np.array([0,0,0])
        self.Parafoil_Moments = np.array([0,0,0])
        self.Cl = 0
        
        #geometric properties
        self.AR = AR
        self.anhedral = epsilon
        self.surface = surface
        self.rigging = rigging
        self.m = m
        
        #airfoil properties
        self.a_0 = a_0*2*pi*AR*tanh(a_0/(2*pi*AR))/a_0 #reduction for LOW AR wing
        self.alpha_0 = alpha_0 #radians
        self.tau = tau #function of Fourier coefficients A_n
        self.a = a_0/(1+(1+tau)*(a_0/(pi*AR)))
        self.Cd_0 = Cd0 #find for airfoil shape
        self.delta = delta #estimate using 5.20 in Anderson, function of taper ratio
        
        #control properties
        self.Right_TE = 0
        self.Left_TE = 0

    #def alpha_opt(): calculate optimal angle of attack, defining rigging etc based on that

    
    def _Calc_Lift(self, alpha, velocity):
        #alpha in radians
        k1 = (3.33-1.33*self.AR) #for 1 < alpha < 2.5
        delta_cl = k1*(sin(alpha-self.alpha_0)**2)*cos(alpha-self.alpha_0)
        self.Cl = self.a * (alpha-self.alpha_0) * cos(self.anhedral)**2 + delta_cl
        #calculate total force
        return 0.5 * density * velocity**2 * self.Cl * self.surface
 
    def _Calc_Drag(self, alpha, velocity):
        #add payload and line drag contribution
        k1 = (3.33-1.33*self.AR) #for 1 < alpha < 2.5
        delta_cd = k1*sin(alpha-self.alpha_0)**3
        Cd = delta_cd + self.Cd_0 + (1+self.delta) * self.Cl**2 / (pi * self.AR)
        return 0.5 * density * velocity**2 * Cd * self.surface
    
    def _Calc_Pitch(self, velocity):
        Cm = 0
        return 0.5 * density * velocity**2 * Cm * self.surface

    def _Parafoil_Forces_Moments(self, alpa, vel, gamma):
        L = self._Calc_Lift(alpa, vel)*cos(gamma) + self._Calc_Drag(alpa, vel)*sin(gamma)
        D = self._Calc_Drag(alpa, vel)*cos(gamma) - self._Calc_Lift(alpa, vel)*sin(gamma)
        self.Parafoil_Forces = np.array([-D ,0, -L])
        Cm_pitch = self._Calc_Pitch(0)
        self.Parafoil_Moments = np.array([0,Cm_pitch,0])

    def _Parafoil_Control(self, turn_velocity):
        #trailing edge deflection in radians
        #see if it is possible to cause an increase in drag if both are deflected
        span = sqrt(self.AR*self.surface)
        if self.Left_TE != 0 and self.Right_TE == 0:
            return np.array([0,0, 0.71 * turn_velocity * self.Left_TE / span])
        elif self.Right_TE != 0 and self.Left_TE == 0:
            return np.array([0,0, -0.71 * turn_velocity * self.Right_TE / span])
        else:
            return np.array([0,0,0])

class Payload():
    def __init__(self, M=10, shape="", payload_cd=0):
        self.M=M
        self.shape = shape
        self.payload_cd = payload_cd

class Variable():
    def __init__(self, var_name=""):
        self.history = []
        self.var_name = var_name
    
    def update_history(self, value):
        self.history.append(value)

    def plot(self, subplot_one, subplot_two, x_or_y, ylabel):
        #input lists of data
        if subplot_one == None:
            fig = plt.figure()
            plt.plot(self.history)
            plt.show()
        elif x_or_y == "x":
            
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(self.history, subplot_one)
            plt.xlabel(self.var_name)
            plt.ylabel(ylabel)
            
            if subplot_two != None:
                ax2 = fig.add_subplot(112)
                ax2.plot(self.history, subplot_two)
            
            plt.show()

        elif x_or_y == "y":
            
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(subplot_one, self.history)
            plt.ylabel(self.var_name)
            plt.xlabel(ylabel)

            if subplot_two != None:
                ax2 = fig.add_subplot(112)
                ax2.plot(subplot_two, self.history)

            plt.show()

        else:
            print("Wrong")

    def plot_3D(self): 
        #plotting more variables
        print("Empty")

class Quaternion():
    def __init__(self, omega=np.array([0,0,0])):
        self.quaternion = np.array([1,0,0,0])
        self.dt = 0.05 #time interval of integration
        self.phi = 0
        self.theta = 0
        self.psi = 0 
        self.omega = omega

    def _to_quaternion(self, euler):
        """
        Set value of attitude quaternion from euler angles.

        :param euler: ([float]) euler angles roll, pitch, yaw.
        """
        self.phi, self.theta, self.psi = euler
        e0 = np.cos(self.psi / 2) * np.cos(self.theta / 2) * np.cos(self.phi / 2) + np.sin(self.psi / 2) * np.sin(self.theta / 2) * np.sin(
            self.phi / 2)
        e1 = np.cos(self.psi / 2) * np.cos(self.theta / 2) * np.sin(self.phi / 2) - np.sin(self.psi / 2) * np.sin(self.theta / 2) * np.cos(
            self.phi / 2)
        e2 = np.cos(self.psi / 2) * np.sin(self.theta / 2) * np.cos(self.phi / 2) + np.sin(self.psi / 2) * np.cos(self.theta / 2) * np.sin(
            self.phi / 2)
        e3 = np.sin(self.psi / 2) * np.cos(self.theta / 2) * np.cos(self.phi / 2) - np.cos(self.psi / 2) * np.sin(self.theta / 2) * np.sin(
            self.phi / 2)
        
        self.quaternion = np.array([e0, e1, e2, e3])

    #def _to_euler(self, attitude):

    def _rot_b_v(self):
        """
        Rotate vector from body frame to vehicle frame.

        :param Theta: ([float]) vector to rotate, either as Euler angles or quaternion
        :return: ([float]) rotated vector
        if len(attitude) == 3:
            phi, th, psi = attitude
            return np.array([
                [np.cos(th) * np.cos(psi), np.cos(th) * np.sin(psi), -np.sin(th)],
                [np.sin(phi) * np.sin(th) * np.cos(psi) - np.cos(phi) * np.sin(psi),
                 np.sin(phi) * np.sin(th) * np.sin(psi) + np.cos(phi) * np.cos(psi), np.sin(phi) * np.cos(th)],
                [np.cos(phi) * np.sin(th) * np.cos(psi) + np.sin(phi) * np.sin(psi),
                 np.cos(phi) * np.sin(th) * np.sin(psi) - np.sin(phi) * np.cos(psi), np.cos(phi) * np.cos(th)]
            ])"""

        
        e0, e1, e2, e3 = self.quaternion
        transfer = np.array([[-1 + 2 * (e0 ** 2 + e1 ** 2), 2 * (e1 * e2 + e3 * e0), 2 * (e1 * e3 - e2 * e0)],
                            [2 * (e1 * e2 - e3 * e0), -1 + 2 * (e0 ** 2 + e2 ** 2), 2 * (e2 * e3 + e1 * e0)],
                            [2 * (e1 * e3 + e2 * e0), 2 * (e2 * e3 - e1 * e0), -1 + 2 * (e0 ** 2 + e3 ** 2)]])
        return transfer

    def _update_quaternion(self):
        def _f_attitude_dot(t, y):
            """
            Right hand side of quaternion attitude differential equation.

            :param t: (float) time of integration
            :param attitude: ([float]) attitude quaternion
            :param omega: ([float]) angular velocity
            :return: ([float]) right hand side of quaternion attitude differential equation.
            """
            p, q, r = self.omega
            T = np.array([[0, -p, -q, -r],
                            [p, 0, r, -q],
                            [q, -r, 0, p],
                            [r, q, -p, 0]
                            ])
            return np.concatenate([0.5 * np.dot(T, y)])
        #print(_f_attitude_dot(0,self.quaternion))

        my_solution = solve_ivp(fun=_f_attitude_dot, t_span=(0, self.dt), y0=self.quaternion)
        self.quaternion = my_solution.y[:, -1]
        #print(self.quaternion)
'''
For now, add yaw as an angular rate
'''

class Dynamics:
    def __init__(self, dt=0.1, I=np.array([[1,0,0],[0,1,0],[0,0,1]]), mass=10, pos=np.array([0,0,100])):
        #properties, I is inertia matrix, mass in Newton
        self.I = I
        self.mass = mass
        #translational
        self.pos = pos #reference system
        self.vel = np.array([20,0,0]) #body system
        self.acc = np.array([0,0,0]) #body system
        #attitude, radians
        self.angular_velocity = np.array([[0,0,0]])
        self.angular_acceleration = np.array([[0,0,0]])
        #position, attitude log
        self.dt = dt
        self.time = 0
        #forces, moments, 3x1 matrices
        self.forces = np.array([0,0,-self.mass*9.81])
        self.moments = np.array([0,0,0])
        #information for turn
        self.vel_mag = 0
        self.gamma = 0
        self.turn_vel = 0

    def update_dynamics(self, rot_bv):
        self.forces = np.dot(rot_bv, self.forces) * np.array([1,1,0]) + self.forces * np.array([0,0,1]) + np.array([0,0,self.mass*9.81])
        #translational acceleration
        self.acc = self.forces * 1/self.mass
        self.vel = self.vel + self.acc * self.dt
        self.vel_mag = np.sqrt(self.vel.dot(self.vel))
        vel_reference = self.vel *  np.array([1,1,-1]) #switch from right hand to attitude positive upwards
        #translation
        self.pos = self.pos + vel_reference*self.dt + 0.5 * (self.dt**2) * self.acc * np.array([1,1,-1])
        #attitude
        self.angular_velocity = np.add(self.angular_velocity, np.dot(self.dt, self.angular_acceleration))
        self.gamma = atan( self.vel[2] / sqrt(self.vel[0]**2+self.vel[1]**2))
        self.turn_vel =  self.vel_mag * cos(self.gamma)

        #update time
        #self.time += self.dt

    def _next_time_step(self):
        self.time += self.dt
        
        #self.log.append(self.forces[1])
        #self.pos_log = np.append(self.pos_log, [self.pos], axis=0)
        #self.vel_log = np.append(self.vel_log, [self.vel], axis=0)
        #self.log.append([self.position, self.velocity, self.acceleration, self.time])
        
    def reset(self):
        self.log = None
        self.time = 0

