import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, tan, pi
from scipy.integrate import solve_ivp

class ParafoilProperties():
    def init(self):
        self.Cl = None
        self.Cd = None
        self.CM_PITCH = None
        self.CM_YAW = None
        self.CM_ROLL = None

    #def set_properties(self, alpha):


class Quaternion():
    def __init__(self, omega=np.array([0,0,0])):
        self.quaternion = np.array([0,0,0,0])
        self.dt = 0.1 #time interval of integration
        self.phi = 0
        self.theta = 0
        self.psi = 0 
        self.omega = omega

    def _f_attitude_dot(self, t, omega_s):
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
        return 0.5 * np.dot(T, omega_s)

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

        if self.quaternion:
            e0, e1, e2, e3 = self.quaternion
            return np.array([[-1 + 2 * (e0 ** 2 + e1 ** 2), 2 * (e1 * e2 + e3 * e0), 2 * (e1 * e3 - e2 * e0)],
                             [2 * (e1 * e2 - e3 * e0), -1 + 2 * (e0 ** 2 + e2 ** 2), 2 * (e2 * e3 + e1 * e0)],
                             [2 * (e1 * e3 + e2 * e0), 2 * (e2 * e3 - e1 * e0), -1 + 2 * (e0 ** 2 + e3 ** 2)]])

        else:
            raise ValueError("Attitude is not Quaternion")    
    
    def _update_quaternion(self):
        my_solution = solve_ivp(fun=lambda t, y: _f_attitude_dot(t,y), t_span=(0, self.dt), y0=self.quaternion)
        self.quaternion = my_solution.y[0]

class Dynamics:
    def __init__(self, dt=0.1, I=np.array([[1,0,0],[0,1,0],[0,0,1]]), mass=1):
        #properties, I is inertia matrix, mass in Newton
        self.I = I
        self.mass = mass
        #translational
        self.pos = np.array([0,0,20]) #reference system
        self.vel = np.array([0,0,0]) #body system
        self.acc = np.array([0,0,0]) #body system
        #attitude, radians
        self.angular_velocity = np.array([[0,0,0]])
        self.angular_acceleration = np.array([[0,0,0]])
        #position, attitude log
        self.log = []
        self.dt = dt
        self.time = 0
        #forces, moments, 3x1 matrices
        self.forces = np.array([0,0,-self.mass*9.81])
        self.moments = np.array([0,0,0])

    def update_dynamics(self, rot_bv):
        #translational acceleration
        self.acc = self.forces * 1/self.mass
        self.vel = self.vel.__add__(self.acc * self.dt)
        vel_reference = rot_bv * self.vel
        #translation
        self.pos = self.pos + vel_reference*self.dt
        #attitude
        self.angular_velocity = np.add(self.angular_velocity, np.dot(self.dt, self.angular_acceleration))

    def _next_time_step(self):
        self.log.append(self.forces[1])
        self.time += self.dt
        #self.log.append([self.position, self.velocity, self.acceleration, self.time])
        
    def reset(self):
        self.log = None
        self.time = 0

parafoil = ParafoilProperties()
parafoil_dynamics = Dynamics()
parafoil_attitude = Quaternion()

t = 0.1
unit_vector = np.array([1,0,0])
omega_sim = np.array([0,0,pi/20])

while t < 0.2:
    parafoil_attitude.omega = omega_sim
    parafoil_attitude._update_quaternion()
    print(parafoil_attitude.quaternion)
    #unitvector = np.dot(parafoil_attitude.quaternion, unit_vector)
    t+= 0.1

"""
while parafoil_dynamics.pos[3] > 0:
    parafoil_dynamics.forces = np.array([0,0,0])"""