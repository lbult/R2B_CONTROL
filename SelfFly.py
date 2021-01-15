import numpy as np
from Vector2 import Vector2
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, tan, pi


'''
class ParafoilProperties:
    def init(self):
        self.Cl = None
        self.Cd = None
        self.CM_PITCH = None
        self.CM_YAW = None
        self.CM_ROLL = None

    def set_properties(self, alpha):
'''

class Quaternion():
    def __init__(self):
        self.quaternion = np.array([0,0,0,0])
        self.delta_q = None        #derivative of quaternion
        self.dt = 0.1 #time interval of integration
        self.phi = 0
        self.theta = 0
        self.psi = 0 

    def _f_attitude_dot(t, attitude, omega):
        """
        Right hand side of quaternion attitude differential equation.

        :param t: (float) time of integration
        :param attitude: ([float]) attitude quaternion
        :param omega: ([float]) angular velocity
        :return: ([float]) right hand side of quaternion attitude differential equation.
        """
        p, q, r = omega
        T = np.array([[0, -p, -q, -r],
                        [p, 0, r, -q],
                        [q, -r, 0, p],
                        [r, q, -p, 0]
                        ])
        self.delta_q = 0.5 * np.dot(T, attitude)

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

    def _rot_b_v(self, attitude):
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

        if len(attitude) == 4:
            e0, e1, e2, e3 = attitude
            return np.array([[-1 + 2 * (e0 ** 2 + e1 ** 2), 2 * (e1 * e2 + e3 * e0), 2 * (e1 * e3 - e2 * e0)],
                             [2 * (e1 * e2 - e3 * e0), -1 + 2 * (e0 ** 2 + e2 ** 2), 2 * (e2 * e3 + e1 * e0)],
                             [2 * (e1 * e3 + e2 * e0), 2 * (e2 * e3 - e1 * e0), -1 + 2 * (e0 ** 2 + e3 ** 2)]])

        else:
            raise ValueError("Attitude is not Quaternion")    

    def update_quaternion(self):
        y0 = self.quaternion
        ivp_vars = np.concatenate([self.delta_q])
        self.quaternion = solve_ivp(fun=lambda t, y: ivp_vars, t_span=(0, self.delta_t), y0=y0)


class Dynamics:
    def __init__(self):
        #properties, I is inertia matrix, mass in Newton
        self.I = np.array([[10,1,3], [0,10,0], [0,0,10]])
        self.mass = 10
        #translational
        self.position = np.array([0,0,0])
        self.velocity = np.array([0,0,0])
        self.acceleration = np.array([0,0,0])
        #attitude, radians
        self.angular_velocity = np.array([[0,0,0]])
        self.angular_acceleration = np.array([[0,0,0]])
        #position, attitude log
        self.log = []
        self.dt = 0.1
        #forces, moments, 3x1 matrices
        self.time = 0
        self.forces = np.array([0,1,0])
        self.moments = np.array([0,0,0])

    def update

    def update_spacecraft(self):
        self.rotation_matrix = np.linalg.multi_dot(self.roll, self.pitch, self.yaw)

    def update_dynamics(self, forces, moments):
        #define forces, moments
        self.forces = forces
        self.moments = moments
        #angular acceleration
        #inverse_I = np.linalg.inv(self.I)
        #self.angular_acceleration = np.dot(inverse_I, moments - np.cross(self.angular_velocity, np.dot(self.I, self.angular_velocity)))
        #translational acceleration
        self.acceleration = self.forces * 1/self.mass

    def update_kinetics(self, dt):
        self.dt = dt
        #translation
        self.position = self.position + self.velocity * self.dt + self.acceleration * 0.5 * (self.dt)**2
        self.velocity = self.velocity.__add__(self.acceleration * self.dt)
        self.time += self.dt
        #attitude
        self.attitude = np.add(self.attitude, np.dot(self.dt, self.angular_velocity, np.dot(self.angular_acceleration, 0.5 * (self.dt)**2)))
        self.angular_velocity = np.add(self.angular_velocity, np.dot(self.dt, self.angular_acceleration))

    def update_log(self):
        self.log.append(self.forces[1])
        #self.log.append([self.position, self.velocity, self.acceleration, self.time])
        
    def reset(self):
        self.log = None
    

my_parafoil = Body()
my_parafoil.mass = 10
my_parafoil.I = np.array([[10,1,3], [0,10,0], [0,0,10]])
o = 1
of = 1
oof = 3
#my_parafoil.update_dynamics(forces=np.array([0, 0, -my_parafoil.mass*9.81]), moments=np.array([10,0,0]))


while o < oof:
    my_parafoil.update_rotation_matrix(3.1415/2)
    my_parafoil.update_log()
    my_parafoil.update_kinetics(0.1)
    o += of

plt.plot(my_parafoil.log)
plt.ylabel('some numbers')
plt.show()