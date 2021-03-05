import numpy as np
from math import sin, cos
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class Quaternion():
    '''
    Keeps track of the angular position of the spacecraft with respect to its starting
        attitude in the reference system. Avoids problems such as gimbal lock.
    '''
    def __init__(self, omega=np.array([0,0,0])):
        self.quaternion = np.array([1,0,0,0])
        self.dt = 0.025 #time interval of integration
        self.phi = 0
        self.theta = 0
        self.psi = 0 
        self.omega = omega
        self.body_g = np.array([0,0,0])

    def _to_quaternion(self):
        """
        Set value of attitude quaternion from euler angles.
        :param psi: ([float]) euler angle yaw
        :param phi: ([float]) euler angle roll
        :param theta: ([float]) euler angle pitch
        :return: ndarray(4), quaternion of form: e0 + e1*i + e2*j + e3*k
        """
        e0 = np.cos(self.psi / 2) * np.cos(self.theta / 2) * np.cos(self.phi / 2) + np.sin(self.psi / 2) * np.sin(self.theta / 2) * np.sin(
            self.phi / 2)
        e1 = np.cos(self.psi / 2) * np.cos(self.theta / 2) * np.sin(self.phi / 2) - np.sin(self.psi / 2) * np.sin(self.theta / 2) * np.cos(
            self.phi / 2)
        e2 = np.cos(self.psi / 2) * np.sin(self.theta / 2) * np.cos(self.phi / 2) + np.sin(self.psi / 2) * np.cos(self.theta / 2) * np.sin(
            self.phi / 2)
        e3 = np.sin(self.psi / 2) * np.cos(self.theta / 2) * np.cos(self.phi / 2) - np.cos(self.psi / 2) * np.sin(self.theta / 2) * np.sin(
            self.phi / 2)
        
        self.quaternion = np.array([e0, e1, e2, e3])

    def _to_euler(self):
        '''
        Refer to _to_quaternion()
        :return self.body_g: rotates gravity vector into body system
        '''
        #set each individual element
        e0, e1, e2, e3 = self.quaternion
        #update roll, pitch, yaw
        self.phi = np.arctan2(2 * (e0 * e1 + e2 * e3), e0 ** 2 + e3 ** 2 - e1 ** 2 - e2 ** 2)
        self.theta = np.arcsin(2 * (e0 * e2 - e1 * e3))
        self.psi = np.arctan2(2 * (e0 * e3 + e1 * e2), e0 ** 2 + e1 ** 2 - e2 ** 2 - e3 ** 2)
        self.body_g = np.array([-sin(self.theta), sin(self.phi)*cos(self.theta), cos(self.phi)*cos(self.theta)])
        return

    def _rot_b_v(self):
        """
        Rotate vector from body frame to vehicle frame.
        :param quaternion: quaternion input to rotation vector
        :return: transfer vector, body --> reference
        """
        e0, e1, e2, e3 = self.quaternion
        transfer = np.array([[-1 + 2 * (e0 ** 2 + e1 ** 2), 2 * (e1 * e2 + e3 * e0), 2 * (e1 * e3 - e2 * e0)],
                            [2 * (e1 * e2 - e3 * e0), -1 + 2 * (e0 ** 2 + e2 ** 2), 2 * (e2 * e3 + e1 * e0)],
                            [2 * (e1 * e3 + e2 * e0), 2 * (e2 * e3 - e1 * e0), -1 + 2 * (e0 ** 2 + e3 ** 2)]])
        return transfer

    def _update_quaternion(self):
        """
        Right hand side of quaternion attitude differential equation.

        :param t: (float) time of integration
        :param attitude: ([float]) attitude quaternion
        :param omega: ([float]) angular velocity (x,y,z)
        :param T: skew symmetric matrix of angular velocities
        :func solve_ivp: solves the initial value problem x = e^(At) * x_0
        :return: ndarray(4), updated quaternion
        """
        def _f_attitude_dot(t, y):
            p, q, r = self.omega
            T = np.array([[0, -p, -q, -r],
                            [p, 0, r, -q],
                            [q, -r, 0, p],
                            [r, q, -p, 0]
                            ])
            return np.concatenate([0.5 * np.dot(T, y)])

        my_solution = solve_ivp(fun=_f_attitude_dot, t_span=(0, self.dt), y0=self.quaternion)
        self.quaternion = my_solution.y[:, -1]

class Variable():
    """
    Logger class, for storing and plotting any and all variables
    """
    def __init__(self, var_name="", limit=10**8):
        self.history = []
        self.var_name = var_name
        self.limit = limit
        
    def update_history(self, value):
        if value < self.limit and value > -self.limit:
            self.history.append(value)
        else:
            value = self.limit
            self.history.append(value)

    def plot(self, subplot_one, subplot_two, x_or_y, ylabel, normalize):
        #input lists of data
        if subplot_one == None:
            fig = plt.figure()
            if normalize:
                plt.gca().set_aspect('equal', adjustable='box')
            plt.plot(self.history)
            plt.xlabel("Time")
            plt.ylabel(self.var_name)
            plt.show()

        elif x_or_y == "x":
            
            fig = plt.figure()
            if normalize:
                plt.gca().set_aspect('equal', adjustable='box')
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
            if normalize:
                plt.gca().set_aspect('equal', adjustable='box')
            ax1 = fig.add_subplot(111)
            ax1.plot(subplot_one, self.history)
            plt.xlabel(self.var_name)
            plt.xlabel(ylabel)

            if subplot_two != None:
                ax2 = fig.add_subplot(112)
                if normalize:
                    plt.gca().set_aspect('equal', adjustable='box')
                ax2.plot(subplot_two, self.history)

            plt.show()

        else:
            print("Wrong")

    def plot_3D(self): 
        #plotting more variables, NOT FINISHED
        print("Empty")