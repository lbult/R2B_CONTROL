import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, tan, pi, tanh
from scipy.integrate import solve_ivp

density = 1.225 #kg/m^3

class ParafoilProperties():
    def __init__(self, AR=2.5, alpha_0=0, a_0=2*pi, tau=0.25, epsilon=5*pi/180, surface=5, Cd0=0.01, delta=0.02):
        #parameters set for thin airfoil theory (see DARE-PRG_R2B Report and Anderson)
        self.Parafoil_Forces = np.array([0,0,0])
        self.Parafoil_Moments = np.array([0,0,0])
        self.AR = AR
        self.alpha_0 = alpha_0 #radians
        self.a_0 = a_0*2*pi*AR*tanh(a_0/(2*pi*AR))/a_0 #reduction for LOW AR wing
        self.tau = tau #function of Fourier coefficients A_n
        self.a = a_0/(1+(1+tau)*(a_0/(pi*AR)))
        self.anhedral = epsilon
        self.surface = surface
        self.Cd_0 = Cd0 #find for airfoil shape
        self.delta = delta #estimate using 5.20 in Anderson, function of taper ratio
        self.Cl = 0
    
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

    def _Parafoil_Forces_Moments(self, alpa, vel):
        L = self._Calc_Lift(alpa, vel)
        D = self._Calc_Drag(alpa, vel)
        self.Parafoil_Forces = np.array([-D ,0, -L])
        Cm_pitch = self._Calc_Pitch(0)
        self.Parafoil_Moments = np.array([0,Cm_pitch,0])

class Quaternion():
    def __init__(self, omega=np.array([0,pi/20,0])):
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

    def _rot_b_v(self, vector):
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
        return np.dot(transfer, vector) 

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
        vel_reference = rot_bv * self.vel * np.array([]) #switch from right hand to attitude positive upwards
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
parafoil._Parafoil_Forces_Moments(1/180*pi, 10)
print(parafoil.Parafoil_Forces)
parafoil_dynamics = Dynamics()
parafoil_attitude = Quaternion()

ts = 0.1
unit_vector = np.array([1,1,1])
mlist = []
klist = []
omega_sim = np.array([0,0,pi/20])

while ts < 1:
    #update parafoil forces and moments
    parafoil._Parafoil_Forces_Moments(1/180*pi, 10)
    # input them into Dynamics
    parafoil_dynamics.forces = parafoil #+add inverse rotation of gravitational vector
    #parafoil_attitude.omega = omega_sim
    parafoil_attitude._update_quaternion()
    unit_vector = parafoil_attitude._rot_b_v(unit_vector)
    #print(parafoil_attitude.quaternion)
    mlist.append(unit_vector[0])
    klist.append(unit_vector[1])
    #unitvector = np.dot(parafoil_attitude.quaternion, unit_vector)
    ts+= 0.05

#print(mlist)
#print(klist)

#plt.plot(mlist)
#plt.plot(klist)
#plt.show()

"""
while parafoil_dynamics.pos[3] > 0:
    parafoil_dynamics.forces = np.array([0,0,0])"""