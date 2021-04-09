import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, tan, pi, tanh, atan, asin, acos, sqrt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d

density = 1.225 #kg/m^3
cloth_density = 0.0678  # kg/m**2

def _Skew_Symmetric_Operator(vector_to_skew):
    return np.array([[0,-vector_to_skew[2],vector_to_skew[1]],
                [vector_to_skew[2],0,-vector_to_skew[0]],
                [-vector_to_skew[1],vector_to_skew[0],0]])

class ParafoilProperties():
    def __init__(self, alpha_0=0, a_0=2*pi, b=43.6, surface=5, Cd0=0.01, delta=0.02, rigging=0, m=10, R=70.3, line_n=697, line_d=2.5, thickness=0.4, ts=0.025):
        #parameters set for thin airfoil theory (see DARE-PRG_R2B Report and Anderson)
        self.Parafoil_Forces = np.array([0,0,0])
        self.Parafoil_Moments = np.array([0,0,0])
        self.Cl = 0
        self.Cd = 0

        #geometric properties
        self.AR = b**2/surface
        self.anhedral = b/(4*R)
        self.b = b
        self.surface = surface
        self.rigging = rigging*pi/180
        self.m = m
        self.t = thickness
        self.c = surface/b # TODO review this later, chord is calculated as surface area divided by span?, so basically assumes shape of parafoil?
        
        #airfoil properties
        self.a_0 = a_0*2*pi*self.AR*tanh(a_0/(2*pi*self.AR))/a_0 #reduction for LOW AR wing
        self.alpha_0 = alpha_0 #radians
        a = 0.03714, 1
        b = 0.14286, 5
        self.tau = (self.AR - a[1])*(b[0]-a[0])/(b[1]-a[1])+a[0] #correction factor for small aspect ratio wings
        self.a = a_0/(1+(1+self.tau)*(a_0/(pi*self.AR)))

        self.Cd_0 = Cd0 #find for airfoil shape
        self.delta = delta #estimate using 5.20 in Anderson, function of taper ratio
        self.app_MMOI_curv = np.array([0,0,0])
        self.gamma = 0

        #line properties
        self.R = R # mean lin length
        self.line_d = line_d/1000 # self.line_d in meters
        self.line_n = line_n
        self.Cdl = 0

        #control properties
        self.TE_DEF = 0 
        self.bank = 0

        # initial condition for angle of attack
        self.alpa = 6.385/180 * pi
        self.alpa_prime = 0

        self.ts = ts  # timestep
        # initialize apparent masses
        self.Parafoil_cloth_mass = 2.1*cloth_density*self.surface


    def _Calc_CG_height(self, payload_m):
        self.cg_height_wrtpayload = (self.Parafoil_cloth_mass*self.R )/(self.Parafoil_cloth_mass+ payload_m)
        return self.cg_height_wrtpayload
    
    def _Calc_Alpha_Trim(self, Cds, x_0):
        Xw, Zw, Zl, Zcg, Zp = self.c / 4, self.R * 0.02, self.R * 0.5, self.R - float(self._Calc_CG_height(9979)), self.R
        #normalise lengths
        Zl = Zl/self.c
        Zp = Zp/self.c
        zw = (Zw-Zcg)/self.c
        xw = (Xw)/self.c
        
        Delta = 1/(self.a) - (1+self.delta)/(pi*self.AR)
        Gamma = 1/(self.a) - 2*(1+self.delta)/(pi*self.AR)
        Cl0 = abs(self.alpha_0) * self.a

        Cm_0 = -Zl*self.Cdl + Zp*Cds + Cl0*xw - self.Cd_0*zw
        Cm_1 = self.a*(xw+Cl0*Delta*zw)
        Cm_2 = Gamma*zw*self.a**2

        def f(x):
            return Cm_0 + Cm_1*x + Cm_2*x**2 +0.09-2.1298*x+22.7498*x**2-649.379*x**3+4915*x**4-10814.1*x**5

        def g(x):
            return Cm_0 + Cm_1*x + Cm_2*x**2 + 69176.3*x**6-31385.3*x**5+4353.3*x**4-102.92*x**3 - 11.5626*x**2-1.04585*x-0.035
        
        def dgdx(x):
            return Cm_1 + 2*Cm_2*x + 6*69176.3*x**5-5*31385.3*x**4+4*4353.3*x**3-3*102.92*x**2-2*11.5626*x -1.04585
        
        def dfdx(x):
            return Cm_1 + 2*Cm_2*x - 2.1298+22.7498*2*x-3*649.379*x**2+4*4915*x**3-5*10814.1*x**4
        
        iteration = 0
        x_tilde = x_0
        while iteration < 100:
            x_tilde = x_tilde - g(x_tilde)/dgdx(x_tilde)
            iteration += 1
            

        x = np.linspace(-pi/38, pi/18, 1000)
        y = np.zeros(1000)
        counter = 0
        for i in x:
            y[counter] = f(i)
            counter+=1

        plt.plot(x,y)
        plt.show()
        return x_tilde


    def _Calc_Lift(self, alpha, velocity):
        """
        Function for calculating the lift force of the parafoil
        :param alpha: angle of attack of the parafoil, 
            angle between the chord line of parafoil 
            (horizontal axis of body reference system, and freestream air) 
            (horizontal axis of environment reference system, float)
        :param velocity: velocity of parafoil in the body reference system, nparray(3)
        :return: lift force of the parafoil, normal to chord (Fn), vertical in vehicle axis system, nparray(3)
        """
        #alpha in radians
        k1 = (3.33-1.33*self.AR) #for 1 < alpha < 2.5
        delta_cl = k1*(sin(alpha-self.alpha_0)**2)*cos(alpha-self.alpha_0)
        self.Cl = self.a * (alpha+self.rigging-self.alpha_0) * cos(self.anhedral)**2 + delta_cl
        
        return 0.5 * density * self.Cl * self.surface * np.array([-velocity[2], 0, velocity[0]*cos(self.bank)]) * sqrt(velocity[2]**2 + velocity[0]**2)
 
    def _Calc_Drag(self, alpha, velocity):
        """
        Function for calculating drag force of the parafoil.
        :param alpha: angle of attack of the parafoil, angle between the chord line of parafoil (horizontal axis of body reference system,
         and freestream air (horizontal axis of environment reference system, float
        :param velocity: velocity of parafoil in the body reference system, nparray(3)
        :return: drag force of the parafoil, horizontal in vehicle axis system, nparray(3)
        """
        #add payload and line drag contribution
        k1 = (3.33-1.33*self.AR) #for 1 < alpha < 2.5
        delta_cd = k1*sin(alpha+self.rigging-self.alpha_0)**3
        self.Cdl = self.line_n*self.R*self.line_d*cos(alpha)**3/self.surface
        self.Cd = delta_cd + self.Cd_0 + (1+self.delta) * self.Cl**2 / (pi * self.AR) + self.Cdl 
        
        return -0.5 * density * self.Cd * self.surface * sqrt(velocity.dot(velocity)) * velocity

    def _Parafoil_Forces_Moments(self, vel):
        """
        :param vel: velocity vector of parafoil, nparray(3)
        :return: Total aerodynamic force vector acting on parafoil, nparray(3)
        """
        self.Parafoil_Forces = self._Calc_Lift(self.alpa, vel) + self._Calc_Drag(self.alpa, vel)
        return self._Calc_Lift(self.alpa, vel) + self._Calc_Drag(self.alpa, vel)

    def _Parafoil_Control(self, turn_velocity):
        """
        :param turn_velocity: velocity tangent to turning arc, nparray(3)
        :return: yaw rate of the parafoil under turning, nparray(3)
        """
        #trailing edge deflection in radians
        if self.TE_DEF != 0:
            return np.array([0,0, 0.71 * turn_velocity * self.TE_DEF / self.b])
        else:
            return np.array([0,0,0])


class Payload():
    """
    :class calculates drag of the paylaod based on payload properties
    :param M: mass of the payload, int/float
    :param shape: shape of the object for Cd, NOT FINISHED
    :param payload_cd: drag coefficient of the object based on payload shape
    :param payload_surface: exposed frontal_area of the payload
    :return: 
    """
    def __init__(self, M=10, shape="", payload_cd=0.02, payload_surface=1):
        self.M=M
        self.shape = shape
        self.payload_cd = payload_cd
        self.payload_surface = payload_surface
        self.Payload_Forces = np.array([])

    def _Calc_Forces(self, velocity):
        self.Payload_Forces = -0.5 * density * self.payload_cd * velocity * sqrt(velocity.dot(velocity))
        return -0.5 * density * self.payload_cd * velocity * sqrt(velocity.dot(velocity)) * np.array([-1,1,1])

class Dynamics:
    def __init__(self, dt=0.025, mass=10, pos=np.array([0,0,500])):
        #translational
        self.pos = pos #reference system
        self.vel_r = np.array([0,0,0]) # reference system
        self.vel = np.array([-17.8,0,6.25]) #body system
        self.acc = np.array([0,0,0]) #body system

        #translational with noise added (Gaussian noise)
        self.pos_noise = 0
        self.vel_noise = 0

        #position, attitude log
        self.dt = dt
        self.time = 0
        self.mass = mass

        #forces, 3x1 matrix
        self.forces = np.array([0,0,0])

        #information for turn
        self.vel_mag = 0
        self.gamma = 0
        self.turn_vel = 0

    def update_dynamics(self, rot_bv): 
        #translational acceleration
        self.acc = self.forces * 1/self.mass
        self.vel = self.vel + self.acc * self.dt
        self.pos = self.pos + self.vel_r*self.dt 

        #derive noisy state from model state
        
        #pos_noise has x,y noise from GPS error
            #z noise from altitude sensor error
        self.pos_noise = self.pos + np.array([np.random.normal(0,3.5),
            np.random.normal(0,3.5), np.random.normal(0,1)])
        self.vel_noise = self.vel + np.array([np.random.normal(0,3.5),
        np.random.normal(0,3.5), 0])

        #reference system
        self.vel_r = np.dot(np.transpose(rot_bv), self.vel) * -1
        self.vel_mag = np.sqrt(self.vel.dot(self.vel))

        #attitude
        self.gamma = np.arctan2(self.vel_r[2], sqrt(self.vel_r[0]**2+self.vel_r[1]**2))
        self.turn_vel =  sqrt(self.vel_r[0]**2+self.vel_r[1]**2)

    def _next_time_step(self):
        self.time += self.dt

    def reset(self):
        self.log = None
        self.time = 0