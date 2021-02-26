import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, tan, pi, tanh, atan, asin, acos, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
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
        print(self.a_0)
        self.alpha_0 = alpha_0 #radians
        a = 0.03714, 1
        b = 0.14286, 5
        self.tau = (self.AR - a[1])*(b[0]-a[0])/(b[1]-a[1])+a[0] #correction factor for small aspect ratio wings
        self.a = a_0/(1+(1+self.tau)*(a_0/(pi*self.AR)))
        print(self.a)
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
        self.Right_TE = 0
        self.Left_TE = 0
        self.bank = 0

        # initial condition for angle of attack
        self.alpa = 7.11/180 * pi
        self.alpa_prime = 0

        self.ts = ts  # timestep

        # initialize apparent masses

        self.Parafoil_cloth_mass = 2.1*cloth_density*self.surface
        #self._Apparent_Masses()
        self._Calc_Alpha_Trim(0.02)


    def _Calc_CG_height(self, payload_m):
        self.cg_height_wrtpayload = (self.Parafoil_cloth_mass*self.R )/(self.Parafoil_cloth_mass+ payload_m)
        return self.cg_height_wrtpayload
    
    def _Calc_Alpha_Trim(self, Cds):
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

        print(-Zl*self.Cdl + Zp*Cds)

        def f(x):
            return Cm_0 + Cm_1*x + Cm_2*x**2 +0.09-2.1298*x+22.7498*x**2-649.379*x**3+4915*x**4-10814.1*x**5
        
        def dfdx(x):
            return Cm_1 + 2*Cm_2*x - 2.1298+22.7498*2*x-3*649.379*x**2+4*4915*x**3-5*10814.1*x**4
        
        def Newton_method(x_0):
            iteration = 0
            x_tilde = x_0
            while iteration < 100:
                x_tilde = x_tilde - f(x_tilde)/dfdx(x_tilde)
                iteration += 1
            return x_tilde

        print(Newton_method(7*pi/180)*180/pi)

        x = np.linspace(-pi/38, pi/18, 1000)
        y = np.zeros(1000)
        counter = 0
        for i in x:
            y[counter] = f(i)
            counter+=1

        plt.plot(x,y)
        plt.show()

    def _Calc_Lift(self, alpha, velocity):
        """
        Function for calculating the lift force of the parafoil
        :param alpha: angle of attack of the parafoil, angle between the chord line of parafoil (horizontal axis of body reference system,
         and freestream air (horizontal axis of environment reference system, float
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
        Cd = delta_cd + self.Cd_0 + (1+self.delta) * self.Cl**2 / (pi * self.AR) + self.Cdl 
        
        return -0.5 * density * Cd * self.surface * sqrt(velocity.dot(velocity)) * velocity

    def _Calc_Pitch(self):
        
        return 

    def _Parafoil_Forces_Moments(self, vel):
        """
        :param vel: velocity vector of parafoil, nparray(3)
        :return: Total aerodynamic force vector acting on parafoil, nparray(3)
        """
        
        self.Parafoil_Forces = self._Calc_Lift(self.alpa, vel) + self._Calc_Drag(self.alpa, vel)
        
        return self._Calc_Lift(self.alpa, vel) + self._Calc_Drag(self.alpa, vel)

    def _Parafoil_Control(self, turn_velocity):
        #trailing edge deflection in radians
        if self.Left_TE != 0 and self.Right_TE == 0:
            return np.array([0,0, 0.71 * turn_velocity * self.Left_TE / self.b])
        elif self.Right_TE != 0 and self.Left_TE == 0:
            return np.array([0,0, -0.71 * turn_velocity * self.Right_TE / self.b])
        else:
            return np.array([0,0,0])


class Payload():
    def __init__(self, M=10, shape="", payload_cd=0.02, payload_surface=1):
        self.M=M
        self.shape = shape
        self.payload_cd = payload_cd
        self.payload_surface = payload_surface
        self.Payload_Forces = np.array([])

    def _Calc_Forces(self, velocity):
        self.Payload_Forces = -0.5 * density * self.payload_cd * velocity * sqrt(velocity.dot(velocity))
        return -0.5 * density * self.payload_cd * velocity * sqrt(velocity.dot(velocity)) * np.array([-1,1,1])

class Variable():
    """
    Logger class, for storing and plotting any and all variables
    """
    def __init__(self, var_name=""):
        self.history = []
        self.var_name = var_name
    
    def update_history(self, value):
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
        #plotting more variables
        print("Empty")

class Quaternion():
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

        :param euler: ([float]) euler angles roll, pitch, yaw.
        """
        #self.phi, self.theta, self.psi = euler
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
        #set each individual element
        e0, e1, e2, e3 = self.quaternion
        #update roll, pitch, yaw
        self.phi = np.arctan2(2 * (e0 * e1 + e2 * e3), e0 ** 2 + e3 ** 2 - e1 ** 2 - e2 ** 2)
        self.theta = np.arcsin(2 * (e0 * e2 - e1 * e3))
        self.psi = np.arctan2(2 * (e0 * e3 + e1 * e2), e0 ** 2 + e1 ** 2 - e2 ** 2 - e3 ** 2)
        #print(self.psi)
        self.body_g = np.array([-sin(self.theta), sin(self.phi)*cos(self.theta), cos(self.phi)*cos(self.theta)])

    def _rot_b_v(self):
        """
        Rotate vector from body frame to vehicle frame.

        :param Theta: ([float]) vector to rotate, either as Euler angles or quaternion
        :return: ([float]) rotated vector
        """
        
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

class Dynamics:
    def __init__(self, dt=0.025, mass=10, pos=np.array([0,0,500])):
        
        #translational
        self.pos = pos #reference system
        self.vel_r = np.array([0,0,0]) # reference system
        self.vel = np.array([-17.8,0,6.25]) #body system
        self.acc = np.array([0,0,0]) #body system

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
        self.pos = self.pos + self.vel_r*self.dt #+ 0.5 * (self.dt**2) * np.dot(np.matrix.transpose(rot_bv), self.acc)
        
        #reference system
        self.vel_r = np.dot(np.transpose(rot_bv), self.vel) * np.array([1,1,-1])
        self.vel_mag = np.sqrt(self.vel.dot(self.vel))

        #attitude
        self.gamma = np.arctan2(self.vel_r[2], sqrt(self.vel_r[0]**2+self.vel_r[1]**2))
        self.turn_vel =  sqrt(self.vel_r[0]**2+self.vel_r[1]**2) #self.vel_mag * cos(self.gamma)

    def _next_time_step(self):
        self.time += self.dt

    def reset(self):
        self.log = None
        self.time = 0


class Wind():
    '''
    Class which sets up the wind coordinate system and wind field
    '''
    def __init__(self, w_x=2.0, w_y=2.0, dt=0.025):
        self.wind_field = np.array([w_x, w_y, 0])
        self.time = 0
        self.dt = dt

    def _Wind_Reference_System(self):
        
        return 
    
    '''
    def _Varying(self):
        Add variations with time or altitude?
        Dryden Gust model?

    def _Update_Wind_Field(self):
        """
        Update wind field and position of parafoil system within 
         wind reference system
        """
        self.time += self.dt
    '''