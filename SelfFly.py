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

def _Calc_CG(l1, w1, l2, w2):
    return (l1*w1+l2*w2)/(w1+w2)

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

        #control properties
        self.Right_TE = 0
        self.Left_TE = 0

        # initial condition for angle of attack
        self.alpa = 0
        self.alpa_prime = 0

        self.ts = ts  # timestep

        # initialize apparent masses

        self.Parafoil_cloth_mass = 2.1*cloth_density*self.surface
        self._Apparent_Masses()


    #def alpha_opt(): calculate optimal angle of attack, defining rigging etc based on that
    def _Calc_CG_height(self, payload_m):
        self.cg_height_wrtpayload = (self.Parafoil_cloth_mass*self.R )/(self.Parafoil_cloth_mass+ payload_m)
        return self.cg_height_wrtpayload

    
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
        #calculate total force
        return 0.5 * density * self.Cl * self.surface * ((velocity[2])**2 + (velocity[0])**2) * np.array([0,0,-1])
        #return 0.5 * density * self.Cl * self.surface * np.array([-velocity[2], 0, velocity[0]]) * sqrt(velocity[2]**2 + velocity[0]**2)
 
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
        Cdl = self.line_n*self.R*self.line_d*cos(alpha)**3/self.surface
        Cd = delta_cd + self.Cd_0 + (1+self.delta) * self.Cl**2 / (pi * self.AR) + Cdl 
        # Cd += 0.15
        Cd = delta_cd + self.Cd_0 + (1+self.delta) * self.Cl**2 / (pi * self.AR)
        print("Cd", Cd)
        return 0.5 * density * Cd * self.surface * velocity.dot(velocity) * np.array([-1,0,0])
        #return -0.5 * density * Cd * self.surface * sqrt(velocity.dot(velocity)) * velocity * np.array

    def _Calc_Pitch(self, velocity, gamma, D_s, y, y_p, m_s):
        """
        TODO: add the airfoil pitching moment term in this pitching EOM, and consider how controls can be integrated
        :param velocity: velocity, nparray(3)
        :param gamma: gamma, angle between velocity vector and ground, float
        :param D_s: Drag force of payload, nparray(3)
        :param y: initial angle of attack, float
        :param y_p: initial rate of change of angle of attack, float
        :param m_s: mass of payload, float
        :return: y_1: change in angle of attack after one time step, float
                 y_prime: change in rate of change of angle of attack after one time step, float
        """
        #assume airfoil pitch coefficient negligible for now
        
        XYZ = -np.matmul(self.app_m_curv,(velocity - np.matmul(_Skew_Symmetric_Operator(np.array([0, 0, self._Calc_CG_height(m_s)])), np.array([0, y_p,0]))))
        LMN = -np.matmul(self.app_MMOI_curv,np.array([0,y_p,0]))
        Parafoil_Forces = self.Parafoil_Forces #self._Calc_Lift(y, velocity) + self._Calc_Drag(y, velocity)

        aa = _Skew_Symmetric_Operator(np.array([0,y_p,0]))
        bb = self.app_MMOI_curv + np.array([[0,0,0],[0,0,0],[0,0,m_s * self._Calc_CG_height(m_s) * self._Calc_CG_height(m_s)]])

        cc = np.matmul(aa, bb)
        angular_rates_matrix = np.array([0, y_p, 0])
        dd = np.matmul(cc, angular_rates_matrix)

        DE_RHS = LMN - \
                 dd + \
                 np.matmul(_Skew_Symmetric_Operator(np.array([0,0,self.R - self._Calc_CG_height(m_s)])),XYZ) + \
                 np.matmul(_Skew_Symmetric_Operator(np.array([0,0,self._Calc_CG_height(m_s)])),D_s) + \
                 np.matmul(_Skew_Symmetric_Operator(np.array([0,0,self.R - self._Calc_CG_height(m_s)])),Parafoil_Forces)
        
        # print(DE_RHS, "egeg")

        y_prime_prime_matrix = np.matmul(np.linalg.inv(bb),DE_RHS)
        
        y_prime_prime = y_prime_prime_matrix[1]
        
        y_prime = y_prime_prime * self.ts
        y_1 = y_prime * self.ts
        
        #print(y_1, "egegeweg")
        '''
        #
        vel = sqrt(velocity.dot(velocity))
        D_s = sqrt(D_s.dot(D_s))
        L_l = -density*vel*vel*self.line_n*self.line_d*self.R*np.cos(y+gamma)**2*np.sin(y+gamma)/2
        D_l = density*vel*vel*self.line_n*self.line_d*self.R*np.cos(y+gamma)**3/2
        y_prime_prime = (self.R * (D_s * np.cos(y + self.rigging)) + self.R * (
                 L_l * np.sin(y + self.rigging) - D_l * np.cos(y + self.rigging)) / 2 - m_s * 9.80665 * self.R * np.sin(
                 y - gamma + self.rigging)- 2000000*y_p) / (m_s * self.R * self.R + self.app_MMOI_curv[1][1])
        y_prime = y_prime_prime * self.ts
        y_1 = y_p * self.ts
        '''
        # print(dd, "aegaewaewbwg")
        return y_1, y_prime

    def _Parafoil_Forces_Moments(self, vel, D_s, m_s, vel_r):
        """

        :param vel: velocity vector of parafoil, nparray(3)
        :param D_s: Drag force vector of payload, nparray(3)
        :param m_s: mass of payload, float
        :return: Total aerodynamic force vector acting on parafoil, nparray(3)
        """
        self.gamma = -abs(np.arctan2(vel_r[2], vel_r[0]))
        a_1, a_p_1 = self._Calc_Pitch(vel, self.gamma, D_s, self.alpa, self.alpa_prime, m_s)
        #self.alpa += a_1
        self.alpa_prime += a_p_1
        self.Parafoil_Forces = self._Calc_Lift(self.alpa, vel) + self._Calc_Drag(self.alpa, vel)
        print(self._Calc_Lift(self.alpa, vel))
        print(self._Calc_Drag(self.alpa, vel))
        return self._Calc_Lift(self.alpa, vel) + self._Calc_Drag(self.alpa, vel), a_p_1

        #Cm_pitch = self._Calc_Pitch(0)
        #self.Parafoil_Moments = np.array([0,Cm_pitch,0])

    def _Parafoil_Control(self, turn_velocity):
        #trailing edge deflection in radians
        #see if it is possible to cause an increase in drag if both are deflected
        #span = sqrt(self.AR*self.surface)
        if self.Left_TE != 0 and self.Right_TE == 0:
            return np.array([0,0, 0.71 * turn_velocity * self.Left_TE / self.b])
        elif self.Right_TE != 0 and self.Left_TE == 0:
            return np.array([0,0, -0.71 * turn_velocity * self.Right_TE / self.b])
        else:
            return np.array([0,0,0])

    def _Apparent_Masses(self):
        epsilon_0 = 2*self.anhedral
        flat_b = 2*np.sin(epsilon_0)*self.R
        k = np.array([0.848, -0.9, self.AR/(1+self.AR)]) # array: [kA, kB, kC]
        k_star = np.array([0.84*self.AR/(1+self.AR), 1.161*self.AR/(1+self.AR), 0.848]) # array: [kA*, kB*, kC*]

        app_m_fl = np.array([density*pi*(self.t**2)*flat_b/4* 0.848, density*pi*(self.t**2)*self.c/4 * (-0.9), density*pi*(self.c**2)*flat_b/4*(self.AR/(1+self.AR))])
        # array: [app_mass_x_fl, app_mass_y_fl, app_mass_z_fl]

        app_MMOI_fl = np.array([density*pi*(self.c**2)*(flat_b**3)/48, density*4*flat_b*(self.c**4)/(48*pi), density*pi*(self.t**2)*(flat_b**3)/48])\
                  *k_star
        # array: [app_MMOIx_fl, app_MMOIy_fl, app_MMOIz_fl]

        a_bar = (self.R-self.R*np.cos(epsilon_0))/(2*self.R*np.sin(epsilon_0))
        a1 = self.R*np.sin(epsilon_0)/epsilon_0
        a2 = a1*app_m_fl[1]/(app_m_fl[1]+app_MMOI_fl[0]/(self.R**2))
        a12 = a1 - a2


        self.app_m_curv = np.array([[app_m_fl[0]*(1+8*(a_bar**2)/3),0,0],
                               [0,((self.R**2)*app_m_fl[1]+app_MMOI_fl[0])/(a1**2),0],
                               [0,0,app_m_fl[2]*np.sqrt(1+2*(a_bar**2)*(1-(self.t/self.c)**2))]])
        # array: [[app_mx_curv,0,0], [0,app_my_curv,0], [0,0,app_mz_curv]]


        self.app_MMOI_curv = np.array([[((a12*self.R)**2*app_m_fl[1]+app_MMOI_fl[0]*a2**2)/(a1**2),0,0],
                                 [0,app_MMOI_fl[1]*(1+pi*(1+self.AR)*self.AR*(a_bar*(self.t/self.c))/6),0],
                                  [0,0,(1+8*a_bar**2)*app_MMOI_fl[2]]])
        # array: [[app_MMOIx_curv, 0, 0], [0,app_MMOIy_curv,0], [0,0,app_MMOIz_curv]]

        # print(self.app_m_curv)


class Payload():
    def __init__(self, M=10, shape="", payload_cd=1, payload_surface=1):
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
'''
For now, add yaw as an angular rate
'''

class Dynamics:
    def __init__(self, dt=0.025, I=np.array([[1,0,0],[0,1,0],[0,0,1]]), mass=10, pos=np.array([0,0,300])):
        #properties, I is inertia matrix, mass in Newton
        self.I = I
        self.mass = mass
        #translational
        self.pos = pos #reference system
        self.vel_r = np.array([0,0,0]) # reference system
        self.vel = np.array([10,0,1]) #body system
        self.acc = np.array([0,0,0]) #body system
        #attitude, radians
        self.angular_velocity = np.array([[0,0,0]])
        self.angular_acceleration = np.array([[0,0,0]])
        #position, attitude log
        self.dt = dt
        self.time = 0
        #forces, moments, 3x1 matrices
        self.forces = np.array([0,0,0])
        self.moments = np.array([0,0,0])
        #information for turn
        self.vel_mag = 0
        self.gamma = 0
        self.turn_vel = 0

    def update_dynamics(self, rot_bv):
        #self.forces = np.dot(np.transpose(rot_bv), self.forces#* np.array([1,1,0]) + self.forces * np.array([0,0,1]) + np.array([0,0,self.mass*9.81])
        #translational acceleration
        self.acc = self.forces * 1/self.mass
        self.vel = self.vel + self.acc * self.dt
        self.vel_r = np.dot(np.transpose(rot_bv), self.vel) * np.array([1,1,-1])
        self.vel_mag = np.sqrt(self.vel.dot(self.vel))
        #vel_reference = self.vel *  np.array([1,1,-1]) #switch from right hand to attitude positive upwards
        #translation
        self.pos = self.pos + self.vel_r*self.dt #+ 0.5 * (self.dt**2) * np.dot(np.matrix.transpose(rot_bv), self.acc)
        #attitude
        #self.angular_velocity = np.add(self.angular_velocity, np.dot(self.dt, self.angular_acceleration))
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