import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np
import scipy


class _All_Dubin_Paths():
    def __init__(self, pos_init=0, pos_final=0, r_traj=10, gamma_g_traj=0,  gamma_traj=0):
        self.pos_init = pos_init
        self.pos_final = pos_final
        self.r_traj = r_traj
        self.gamma_traj = gamma_traj
        self.gamma_g_traj = gamma_g_traj

        #initiate the cost of all paths
        self.tau_rsl = 0
        self.tau_rsr = 0
        self.tau_lsr = 0
        self.tau_lsl = 0
        self.tau_rlr = 0
        self.tau_lrl = 0
        
        self.tau_min = 0
        
        #initiate all lengths t, p, q
        self.rsl_traj = np.array([0,0,0])
        self.rsr_traj = np.array([0,0,0])
        self.lsr_traj = np.array([0,0,0])
        self.lsl_traj = np.array([0,0,0])
        self.rlr_traj = np.array([0,0,0])
        self.lrl_traj = np.array([0,0,0])

        #trajectory as x, y, heading, altitude lists
        self.pos_x = [0]
        self.pos_y = [0]
        self.heading = [0]
        self.alt = [0]


    def _RSL(self):
        try:
            Lcc = sqrt((self.pos_final[0]-self.r_traj*sin(self.pos_final[2])-self.r_traj)**2+ (self.pos_final[1]+self.r_traj*cos(self.pos_final[2]))**2)
            Ls = sqrt(Lcc**2 - 4*self.r_traj**2)
            phi_1 = -np.arctan2(self.pos_final[1]+self.r_traj*cos(self.pos_final[2]), self.pos_final[0]-self.r_traj*sin(self.pos_final[2])-self.r_traj) + np.arctan2(2*self.r_traj, Ls) + pi/2
            if phi_1 < 0:
                phi_1 += 2*pi
            phi_2 = self.pos_final[2] + phi_1 - pi/2
            if phi_2 < 0:
                phi_2 += 2*pi
            
            self.tau_rsl = abs((abs(phi_1)+abs(phi_2))*self.r_traj + Ls) * tan(self.gamma_traj)    
            self.rsl_traj = np.array([phi_1, Ls, phi_2])
        except:
            print("Math Domain Error")


    def _LSR(self):
        try:
            Lcc = sqrt((self.pos_final[0]-self.r_traj*sin(self.pos_final[2])+self.r_traj)**2+ (self.pos_final[1]-self.r_traj*cos(self.pos_final[2]))**2)
            Ls = sqrt(Lcc**2 - 4*self.r_traj**2)
            phi_1 = -pi/2 + abs(np.arctan2(self.pos_final[1]-self.r_traj*cos(self.pos_final[2]), self.pos_final[0]-self.r_traj*sin(self.pos_final[2]) +self.r_traj)) + np.arctan2(2*self.r_traj, Ls)
            if phi_1 < 0:
                phi_1 += 2*pi
            phi_2 = -self.pos_final[2] + phi_1 + pi/2
            if phi_2 < 0:
                phi_2 += 2*pi

            self.tau_lsr = abs((Ls + self.r_traj*(phi_1+phi_2))*tan(self.gamma_traj))
            self.lsr_traj = np.array([phi_1, Ls, phi_2])

        except:
            print("Math Domain Error")

    def _LSL(self):
        try:
            phi_1 = np.arctan2( self.pos_final[1]+self.r_traj*cos(self.pos_final[2])-self.r_traj, self.pos_final[0]-self.r_traj*sin(self.pos_final[2])) 
            if phi_1 < 0:
                phi_1 += 2*pi
            phi_2 = self.pos_final[2]-phi_1
            if phi_2 < 0:
                phi_2 += 2*pi
            Ls = sqrt((self.pos_final[1]+self.r_traj*cos(self.pos_final[2])-self.r_traj)**2 + (self.pos_final[0]-self.r_traj*sin(self.pos_final[2]))**2)

            self.tau_lsl = abs((Ls + self.r_traj*(phi_1+phi_2))*tan(self.gamma_traj))
            self.lsl_traj = np.array([phi_1, Ls, phi_2])
        
        except:
            print("Math Domain Error")

    def _RSR(self):
        try:
            phi_1 = np.arctan2( self.pos_final[1]-self.r_traj*cos(self.pos_final[2])+self.r_traj, self.pos_final[0]-self.r_traj*sin(self.pos_final[2]))
            phi_1 -= 2*pi
            phi_1 = abs(phi_1)
            
            Ls = sqrt((self.pos_final[1]-self.r_traj*cos(self.pos_final[2])+self.r_traj)**2 + (self.pos_final[0]-self.r_traj*sin(self.pos_final[2]))**2)
            phi_2 = -phi_1 + self.pos_final[2]
            if phi_2 < 0:
                phi_2 += 2*pi

            self.tau_rsr = abs((abs(phi_1)+abs(phi_2))*self.r_traj + Ls) * tan(self.gamma_traj)    
            self.rsr_traj = np.array([phi_1, Ls, phi_2])
        
        except:
            print("Math Domain Error")
        
        

    def _Minimum_Tau(self):
        self._RSR()
        self._LSR()
        self._RSR()
        self._LSL()

        min_array = np.array([self.tau_rsl,self.tau_rsr,self.tau_lsr,self.tau_lsl])
        print(min_array)
        self.tau_min = min(min_array[min_array != 0])

    def _Go_Left(self, rotate):
        x_i = self.pos_x[-1]
        y_i = self.pos_y[-1]
        self.heading.append(self.heading[-1]) 
        dtheta = rotate/100
        i = 0
        while i < 100:
            self.pos_x.append(x_i - self.r_traj*sin(self.heading[-1]) + self.r_traj*sin(self.heading[-1] + dtheta*i))
            self.pos_y.append(y_i + self.r_traj*cos(self.heading[-1]) - self.r_traj*cos(self.heading[-1] + dtheta*i))
            i += 1
        self.heading.append(self.heading[-1]+rotate)

    def _Go_Right(self, rotate):
        x_i = self.pos_x[-1]
        y_i = self.pos_y[-1]
        dtheta = rotate/100 
        i = 0
        while i < 100:
            self.pos_x.append(x_i + self.r_traj*sin(self.heading[-1]) - self.r_traj*sin(self.heading[-1] - dtheta*i))
            self.pos_y.append(y_i - self.r_traj*cos(self.heading[-1]) + self.r_traj*cos(self.heading[-1] - dtheta*i))
            i += 1
        self.heading.append(self.heading[-1]-rotate)

    def _Straight(self, length):
        self.pos_x.append(self.pos_x[-1]+length*cos(self.heading[-1]))
        self.pos_y.append(self.pos_y[-1]+length*sin(self.heading[-1]))
        

    def _Remove_Path(self):
        self.pos_x = [0]
        self.pos_y = [0]
        self.heading = [0]
        self.alt = [0]


'''

def _Gen_Wind_Field(x_speed, y_speed):
    #try later with periodical disturbances at time tau
    return np.array([x_speed, y_speed, 0])


def _Return_Minimum_Value(values):
    return min(values)

#initial values of time changing variables
gamma = 0
sigma = pi/6
V = 0
#tau = init_height
x_dot = 0
y_dot = 0
phi_dot = 0

#assume continuous density over h
density = 1.225 #kg/m^3

dt = 0.01
tau_min = 0.01
sim = True

gamma_min = atan(tan(gamma_g)/cos(pi/6))
V_min = sqrt(V_g**2 * cos(gamma_min) / (cos(gamma_g)*cos(pi/6)))
r_min = V_min**2 * cos(gamma_min) / (9.81* tan(pi/6))
tau_full = abs(2*pi*r_min/(V_min*cos(gamma_min)) * V_min * sin(gamma_min))

print(r_min)

tau_min = _RSL(atan(init_conditions[1]/init_conditions[0]), 
    sqrt((init_conditions[1])**2 + (init_conditions[0])**2),
    r_min, gamma_min, gamma_g, init_conditions[1]/init_conditions[0]
    )[0]

print(tau_min)
print(tau_full)
print(tau_f)

eta = (tau_f - tau_min)/tau_full

print(eta)

if eta > 0 and tau_min>0:
    while sim:
        
        gamma = atan(tan(gamma_g)/cos(sigma))
        V = sqrt(V_g**2 * cos(gamma) / (cos(gamma_g)*cos(sigma)))

        r = V**2 * cos(gamma) / (9.81* tan(sigma))
        
        tau = _RSL(atan(init_conditions[1]/init_conditions[0]), 
        sqrt(abs(init_conditions[1])**2 + abs(init_conditions[0])**2),
        r, gamma, gamma_g, init_conditions[1]/init_conditions[0]
        )[0]

        arclengths = _RSL(atan(init_conditions[1]/init_conditions[0]), 
        sqrt(abs(init_conditions[1])**2 + abs(init_conditions[0])**2),
        r, gamma, gamma_g, init_conditions[1]/init_conditions[0]
        )[1]

        #parameters that change with time
        if tau > 1.01*tau_f:
            sigma += 0.0005
        elif tau < 0.99*tau_f:
            sigma -= 0.0005
        else:
            sim = False 
            print(sigma)
            print(r)

        
        #tau += -V*sin(gamma) * dt
        
        a_tilda = -1/(tan(gamma))
        u = 1/(r*cos(gamma))

        init_conditions = init_conditions + np.array([
            a_tilda * cos(gamma),
            a_tilda * sin(gamma),
            a_tilda * u
        ])

        
        init_conditions = init_conditions + np.array([
            -(cos(init_conditions[2])/tan(gamma) + wind[0]/V*sin(gamma)),
            -(sin(init_conditions[2])/tan(gamma) + wind[1]/V*sin(gamma)),
            -9.81*tan(sigma)/(V**2*sin(gamma))
        ])

        phi_dot = 9.81 * tan(sigma) / V

    #total_length = V * cos(gamma) * init_height #/(V*sin(gamma))

    arclengths = arclengths

    print(arclengths)

    r = 172.8
    
    #x_1l = init_conditions[0] + r*sin(init_conditions[2]+arclengths[0]) - r*sin(init_conditions[2])
    #y_1l = init_conditions[1] - r*cos(init_conditions[2]+arclengths[0]) + r*cos(init_conditions[2])
    #phi_1l = init_conditions[2]+arclengths[0]
    init_conditions = np.array([-10, -10, 0])

    x_1l =  r*sin(init_conditions[2]+arclengths[0]) - r*sin(init_conditions[2])
    y_1l =  r*cos(init_conditions[2]+arclengths[0]) + r*cos(init_conditions[2])
    phi_1l = init_conditions[2]+arclengths[0]

    x_2l = x_1l+arclengths[1]*cos(phi_1l)
    y_2l = y_1l+arclengths[1]*sin(phi_1l)
    phi_2l = phi_1l

    x_3l = x_2l - r*sin(init_conditions[2]-arclengths[2]) + r*sin(init_conditions[2])
    y_3l = y_2l + r*cos(init_conditions[2]-arclengths[2]) - r*cos(init_conditions[2])
    phi_3l = phi_1l - arclengths[2]

    update = [init_conditions[0],x_1l, x_2l, x_3l]

    update_2 = [init_conditions[1],y_1l, y_2l, y_3l]

    print(update, update_2)

    #plt.plot(init_conditions[0], init_conditions[1], 'bo')
    plt.plot(update,update_2, 'r')
    
    #plt.plot(x_2l,y_2l, 'bo')
    #plt.plot(x_3l,y_3l, 'bo')

    # show the plot
    plt.show()'''