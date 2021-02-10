import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np
import scipy


class _All_Dubin_Paths():
    def __init__(self, pos_init=0, pos_final=0, r_traj=100, gamma_g_traj=0,  gamma_traj=0):
        self.pos_init = pos_init
        self.pos_final = pos_final
        self.r_traj = r_traj
        self.gamma_traj = gamma_traj
        self.gamma_g_traj = gamma_g_traj
        self.gammas = atan(abs(self.pos_init[1])/abs(self.pos_init[0]))
        self.mus = self.gammas-self.pos_init[2]
        self.lambdas = self.gammas-self.pos_final[2]
        self.d = sqrt((pos_init[0])**2 + (pos_init[1])**2)

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


    def _RSL(self):
        p_rsl = sqrt(abs(-2*self.r_traj**2 + self.d**2 + cos(self.lambdas-self.mus)*2*self.r_traj**2 + 2*self.r_traj*self.d*(sin(self.lambdas)+ sin(self.mus))))
        q_rsl = - self.mus - atan(2*self.r_traj/p_rsl) + atan((-self.r_traj*cos(self.lambdas)-self.r_traj*cos(self.mus))/(self.d+self.r_traj*sin(self.lambdas)+self.r_traj*sin(self.mus)))
        t_rsl = - self.lambdas - atan(2*self.r_traj/p_rsl) + atan((-self.r_traj*cos(self.lambdas)-self.r_traj*cos(self.mus))/(self.d+self.r_traj*sin(self.lambdas)+self.r_traj*sin(self.mus)))
        
        tau_1 = self.r_traj*tan(self.gamma_traj)*(self.lambdas-self.mus+2*t_rsl) + p_rsl*tan(self.gamma_g_traj)
        if tau_1 > 0:
            self.tau_rsl = tau_1    
            self.rsl_traj = np.array([t_rsl, p_rsl, q_rsl])

    def _LSR(self):
        p_lsr = sqrt(abs(-2*self.r_traj**2 + self.d**2 + cos(self.lambdas-self.mus)*2*self.r_traj**2 - 2*self.r_traj*self.d*(sin(self.lambdas)+ sin(self.mus))))
        q_lsr = self.mus + atan(2*self.r_traj/p_lsr) - atan((self.r_traj*cos(self.lambdas)+self.r_traj*cos(self.mus))/(self.d-self.r_traj*sin(self.lambdas)-self.r_traj*sin(self.mus)))
        t_lsr = self.lambdas + atan(2*self.r_traj/p_lsr) - atan((self.r_traj*cos(self.lambdas)+self.r_traj*cos(self.mus))/(self.d-self.r_traj*sin(self.lambdas)-self.r_traj*sin(self.mus)))

        tau_2 = self.r_traj*tan(self.gamma_traj)*(self.lambdas-self.mus+2*t_lsr) + p_lsr*tan(self.gamma_g_traj)
        if tau_2 > 0:
            self.tau_lsr = tau_2    
            self.lsr_traj = np.array([t_lsr, p_lsr, q_lsr])

    def _RSR(self):
        p_rsr = sqrt(abs(2*self.r_traj**2 + self.d**2 - cos(self.lambdas-self.mus)*2*self.r_traj**2 + 2*self.r_traj*self.d*(-sin(self.lambdas)+ sin(self.mus))))
        q_rsr = self.mus - atan((self.r_traj*cos(self.lambdas)-self.r_traj*cos(self.mus))/(self.d+self.r_traj*sin(self.lambdas)-self.r_traj*sin(self.mus)))
        t_rsr = - self.lambdas  - atan((self.r_traj*cos(self.lambdas)-self.r_traj*cos(self.mus))/(self.d+self.r_traj*sin(self.lambdas)-self.r_traj*sin(self.mus)))

        tau_3 = self.r_traj*tan(self.gamma_traj)*(-self.lambdas+self.mus+2*t_rsr) + p_rsr*tan(self.gamma_g_traj)
        if tau_3 > 0:
            self.tau_rsr = tau_3    
            self.rsr_traj = np.array([t_rsr, p_rsr, q_rsr])

    def _LSL(self):
        p_lsl = sqrt(abs(2*self.r_traj**2 + self.d**2 - cos(self.lambdas-self.mus)*2*self.r_traj**2 + 2*self.r_traj*self.d*(-sin(self.lambdas)+ sin(self.mus))))
        q_lsl = - self.mus + atan((self.r_traj*cos(self.lambdas)-self.r_traj*cos(self.mus))/(self.d-self.r_traj*sin(self.lambdas)+self.r_traj*sin(self.mus)))
        t_lsl = self.lambdas  - atan((self.r_traj*cos(self.lambdas)-self.r_traj*cos(self.mus))/(self.d-self.r_traj*sin(self.lambdas)+self.r_traj*sin(self.mus)))

        tau_4 = self.r_traj*tan(self.gamma_traj)*(-self.lambdas+self.mus) + p_lsl*tan(self.gamma_g_traj)
        if tau_4 > 0:
            self.tau_lsl = tau_4   
            self.lsl_traj = np.array([t_lsl, p_lsl, q_lsl])

    def _LRL(self):
        p_lrl = acos(0.125 * (6 + 2*cos(self.lambdas - self.mus)+ 2*self.d * (sin(self.lambdas) - sin(self.mus))/self.r_traj - self.d**2/self.r_traj**2))
        t_lrl  = self.lambdas + p_lrl/2 - atan((self.r_traj*cos(self.lambdas)-self.r_traj*cos(self.mus))/(self.d-self.r_traj*sin(self.lambdas)+self.r_traj*sin(self.mus)))
        q_lrl  = self.lambdas - t_lrl + p_lrl - self.mus

        tau_5 = self.r_traj * tan(self.gamma_traj) * (self.lambdas-self.mus+p_lrl)
        if tau_5 > 0:
            self.tau_lrl = tau_5
            self.lrl_traj = np.array([t_lrl, p_lrl, q_lrl])

    def _RLR(self):
        p_rlr = acos(0.125 * (6 + 2*cos(self.lambdas - self.mus)- 2*self.d * (sin(self.lambdas) - sin(self.mus))/self.r_traj - self.d**2/self.r_traj**2))
        t_rlr  = -self.lambdas + p_rlr/2 + atan(-(self.r_traj*cos(self.lambdas)+self.r_traj*cos(self.mus))/(self.d+self.r_traj*sin(self.lambdas)-self.r_traj*sin(self.mus)))
        q_rlr  = -self.lambdas - t_rlr + p_rlr

        tau_6 = self.r_traj * tan(self.gamma_traj) * (-self.lambdas+self.mus+2*p_rlr)
        if tau_6 > 0:
            self.tau_rlr = tau_6
            self.rlr_traj = np.array([t_rlr, p_rlr, q_rlr])


    def _Minimum_Tau(self):
        self._RSL()
        self._LSR()
        self._RSR()
        self._LSL()
        self._LRL()
        self._RLR()

        min_array = np.array([self.tau_rsl,self.tau_rsr,self.tau_lsr,self.tau_lsl,self.tau_rlr,self.tau_lrl])
        print(min_array)
        self.tau_min = min(min_array[min_array != 0])

    #def _Calc_Taus():




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