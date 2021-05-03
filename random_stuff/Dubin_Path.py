import matplotlib.pyplot as plt
from math import tan, cos, pi, sin, atan, asin, acos,sqrt
import numpy as np
import scipy, timeit


class _All_Dubin_Paths():
    def __init__(self, pos_init=0, pos_final=0, gamma_g_traj=0, altitude=0, v_g=0, sigma_max=0, wind_direction=0, wind_magnitude=1, monte_carlo=100):
        """

        :param pos_init:
        :param pos_final:
        :param gamma_g_traj:
        :param altitude:
        :param v_g:
        :param sigma_max:
        :param wind_direction: measured clockwise from true north starting at 0 degrees up until 359 degrees, 0 degrees = wind
        blowing north, 90 degrees means wind is blowing east etc.
        :param wind_magnitude: speed of wind in meters per second at a height of 10 meters above ground (10 meters above ground
        is the height at which weather stations take their wind speed measurements usually)
        :param monte_carlo: number of iterations to do in a monte carlo loop
        """
        self.pos_init = pos_init
        self.pos_final = pos_final
        self.pos_final_o = np.array([-self.pos_final[1], self.pos_final[0], self.pos_final[2] - pi/2])

        self.gamma_g_traj = gamma_g_traj
        self.altitude = altitude
        self.alt_init = altitude
        self.sigma_max = sigma_max
        self.v_g = v_g
        self.sigma_max_init = sigma_max

        self.gamma_traj = np.arctan2(tan(self.gamma_g_traj), cos(self.sigma_max))
        self.v_min = sqrt(self.v_g**2 * cos(self.gamma_traj) / (cos(self.gamma_traj)*cos(self.sigma_max)))
        self.r_traj = (self.v_min)**2 * cos(self.gamma_traj) / (9.81* tan(self.sigma_max))

        #wind field characteristics & monte carlo analysis
        weibull_shape = 2.2 ## Setting the wind speed distribution parameters, a weibull distribution is the normal distribution but skewed to the lower end
        self.monte_loops = monte_carlo
        heading_sigma = 9   ## standard deviation of wind heading angle in degrees
        wind_heading_distribution = np.random.normal(wind_direction, heading_sigma, self.monte_loops)
        self.wind_vector_distribution = np.array([np.sin(np.pi*wind_heading_distribution/180), np.cos(np.pi*wind_heading_distribution/180)])            ## montecarlo, 2D array of wind headings, with monte_loops number of columns
        self.wind_vector_distribution[0][0], self.wind_vector_distribution[1][0] = np.sin(np.pi*wind_direction/180), np.cos(np.pi*wind_direction/180)   ## simply setting the first wind heading to be the avg heading

        self.wind_speed_distribution = wind_magnitude * np.random.weibull(weibull_shape, self.monte_loops)  ## montecarlo, 2D array of wind velocity magnitudes
        self.wind_speed_distribution[0] = wind_magnitude    ## simply setting the first wind velocity to be the avg heading

        #plotting all the montecarlo distribution stuff to give user info, remove if its annoying
        count, bins, ignored = plt.hist(self.wind_speed_distribution, 30, density=True)
        plt.ylabel("frequency (i think?)")
        plt.xlabel("speed [m/s]")
        plt.title("Wind speed distribution")
        plt.show()
        plt.cla()
        count, bins, ignored = plt.hist(wind_heading_distribution, 30, density=True)
        plt.ylabel("frequency (i think?)")
        plt.xlabel("heading [degree]")
        plt.title("Wind vector distribution")
        plt.show()

        #initiate the cost of all paths
        self.tau_rsl = 0
        self.tau_rsr = 0
        self.tau_lsr = 0
        self.tau_lsl = 0
        self.tau_rlr = 0
        self.tau_lrl = 0
        
        self.tau_min = 0
        self.eta = 0
        
        #initiate all lengths t, p, q
        self.rsl_traj = np.array([0,0,0])
        self.rsr_traj = np.array([0,0,0])
        self.lsr_traj = np.array([0,0,0])
        self.lsl_traj = np.array([0,0,0])
        self.rlr_traj = np.array([0,0,0])
        self.lrl_traj = np.array([0,0,0])

        self.chosen_traj = np.array([0,0,0])

        # shortest trajectory
        # trajectory as x, y, heading, altitude lists
        self.pos_xs = [0]
        self.pos_ys = [0]
        self.headings = [0]
        self.alts = [altitude]

        # minimum control trajectory
        # trajectory as x, y, heading, altitude lists
        self.pos_x = [0]
        self.pos_y = [0]
        self.heading = [0]
        self.alt = [altitude]
        self.control = [0]

    def _RSL(self):
        try:
            Lcc = sqrt((self.pos_final_o[0]-self.r_traj*sin(self.pos_final_o[2])-self.r_traj)**2+ (self.pos_final_o[1]+self.r_traj*cos(self.pos_final_o[2]))**2)
            Ls = sqrt(Lcc**2 - 4*self.r_traj**2)
            phi_1 = -np.arctan2(self.pos_final_o[1]+self.r_traj*cos(self.pos_final_o[2]), self.pos_final_o[0]-self.r_traj*sin(self.pos_final_o[2])-self.r_traj) + np.arctan2(2*self.r_traj, Ls) + pi/2
            if phi_1 < 0:
                phi_1 += 2*pi
            phi_2 = self.pos_final_o[2] + phi_1 - pi/2
            if phi_2 < 0:
                phi_2 += 2*pi
            
            self.tau_rsl = abs((abs(phi_1)+abs(phi_2))*self.r_traj*tan(self.gamma_traj)) + abs( Ls * tan(self.gamma_g_traj))    
            self.rsl_traj = np.array([phi_1, Ls, phi_2])
        except:
            print("Math Domain Error")


    def _LSR(self):
        try:
            Lcc = sqrt((self.pos_final_o[0]-self.r_traj*sin(self.pos_final_o[2])+self.r_traj)**2+ (self.pos_final_o[1]-self.r_traj*cos(self.pos_final_o[2]))**2)
            Ls = sqrt(Lcc**2 - 4*self.r_traj**2)
            phi_1 = -pi/2 + abs(np.arctan2(self.pos_final_o[1]-self.r_traj*cos(self.pos_final_o[2]), self.pos_final[0]-self.r_traj*sin(self.pos_final_o[2]) +self.r_traj)) + np.arctan2(2*self.r_traj, Ls)
            if phi_1 < 0:
                phi_1 += 2*pi
            phi_2 = -self.pos_final_o[2] + phi_1 + pi/2
            if phi_2 < 0:
                phi_2 += 2*pi

            self.tau_lsr = abs((abs(phi_1)+abs(phi_2))*self.r_traj*tan(self.gamma_traj)) + abs( Ls * tan(self.gamma_g_traj))    
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

            self.tau_lsl = abs((abs(phi_1)+abs(phi_2))*self.r_traj*tan(self.gamma_traj)) + abs( Ls * tan(self.gamma_g_traj))    
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

            self.tau_rsr = abs((abs(phi_1)+abs(phi_2))*self.r_traj*tan(self.gamma_traj)) + abs( Ls * tan(self.gamma_g_traj))      
            self.rsr_traj = np.array([phi_1, Ls, phi_2])
        
        except:
            print("Math Domain Error")
        
        

    def _Minimum_Tau(self):
        self._RSR()
        self._LSR()
        self._RSR()
        self._LSL()
        min_array = [self.tau_rsl,self.tau_rsr,self.tau_lsr,self.tau_lsl]
        min_arrays = []
        for taus in min_array:
            if taus != 0:
                min_arrays.append(taus)
        
        self.tau_min = min(min_arrays)
        tau_place = min_array.index(self.tau_min)
        tau_full = abs(2*pi*self.r_traj*tan(self.gamma_traj))

        
        self.eta = (self.altitude - self.tau_min)/tau_full

        iteration_1 = True
        not_converged = True
        tau_place = 0

        while not_converged and self.eta > 0:
            
            if self.sigma_max> 0:
                self.sigma_max -= 0.005
                self._Remove_Path()
            elif self.sigma_max<0:
                tau_place+=1
                self.sigma_max = self.sigma_max_init
                self._Remove_Path()
            elif tau_place >= 4:
                print("Oof")
                not_converged = False


            self.gamma_traj = np.arctan2(tan(self.gamma_g_traj), cos(self.sigma_max))
            self.v_min = sqrt(self.v_g**2 * cos(self.gamma_traj) / (cos(self.gamma_traj)*cos(self.sigma_max)))
            self.r_traj = (self.v_min)**2 * cos(self.gamma_traj) / (9.81* tan(self.sigma_max))

            if iteration_1:
                if tau_place == 0:
                    self._RSL() # i dont think that these dubin paths need to be reinitialized everytime right?
                    #self.pos_final = np.array([self.pos_final[1], -self.pos_final[0], self.pos_final[2] - pi/2])
                    self._Go_Right(self.rsl_traj[0])
                    self._Straight(self.rsl_traj[1])
                    self._Go_Left(self.rsl_traj[2])
                elif tau_place == 1:
                    self._RSR()
                    self._Go_Right(self.rsr_traj[0])
                    self._Straight(self.rsr_traj[1])
                    self._Go_Right(self.rsr_traj[2])
                elif tau_place == 2:
                    self._LSR()
                    #self.pos_final = np.array([-self.pos_final[1], self.pos_final[0], self.pos_final[2] - pi/2])
                    self._Go_Left(self.lsr_traj[0])
                    self._Straight(self.lsr_traj[1])
                    self._Go_Right(self.lsr_traj[2])
                elif tau_place == 3:
                    self._LSL()
                    self._Go_Left(self.lsl_traj[0])
                    self._Straight(self.lsl_traj[1])
                    self._Go_Left(self.lsl_traj[2])

                self.pos_xs = self.pos_x
                self.pos_ys = self.pos_y
                self.headings = self.heading
                self.alts = self.alt
                
                self._Remove_Path()
                iteration_1 = False

            if tau_place == 0:
                self._RSL()
                self._Go_Right(self.rsl_traj[0])
                self._Straight(self.rsl_traj[1])
                self._Go_Left(self.rsl_traj[2])
                if self.tau_rsl < 1.02*self.altitude and self.tau_rsl > 0.98*self.altitude:
                    self.chosen_traj = self.rsl_traj * np.array([1*self.r_traj, 1, -1*self.r_traj]) * abs(tan(self.gamma_g_traj))
                    not_converged = False
                else:
                    self._Remove_Path()
            elif tau_place == 1:
                self._RSR()
                self._Go_Right(self.rsr_traj[0])
                self._Straight(self.rsr_traj[1])
                self._Go_Right(self.rsr_traj[2])
                if self.tau_rsr < 1.02*self.altitude and self.tau_rsr > 0.98*self.altitude:
                    self.chosen_traj = self.rsr_traj * np.array([1*self.r_traj, 1, 1*self.r_traj]) * abs(tan(self.gamma_g_traj))
                    not_converged = False
                else:
                    self._Remove_Path()
            elif tau_place == 2:
                self._LSR()
                self._Go_Left(self.lsr_traj[0])
                self._Straight(self.lsr_traj[1])
                self._Go_Right(self.lsr_traj[2])
                if self.tau_lsr < 1.02*self.altitude and self.tau_lsr > 0.98*self.altitude:
                    self.chosen_traj = self.lsr_traj * np.array([-1*self.r_traj, 1, 1*self.r_traj]) * abs(tan(self.gamma_g_traj))
                    not_converged = False
                else:
                    self._Remove_Path()
            elif tau_place == 3:
                self._LSL()
                self._Go_Left(self.lsl_traj[0])
                self._Straight(self.lsl_traj[1])
                self._Go_Left(self.lsl_traj[2])
                if self.tau_lsl < 1.02*self.altitude and self.tau_lsl > 0.98*self.altitude:
                    self.chosen_traj = self.lsl_traj * np.array([-1*self.r_traj, 1, -1*self.r_traj]) * abs(tan(self.gamma_g_traj))
                    not_converged = False
                else:
                    self._Remove_Path()
            else:
                print("idk")
                not_converged = False
        self.pos_xs_w, self.pos_ys_w = self._Wind_coordinate_Transform(self.pos_xs, self.pos_ys, self.alt)
        self.pos_x_w, self.pos_y_w = self._Wind_coordinate_Transform(self.pos_x, self.pos_y, self.alt)

    def _Go_Left(self, rotate):
        x_i = self.pos_x[-1]
        y_i = self.pos_y[-1]
        this_heading = self.heading[-1]
        dtheta = rotate/100
        i = 0
        while i < 100:
            self.pos_x.append(x_i - self.r_traj*sin(this_heading) + self.r_traj*sin(this_heading + dtheta*i))
            self.pos_y.append(y_i + self.r_traj*cos(this_heading) - self.r_traj*cos(this_heading + dtheta*i))
            self.heading.append(self.heading[-1]+dtheta)
            self.alt.append(self.alt[-1]-abs(dtheta*self.r_traj*tan(self.gamma_traj)))
            self.control.append(-1)
            i += 1

    def _Go_Right(self, rotate):
        x_i = self.pos_x[-1]
        y_i = self.pos_y[-1]
        dtheta = rotate/100 
        this_heading = self.heading[-1]
        i = 0
        while i < 100:
            self.pos_x.append(x_i + self.r_traj*sin(this_heading) - self.r_traj*sin(this_heading - dtheta*i))
            self.pos_y.append(y_i - self.r_traj*cos(this_heading) + self.r_traj*cos(this_heading - dtheta*i))
            self.heading.append(self.heading[-1]-dtheta)
            self.alt.append(self.alt[-1]-abs(dtheta*self.r_traj*tan(self.gamma_traj)))
            self.control.append(1)
            i += 1

    def _Straight(self, length):
        #self.pos_x.append(self.pos_x[-1]+length*cos(self.heading[-1]))
        #self.pos_y.append(self.pos_y[-1]+length*sin(self.heading[-1]))
        dlength = length/100
        i = 0
        while i < 100:
            self.pos_x.append(self.pos_x[-1]+dlength*cos(self.heading[-1]))
            self.pos_y.append(self.pos_y[-1]+dlength*sin(self.heading[-1]))
            self.heading.append(self.heading[-1])
            self.alt.append(self.alt[-1]-abs(dlength*tan(self.gamma_g_traj)))
            self.control.append(0)
            i += 1

    def _Remove_Path(self):
        self.pos_x = [0]
        self.pos_y = [0]
        self.heading = [0]
        self.alt = [self.alt_init]
        self.control = [0]

    def _Wind_Vector_Field(self, z, wind, wind_vector):
        """

        :param z:
        :param wind:
        :param wind_vector:
        :return: i've implemented a form of log rule to make the wind velocity as a function of altitude, plot the return formula on desmos to see how the distribution looks like
        """
        k = 0.3     ## exponent factor that changes the uniformity of wind speeds at higher altitudes
        if z <= 0:
            return (wind/np.log((10**k)/10+1)*np.log((0**k)/10+1)*wind_vector)
        else:
            return (wind / np.log((10 ** k) / 10 + 1) * np.log((z ** k) / 10 + 1) * wind_vector)

    def _Wind_coordinate_Transform(self, x_list, y_list, alt_list):

        x = np.flip(np.array(x_list))
        y = np.flip(np.array(y_list))
        alt = np.flip(np.array(alt_list))
        kappa_g = -1/(self.v_g*sin(self.gamma_g_traj))

        x_w = np.zeros((len(x), self.monte_loops))
        y_w = np.zeros((len(y), self.monte_loops))

        tic = timeit.default_timer()

        # performing the wind coordinate transform
        print(f"Calculating monte carlo simulations for {len(self.wind_speed_distribution)} rounds of iterations...")
        w_vector_k = self.wind_vector_distribution
        w_speed_k = self.wind_speed_distribution
        for i in range(len(x)):
            x_temp = x[i]
            y_temp = y[i]
            for j in range(i, len(x) - 1):
                dtau = alt[j] - alt[j + 1]
                wind = self._Wind_Vector_Field(alt[j], w_speed_k, w_vector_k)
                # simpson's rule of numerical integration
                if j == i or j == len(x) - 2:
                    x_temp -= kappa_g * wind[0] * dtau / 3
                    y_temp -= kappa_g * wind[1] * dtau / 3
                elif j % 2 == 0:
                    x_temp -= 4 * kappa_g * wind[0] * dtau / 3
                    y_temp -= 4 * kappa_g * wind[1] * dtau / 3
                elif j % 2 == 1:
                    x_temp -= 2 * kappa_g * wind[0] * dtau / 3
                    y_temp -= 2 * kappa_g * wind[1] * dtau / 3

            x_w[i] = x_temp
            y_w[i] = y_temp


        toc = timeit.default_timer()
        print(f'monte carlo time = {toc-tic}')
        return np.flip(x_w).transpose(), np.flip(y_w).transpose()