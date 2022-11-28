import math
import numpy as np

ts = 0.03033

def _Controller(desired_heading, actual_position, 
    desired_pos_x, desired_pos_y, actual_heading, 
    error_var, control_var, desired_heading_var, dynamics_prop):
    '''
    :param left_or_right: negative if left turn, positive if right turn
    :param actual_position: actual position ([0,1]) and heading ([2]), nparray(3)
    :param arc_center: desired arc center of turn, list(2)
    :param desired_radius: desired turn radius, float
    '''
    dx = actual_position[0] - desired_pos_x
    dy = actual_position[1] - desired_pos_y
    
    dphi = actual_heading - desired_heading
    e_along = dx*math.cos(desired_heading) + dy*math.sin(desired_heading)  
    dphi_dt = (desired_heading - desired_heading_var.history[-1])/ts
    
    ''' X-38 PD Gains
    kp = 0.012
    kd = 0.0665'''

    #Drop test PD Gains
    kp = 0.026
    kd = 0.04
    
    #define proportional and derivate error
    e = dx*math.sin(desired_heading) - dy*math.cos(desired_heading)
    de_dt = - dphi_dt * e_along - dynamics_prop.vel_mag * math.cos(dynamics_prop.gamma) * math.sin(dphi)
    
    #define control input
    control_var.update_history(kp * e + kd * de_dt)
    error_var.update_history(e)
    return