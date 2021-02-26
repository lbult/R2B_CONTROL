import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi

#theta = np.linspace(-np.pi,np.pi,100)

# the function, which is y = x^2 here
def answer(theta):
    return np.cos(theta-pi/2) - np.sin(theta-pi/2)*(4-abs(theta))

def second(theta):
    return np.sin(theta-pi/2) + 1 + np.cos(theta-pi/2)*(4-abs(theta))

def one(theta):
    return np.cos(theta-pi/2) - np.sin(theta-pi/2)*(4-abs(theta))

def two(theta):
    return -np.sin(theta-pi/2) - 1 - np.cos(theta-pi/2)*(4-abs(theta))

zeta = 0
f = 0
x =[]
y =[]

k = []
r = []

while f < 1000:
    x.append(answer(zeta))
    y.append(second(zeta))
    k.append(one(zeta))
    r.append(two(zeta))
    zeta += pi/500
    f+=1

plt.plot(x,y, 'r')
plt.plot(k,r, 'r')

# show the plot
plt.show()

"""
Store:
dtheta_r = -arclengths[0]
    theta_r = pi + dtheta_r
    dtheta_l = -arclengths[2]
    theta_l= -pi/2 + dtheta_l

    print(dtheta_l)
    print(dtheta_r)

    def func_1():
        my_one = 0
        interval_one = arclengths[0]/(1000)
        my_theta_r = pi
        while my_one < 1000:
            x_1l.append(init_conditions[0] + r + r*cos(my_theta_r))
            y_1l.append(init_conditions[1] + r*sin(my_theta_r))
            my_theta_r += interval_one
            my_one += 1

    def func_2():
        my_two = 0
        interval_two = arclengths[2]/(1000)
        my_theta_l = -pi/2
        while my_two < 1000:
            x_2l.append(r*cos(my_theta_l ))
            y_2l.append(r + r*sin(my_theta_l))
            my_theta_l += interval_two
            my_two += 1

    def func_3():
        dy = y_2l[-1] - y_1l[-1]
        dx = x_2l[-1] - x_1l[-1]
        a = dy/dx
        b = y_1l[-1] - dy/dx * x_1l[-1]

        x_33 = x_1l[-1]

        my_three = 0
        interval_three = dx/1000
        while my_three < 1000:
            y_3l.append(+dy/dx*x_33+b)
            x_3l.append(x_33)
            x_33 += interval_three
            my_three += 1

    func_1()
    func_2()
    func_3()
"""