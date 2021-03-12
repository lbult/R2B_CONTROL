import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi

#Reachable One
thetas = np.linspace(0,np.pi*1.3,100)

turn_radius = 15
length = 120

def reach_able(theta):
    return turn_radius*sin(theta) + (length-turn_radius*theta)*cos(theta), turn_radius*cos(0) - turn_radius*cos(theta) + (length-turn_radius*theta)*sin(theta)

xs_1 = []
ys_1 = []
ys_negative_1 = []


for angle in thetas:
    if reach_able(angle)[1] >= 0:
        xs_1.append(reach_able(angle)[0])
        ys_1.append(reach_able(angle)[1])
        ys_negative_1.append(-reach_able(angle)[1])

xs = []
ys = []
ys_negative = []

def reach_able_2(theta):
    x = turn_radius*sin(theta) + (length-turn_radius*theta-turn_radius*pi)*cos(theta) - turn_radius*sin(theta) + turn_radius*sin(theta + pi)
    y = turn_radius*cos(0) - turn_radius*cos(theta) + (length-turn_radius*theta-turn_radius*pi)*sin(theta) + turn_radius*cos(theta) - turn_radius*cos(theta+pi)
    print(length-turn_radius*theta-turn_radius*pi)
    return x,y

for angle in thetas:
    if reach_able_2(angle)[1] >= 0:
        xs.append(reach_able_2(angle)[0])
        ys.append(reach_able_2(angle)[1])
        ys_negative.append(-reach_able_2(angle)[1])

plt.plot(xs, ys, 'b', label="Total Reachable Area")
plt.plot(xs, ys_negative, 'b')
plt.plot(xs_1, ys_1, 'r', label="180 degree turn into the wind")
plt.plot(xs_1, ys_negative_1, 'r')
plt.gca().set_aspect('equal', adjustable='box')

plt.legend(loc='lower center')

plt.show()