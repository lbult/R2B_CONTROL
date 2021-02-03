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
    zeta += pi/1000
    f+=1

plt.plot(x,y, 'r')
plt.plot(k,r, 'r')

# show the plot
plt.show()