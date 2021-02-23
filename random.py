from math import pi, sqrt

AR = 2.5
Cdl = 0.1373
Cd_s = 0.02
delta = 0.02
Cl0 = 3.5*pi/180 * 3.5413
Cd0 = 0.005
Cl_alpha = 2*pi
Delta = 1/(2*pi) - (1+delta)/(pi*AR)
Gamma = 1/(2*pi) - 2*(1+delta)/(pi*AR)
chord = 16.1
Xw = 3.996559633027523
Zw = 1.406
Zl = 35.15/chord
Zcg = 69.60776676809186
Zp = 70.3/chord

zw = (Zw-Zcg)/chord
xw = (Xw)/chord

Cm_0 = -Zl*Cdl - Zp*Cd_s + Cl0*xw - Cd0*zw
Cm_1 = Cl_alpha*(xw+Cl0*Delta*zw)
Cm_2 = Gamma*zw*Cl_alpha**2

print(Cm_0)
print(Cm_1)
print(Cm_2)


def find_Cm(x):
    f = -297/68000 * x**5 + 1/1024*x**4 - 13/1536*x**3 + 179/6400*x**2 - 1303/24000*x
    return f


def Cm_is_zero():
    return find_Cm(-9.61) + Cm_0 + Cm_1 * -9.61 + Cm_2 * (-9.61)**2

print(180/pi*(-Cm_1 + sqrt(Cm_1**2-4*(Cm_0-0.12)*Cm_2))/(2*Cm_2))

