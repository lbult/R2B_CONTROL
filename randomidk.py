# hello
import math
class parafoil:
    def __init__(self, alpha, V, S, b):
        self.alpha = alpha
        self.V = V
        self.S = S
        self.b = b
        self.AR = (b**2)/S
        self.alpha_zl = -7
        self.cl_alpha = self.get_cl_alpha("deg")

    def get_cl_alpha(self,unit):
        a = 0.03714, 1
        b = 0.14286, 5
        tau = (self.AR - a[1])*(b[0]-a[0])/(b[1]-a[1])+a[0]

        cl_alpha_0 = 6.89
        coeff = 2*math.pi*self.AR/cl_alpha_0
        k = coeff*math.tanh(1/coeff)
        cl_alpha_0_corr = cl_alpha_0*k

        cl_alpha = (math.pi*self.AR*cl_alpha_0_corr)/(math.pi*self.AR+cl_alpha_0_corr*(1+tau))

        if unit == "rad":
            return cl_alpha
        if unit == "deg":
            return cl_alpha*math.pi/180


dt = 10**(-5)

a = parafoil(0, 0, 9, 3)
print(a.cl_alpha)

