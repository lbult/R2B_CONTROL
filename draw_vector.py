import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.text import Annotation
from matplotlib.patches import FancyArrowPatch
import matplotlib.animation as anim
import time

class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._xyz = (x,y,z)
        self._dxdydz = (dx,dy,dz)

    def draw(self, renderer):
        x1,y1,z1 = self._xyz
        dx,dy,dz = self._dxdydz
        x2,y2,z2 = (x1+dx,y1+dy,z1+dz)

        xs, ys, zs = proj_transform((x1,x2),(y1,y2),(z1,z2), renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        super().draw(renderer)
    
def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''
    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)

setattr(Axes3D,'arrow3D',_arrow3D)

mlist = [0.99996915764479, 0.999722430218001, 0.9988898749619708, 0.9969173337331292, 0.9930684569549278, 0.9864292571764973, 0.9759167619387495, 0.9602936856769454, 0.9381913359224867, 0.908143173825084, 0.8686315144381939, 0.8181497174250261, 0.7552818122851862, 0.6788007455329443, 0.5877852522924757, 0.48175367410171754, 0.3608108264876437, 0.22580126686910487]
klist = [-0.007853900888711334, -0.023559764833610164, -0.04710645070964269, -0.07845909572784504, -0.11753739745783782, -0.16418684656886326, -0.21814324139654304, -0.27899110603922994, -0.3461170570774939, -0.41865973753742936, -0.49545866843240904, -0.5750052520432805, -0.6554001709117963, -0.7343225094356884, -0.8090169943749511, -0.8763066800438679, -0.9326390231430992, -0.9741733869698548]

#plt.ion()

def init():
    arrow_one = Arrow3D(0, 0, 0, 0, 0, 0, mutation_scale=20,
            arrowstyle="-|>",
            linestyle='dashed')
    arrow_two = Arrow3D(0, 0, 0, 0, 0, 0, mutation_scale=20,
            arrowstyle="-|>",
            linestyle='dashed')
    ax.add_artist(arrow_one)
    ax.add_artist(arrow_two)

'''def animate(i):
    arrow_one = Arrow3D(0, 0, 0, 
            mlist[i], 0, 0, 
            mutation_scale=20,
            arrowstyle="-|>",
            linestyle='dashed')
    arrow_two = _(0, 0, 0,    
            0, -klist[i], 0, 
            mutation_scale=20,
            arrowstyle="-|>",
            linestyle='dashed')
    ax.add_artist(arrow_one)
    ax.add_artist(arrow_two)
    plt.plot()'''

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('3D Arrows Demo')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_xlim(0,2)

f = 1

while f < len(mlist):
    ax.quiver(0, 0, 0, 1, 1, 1, 
    length=f/2, arrow_length_ratio=0.3, pivot='tail', normalize=False)
    plt.show(block=False)
    #plt.pause(3)
    time.sleep(5)
    plt.close("all")
    f+=1


#ani = anim.FuncAnimation(
    #fig, animate, frames= len(mlist), init_func=init, interval=1, blit=True)




s = 1