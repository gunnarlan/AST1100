from matplotlib.pyplot import *
from mpl_toolkits.mplot3d.axes3d import *
import matplotlib.animation as animation
import numpy as np

L=1e-6
T=1e5
k=1.38064852e-23
m_H2=2.01588 #g/mol
Avogadro=6.022140e23
m_molecule=m_H2/(float(Avogadro)*1000)
sigma=np.sqrt(k*T/float(m_molecule))

def update_lines(num, dataLines, lines) :
    for line, data in zip(lines, dataLines) :
        line.set_data(data[0:2, num-1:num])
        line.set_3d_properties(data[2,num-1:num])
    return lines

#fig = figure()
#ax = fig.add_subplot(111, projection='3d')

def plot():
    def update_lines(num, dataLines, lines) :
        for line, data in zip(lines, dataLines) :
            line.set_data(data[0:2, num-1:num])
            line.set_3d_properties(data[2,num-1:num])
        return lines

    # Attach 3D axis to the figure
    fig = figure()
    ax = Axes3D(fig)
    m = 100   # number of frames in the animation
    n = 25 # number of particles you want to animate
    N = 1000 # number of time steps in your data
    dt=1e-9
    pos=np.zeros((n,3,N))
    pos[:,:,0]=np.random.uniform(0, L, size=(n,3))
    v=np.random.normal(0, sigma, size=(n,3))
    for k in range(1,N):
	pos[:,:,k]=pos[:,:,k-1]+v*dt/float(N)
	particles_passed_left=pos[:,:,k] < 0
	pos[particles_passed_left, k]=0
	particles_passed_right=pos[:,:,k]>L
	pos[particles_passed_right,k]=L
	v[particles_passed_left] = -1*v[particles_passed_left]
	v[particles_passed_right] = -1*v[particles_passed_right]
    

    data = pos # this is the positions of the particles
    # to be animated. In this code it should be an array of shape (n, 3, N)

    # ...Fill out the data array here! data[p, :, i] is particle p, time step i.

    # creates animation data for all your different particles
    lines = [i for i in range(n)]
    for i in range(n):
        lines[i] = [ax.plot(data[i][0,0:1],
        data[i][1,0:1], data[i][2,0:1], 'o')[0]]

    # Set the axes properties
    ax.set_xlim3d([0,L])
    ax.set_xlabel('X')

    ax.set_ylim3d([0,L])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0,L])
    ax.set_zlabel('Z')

    ax.set_title('3D random particles')

    # Creating the Animation object
    ani = [i for i in range(n)]
    for i in range(n):
        ani[i] = animation.FuncAnimation(fig,
        update_lines, m, fargs=([data[i]], lines[i]),
        interval=50, blit=False)
    show()

plot()
