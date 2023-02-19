import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv
import warnings; warnings.simplefilter('ignore')
import imageio

plt.rcParams['font.size'] = 16

# setting physical parameters for advection and diffusion
R0 = 3
m = 3
nu = 1

# setting parameters for numerical solving
Ngrid = 500
Nsteps = 1000
dt = 0.01
dx = 0.1

# initiating spatial array
x = np.arange(0, 50*1.0, dx) / 50
r = x[1:]*R0 

# getting the factors to be used in the solving
#u = -np.mean(3*nu/r)
u = -1
alpha = u*dt/(2*dx)
D = 3*nu
beta = D*dt/dx**2

# setting the initial state of the distribution as a sharp gaussian at 0.5
std = 0.01
f = (1/np.sqrt(2*np.pi)/std) * np.exp(-(x-0.5)**2/(2*std**2))   # array to be modified through solving
f0 = (1/np.sqrt(2*np.pi)/std) * np.exp(-(x-0.5)**2/(2*std**2))  # array to compare initial distribution

# looking at the initial distribution of surface density of the disk
plt.plot(x, f)
plt.title('Initial distribution')
plt.ylabel(r'$\Sigma$')
plt.xlabel('x')
plt.show()

# updating f to show the time evolution of the disk

# initializing the figure
plt.ion()
fig, ax = plt.subplots(1,1)
line, = ax.plot(x, f, '.')

# setting axis limits for consistency
ax.set_xlim([0,1])
ax.set_ylim([0,10])

# update f
for n in range(Nsteps):

    # advection step
    f[1:Ngrid-1] = 0.5*(f[2:] + f[:Ngrid-2]) - alpha*(f[2:] - f[:Ngrid-2])

    # diffusion step
    A = np.eye(Ngrid)*(1.0 + 2.0*beta) + np.eye(Ngrid, k=1)*-beta + np.eye(Ngrid, k=-1)*-beta
    f = np.linalg.solve(A, f)

    # plotting 
    if n % 10 == 0:
        line.set_ydata(f)
        fig.canvas.draw()
        plt.plot(x, f0, color='black', linestyle='dashed')
        plt.xlabel('x')
        plt.ylabel(r'$\Sigma$')
        plt.pause(0.001) 

plt.show()


# making an animation to save
# commented out because saving takes some time
# the animation can be seen by running the file as it is with the plotting above

'''for n in range(Nsteps):

    # advection step
    f[1:Ngrid-1] = 0.5*(f[2:] + f[:Ngrid-2]) - alpha*(f[2:] - f[:Ngrid-2])

    # diffusion step
    A = np.eye(Ngrid)*(1.0 + 2.0*beta) + np.eye(Ngrid, k=1)*-beta + np.eye(Ngrid, k=-1)*-beta
    f = np.linalg.solve(A, f)

    if n % 10 == 0:
        plt.figure(figsize=(10,7))
        plt.plot(x, f0, color='black', linestyle='dashed')
        plt.plot(x, f, '.')
        plt.title('Time evolution of disk')
        plt.ylabel(r'$\Sigma$')
        plt.xlabel('x')
        plt.xlim([0,1])
        plt.ylim([0,10])
        plt.savefig('distribution{}.png'.format(n))
        plt.close()

frames = []
for n in range(Nsteps):
    if n % 10 == 0:
        image = imageio.imread(f'distribution{n}.png')
        frames.append(image)

imageio.mimsave('animation.gif', # output gif
                frames,          # array of input frames
                fps=5)           # optional: frames per second'''