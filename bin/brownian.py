#!/usr/bin/env python3
import numpy as np
from pylab import plot, show, grid, axis, xlabel, ylabel, title
from math import sqrt
from scipy.stats import norm
from matplotlib import pyplot as plt
from matplotlib import animation
#########################################################################################
def brownian(x0, n, dt, delta, out=None):
    """
    Generate an instance Brownian motion based on:
    x(t) = x(0) + N(0, delta**2*t; 0, t)
    where N(a,b;t0,t1) is a normally distributed random variable with mean a and variance b.
    """
    x0 = np.asarray(x0)

    r = norm.rvs(size=x0.shape + (n,), scale=delta*sqrt(dt))

    if out is None:
        out = np.empty(r.shape)
    
    #Compute Brownian motion by forming the cumulative sum of the random samples
    np.cumsum(r, axis=-1, out=out)

    #Add the initial condition.
    out += np.expand_dims(x0, axis=-1)
    return out
#######################################################################################
#The Wiener process parameter.
delta = 0.25

T = 10.0

N = 500
#Time step size
dt = T/N

# Initial values for x
x = np.empty((2, N+1))
x[:,0] = 0.0

brownian(x[:,0], N, dt, delta, out=x[:,1:])

plot(x[0], x[1])

plot(x[0,0], x[1,0], 'go')
plot(x[0,-1],x[1,-1], 'ro')

title('2D Brownian Motion')
xlabel('x', fontsize=20)
ylabel('y', fontsize=20)
axis('equal')
grid(True)
show()

# fig = plt.figure()
# ax = plt.axes(xlim=(0,1), ylim=(-1,1))
# line, = ax.plot([],[], lw=2)

# def init():
#     line.set_data([],[])
#     return line

# def animate(i):
#     x=np.array(x[0])
#     y=np.array(x[1])
#     line.set_data(x,y)
#     return line

# anim = animation.FuncAnimation(fig, animate, init_func=
# init, frames=200, interval=20, blit=True)
    
# anim.save('brownian.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
# plot.show()