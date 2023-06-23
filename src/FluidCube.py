import numpy as np
import pygame
from math import floor

N = 80
SCALE = 10

class FluidCube:
    def __init__(self, diffusion, viscosity, dt):
        self.size = N
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity
        n = N

        self.s = np.zeros((n, n))
        self.density = np.zeros((n, n))

        self.Vx = np.zeros((n, n))
        self.Vy = np.zeros((n, n))

        self.Vx0 = np.zeros((n, n))
        self.Vy0 = np.zeros((n, n))


    def FluidCubeAddDensity(self, x, y, amount):
        '''Adding some dye density'''
        self.density[x, y] += amount
        self.density[x, y] = self.density[x, y] if self.density[x, y] < 255 else 255
        


    def FluidCubeAddVelocity(self, x, y, amountX, amountY):
        '''Add some velocity'''
        self.Vx[x, y] += amountX
        self.Vy[x, y] += amountY

    
    def renderD(self):
        '''Density rendering'''
        s = pygame.display.get_surface()
        for i in range(N):
            for j in range(N):
                if (self.density[i, j] > 0):
                    x = i * SCALE
                    y = j * SCALE
                    d = self.density[i, j]
                    d = 255 if d > 255 else d
                    pygame.draw.rect(s, (0, 0, d), (x, y, SCALE, SCALE))


    def renderV(self):
        '''Velocity rendering'''
        s = pygame.display.get_surface()
        for i in range(N):
            for j in range(N):
                if (self.density[i, j] > 0):
                    x = i * SCALE
                    y = j * SCALE
                    vx = self.Vx[i, j]
                    vy = self.Vy[i, j]
                    if (abs(vx) >= 0.1 and abs(vy) > 0.1):
                        pygame.draw.line(s, (255, 255, 255), (x, y), (x + vx * SCALE, y + vy * SCALE))


    def fadeD(self):
        '''Disappearance of the density or dissolution of the dye'''
        self.density = np.where(self.density > 0.02, self.density - 0.02, 0)


    def FluidCubeStep(self):
        self.density = np.where(self.density < 255, self.density, 255)
        
        self.Vx0 = diffuse(1, self.Vx0, self.Vx, self.visc, self.dt, 4, N)
        self.Vy0 = diffuse(2, self.Vy0, self.Vy, self.visc, self.dt, 4, N)

        self.Vx0, self.Vy0, self.Vx, self.Vy = project(self.Vx0, self.Vy0, self.Vx, self.Vy, 4, N)

        self.Vx = advect(1, self.Vx, self.Vx0, self.Vx0, self.Vy0, self.dt, N)
        self.Vy = advect(2, self.Vy, self.Vy0, self.Vy0, self.Vx0, self.dt, N)

        self.Vx, self.Vy, self.Vx0, self.Vy0 = project(self.Vx, self.Vy, self.Vx0, self.Vy0, 4, N)

        self.s = diffuse(0, self.s, self.density, self.diff, self.dt, 4, N)
        self.density = advect(0, self.density, self.s, self.Vx, self.Vy, self.dt, N)

        


def set_bnd(b, x, n):
    '''This is short for "set bounds", and it's a way to keep
       fluid from leaking out of your box. Walls are added by
       treating the outer layer of cells as the wall. Basically,
       every velocity in the layer next to this outer layer is
       mirrored. So when you have some velocity towards the
       wall in the next-to-outer layer, the wall gets a 
       velocity that perfectly counters it.'''
    x[1:-1, 0 ] = -x[1:-1, 1 ] if b == 2 else x[1:-1, 1 ]
    x[1:-1, -1] = -x[1:-1, -2] if b == 2 else x[1:-1, -2]

    x[ 0, 1:-1] = -x[ 1, 1:-1] if b == 1 else x[ 1, 1:-1]
    x[ -1,1:-1] = -x[-2, 1:-1] if b == 1 else x[-2, 1:-1]

    x[0,  0 ] = 0.5 * (x[1,  0 ] + x[0,  1 ])
    x[0,  -1] = 0.5 * (x[1,  -1] + x[0,  -2])
    x[-1, 0 ] = 0.5 * (x[-2, 0 ] + x[-1, 1 ])
    x[-1, -1] = 0.5 * (x[-2, -1] + x[-1, -2])
    return x
   


def lin_solve(b, x, x0, a, c, iter, n):
    '''It's solving a linear differential equation of some sort.
    This is done by running through the whole array and setting 
    each cell to a combination of its neighbors. It does this 
    several times; the more iterations it does, the more accurate 
    the results, but the slower things run. After each iteration, 
    it resets the boundaries so the calculations don't explode.'''
    cRecip = 1.0 / c
    for k in range(iter):
        x[1:-1, 1:-1] = (x0[1:-1, 1:-1]
                        + a*(    x[2:,   1:-1]
                                +x[0:-2, 1:-1]
                                +x[1:-1,  2: ]
                                +x[1:-1, 0:-2]
                        )) * cRecip
        x = set_bnd(b, x, n)
    return x


def diffuse (b, x, x0, diff, dt, iter, n):
    '''We use diffusion both in the obvious case
       of making the dye spread out, and also in 
       the less obvious case of making the velocities
       of the fluid spread out'''
    a = dt * diff * (n - 2) * (n - 2)
    x = lin_solve(b, x, x0, a, 1 + 6 * a, iter, n)
    return x


def project(velocX, velocY, p, div, iter, n):
    '''The amount of fluid in each box has to stay constant.
       The amount of fluid going in has to be exactly equal
       to the amount of fluid going out. The other operations
       tend to screw things up so that you get some boxes
       with a net outflow, and some with a net inflow. 
       This operation runs through all the cells and fixes
       them up so everything is in equilibrium.'''
    div[1:-1, 1:-1] = -0.5 * (
            velocX[ 2:, 1:-1]
           -velocX[0:-2,1:-1]
           +velocY[1:-1, 2: ]
           -velocY[1:-1,0:-2]
            )/n
    p[1:-1, 1:-1] = 0
    div = set_bnd(0, div, n); 
    p = set_bnd(0, p, n)
    p = lin_solve(0, p, div, 1, 6, iter, n)
 
    velocX[1:-1, 1:-1] -= 0.5 * (  p[2:, 1:-1]
                                    -p[0:-2, 1:-1]) * n
    velocY[1:-1, 1:-1] -= 0.5 * (  p[1:-1, 2:]
                                    -p[1:-1, 0:-2]) * n
    velocX = set_bnd(1, velocX, n)
    velocY = set_bnd(2, velocY, n)
    return velocX, velocY, p, div


def advect(b, d, d0,  velocX, velocY, dt, n):
    '''Every cell has a set of velocities, and these
       velocities make things move. This is called
       advection. As with diffusion, advection applies
       both to the dye and to the velocities themselves.
       This function is responsible for actually moving 
       things around. To that end, it looks at each cell
       in turn. In that cell, it grabs the velocity, 
       follows that velocity back in time, and sees where
       it lands. It then takes a weighted average of the
       cells around the spot where it lands, then applies
       that value to the current cell.'''
    
    dtx = dt * (n - 2)
    dty = dt * (n - 2)
    
    tmp1 = dtx * velocX
    tmp2 = dty * velocY
    i = np.arange(0, n)
    j = np.arange(0, n)
    j, i = np.meshgrid(i, j)

    x = i - tmp1
    y = j - tmp2
    x[x < .5] = .5
    x[x > n + .5] = n + .5
    y[y < .5] = .5
    y[y > n + .5] = n + .5
    i0 = np.floor(x).astype(np.int32)
    i1 = i0 + 1
    j0 = np.floor(y).astype(np.int32)
    j1 = j0 + 1
    i0 = np.clip(i0, 0 , n-1)
    i1 = np.clip(i1, 0 , n-1)
    j0 = np.clip(j0, 0 , n-1)
    j1 = np.clip(j1, 0 , n-1)
    
    s1 = x - i0
    s0 = 1 - s1
    t1 = y - j0
    t0 = 1 - t1
    d = (s0 *  (t0 * d0[i0, j0]
                            +t1 * d0[i0, j1])
                    +s1 * ( t0 * d0[i1, j0]
                            +t1 * d0[i1, j1]))
    d = set_bnd(b, d, n)
    return d
