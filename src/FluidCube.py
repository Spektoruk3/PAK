import numpy as np
import pygame
from math import floor

N = 128
SCALE = 4

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
        self.density[x, y] += amount
        self.density[x, y] = self.density[x, y] if self.density[x, y] < 255 else 255
        


    def FluidCubeAddVelocity(self, x, y, amountX, amountY):
        self.Vx[x, y] += amountX
        self.Vy[x, y] += amountY

    
    def renderD(self):
        s = pygame.display.get_surface()
        for i in range(N):
            for j in range(N):
                if (self.density[i, j] > 0):
                    x = i * SCALE
                    y = j * SCALE
                    d = self.density[i, j]
                    d = 255 if d > 255 else d
                    pygame.draw.rect(s, (d, d, d), (x, y, SCALE, SCALE))


    def renderV(self):
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
        self.density = np.where(self.density > 1.8, self.density - 0.2, 0)


    def FluidCubeStep(self):
        n          = self.size
        visc     = self.visc
        diff     = self.diff
        dt       = self.dt
        Vx      = self.Vx
        Vy      = self.Vy
        Vx0     = self.Vx0
        Vy0     = self.Vy0
        s       = self.s
        density = self.density

        diffuse(1, Vx0, Vx, visc, dt, 4, n)
        diffuse(2, Vy0, Vy, visc, dt, 4, n)

        project(Vx0, Vy0, Vx, Vy, 4, n)
        
        advect(1, Vx, Vx0, Vx0, Vy0, dt, n)
        advect(2, Vy, Vy0, Vx0, Vy0, dt, n)

        project(Vx, Vy, Vx0, Vy0, 4, n)

        diffuse(0, s, density, diff, dt, 4, n)
        advect(0, density, s, Vx, Vy, dt, n)


def set_bnd(b, x, n):
    x[1:-1, 0 ] = -x[1:-1, 1 ] if b == 2 else x[1:-1, 1 ]
    x[1:-1, -1] = -x[1:-1, -2] if b == 2 else x[1:-1, -2]

    x[ 0, 1:-1] = -x[ 1, 1:-1] if b == 1 else x[ 1, 1:-1]
    x[ -1,1:-1] = -x[-2, 1:-1] if b == 1 else x[-2, 1:-1]

    x[0,  0 ] = 0.5 * (x[1,  0 ] + x[0,  1 ])
    x[0,  -1] = 0.5 * (x[1,  -1] + x[0,  -2])
    x[-1, 0 ] = 0.5 * (x[-2, 0 ] + x[-1, 1 ])
    x[-1, -1] = 0.5 * (x[-2, -1] + x[-1, -2])
   


def lin_solve(b, x, x0, a, c, iter, n):
    cRecip = 1.0 / c
    for k in range(iter):
        x[1:-1, 1:-1] = (x0[1:-1, 1:-1]
                        + a*(    x[2:,   1:-1]
                                +x[0:-2, 1:-1]
                                +x[1:-1,  2: ]
                                +x[1:-1, 0:-2]
                        )) * cRecip
        set_bnd(b, x, n)


def diffuse (b, x, x0, diff, dt, iter, n):
    a = dt * diff * (n - 2) * (n - 2)
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, n)


def project(velocX, velocY, p, div, iter, n):
    div[1:-1, 1:-1] = -0.5 * (
            velocX[ 2:, 1:-1]
           -velocX[0:-2,1:-1]
           +velocY[1:-1, 2: ]
           -velocY[1:-1,0:-2]
            )/n
    p[1:-1, 1:-1] = 0
    set_bnd(0, div, n); 
    set_bnd(0, p, n)
    lin_solve(0, p, div, 1, 6, iter, n)
 
    velocX[1:-1, 1:-1] -= 0.5 * (  p[2:, 1:-1]
                                    -p[0:-2, 1:-1]) * n
    velocY[1:-1, 1:-1] -= 0.5 * (  p[1:-1, 2:]
                                    -p[1:-1, 0:-2]) * n
    set_bnd(1, velocX, n)
    set_bnd(2, velocY, n)


def advect(b, d, d0,  velocX, velocY, dt, n):
    
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
    set_bnd(b, d, n)

