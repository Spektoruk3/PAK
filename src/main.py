import numpy as np
from math import floor


class FluidCube:
    def __init__(self, size, diffusion, viscosity, dt):
        self.size = size
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity
        n = size

        self.s = np.zeros((n, n))
        self.density = np.zeros((n, n))

        self.Vx = np.zeros((n, n))
        self.Vy = np.zeros((n, n))

        self.Vx0 = np.zeros((n, n))
        self.Vy0 = np.zeros((n, n))


    def FluidCubeAddDensity(self, x, y, amount):
        self.density[x, y] += amount


    def FluidCubeAddVelocity(self, x, y, amountX, amountY):
        self.Vx[x, y] += amountX
        self.Vy[x, y] += amountY

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
    
    Nfloat = n
    

    for j, jfloat in range(1,n - 1): 
        for i, ifloat in range(1,n - 1):
            tmp1 = dtx * velocX[i, j]
            tmp2 = dty * velocY[i, j]
            x    = ifloat - tmp1 
            y    = jfloat - tmp2
               
            if x < 0.5 : x = 0.5 
            if x > Nfloat + 0.5 : x = Nfloat + 0.5 
            i0 = floor(x)
            i1 = i0 + 1.0
            if(y < 0.5): y = 0.5 
            if(y > Nfloat + 0.5): y = Nfloat + 0.5 
            j0 = floor(y)
            j1 = j0 + 1.0 
               
            s1 = x - i0
            s0 = 1.0 - s1 
            t1 = y - j0
            t0 = 1.0 - t1
             
            i0i = i0
            i1i = i1
            j0i = j0
            j1i = j1
                
            d[i, j] = s0 * (t0 * d0[i0i, j0i] + t1 * d0[i0i, j1i]) +s1 * (t0 * d0[i1i, j0i] + t1 * d0[i1i, j1i])
    set_bnd(b, d, n)




arr = np.array([2,3,4,5,6])
print(arr[1:-1])