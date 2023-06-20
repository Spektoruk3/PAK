import FluidCube as fc
import pygame
import sys
import cv2

N = fc.N
SCALE = fc.SCALE
       

def run():
    pygame.init()
    screen = pygame.display.set_mode((N*SCALE, N*SCALE))
    pygame.display.set_caption('FluidCubeGame')
    fluid = fc.FluidCube(1e-2, 1e-9 / 2, 1e-7)
    surf = pygame.surfarray.make_surface(cv2.resize(fluid.density, (N*SCALE, N*SCALE)))

    lx, ly, x, y = 0, 0, 0, 0
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                sys.exit()
        
        screen.fill(0)
        fluid.FluidCubeStep()
        lx, ly = x, y 
        x, y = pygame.mouse.get_pos()
        fluid.FluidCubeAddDensity(x//SCALE, y//SCALE, 250)
        amtX = float(x - lx) 
        amtY = float(y - ly)
        fluid.FluidCubeAddVelocity(x//SCALE, y//SCALE, amtX, amtY)
        fluid.fadeD()
        pygame.surfarray.blit_array(surf, cv2.resize(fluid.density, (N*SCALE, N*SCALE)))
        screen.blit(surf, (0, 0))
        pygame.display.flip()

run()