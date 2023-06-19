import FluidCube as fc
import pygame
import sys

N = fc.N
SCALE = fc.SCALE
       

def run():
    pygame.init()
    screen = pygame.display.set_mode((N*SCALE, N*SCALE))
    pygame.display.set_caption('FluidCubeGame')
    fluid = fc.FluidCube(0, 2, 0.2)
    lx, ly, x, y = 0, 0, 0, 0
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                sys.exit()
        
        screen.fill(0)
        fluid.FluidCubeStep()
        lx, ly = x, y 
        x, y = pygame.mouse.get_pos()
        fluid.FluidCubeAddDensity(x//SCALE, y//SCALE, 100)
        amtX = float(x - lx) 
        amtY = float(y - ly)
        fluid.FluidCubeAddVelocity(x//SCALE, y//SCALE, amtX, amtY)
        fluid.renderD()
        fluid.renderV()
        fluid.fadeD()
        pygame.display.flip()

run()
