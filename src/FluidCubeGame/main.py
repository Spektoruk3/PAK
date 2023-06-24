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
    fluid = fc.FluidCube(1e-9 / 2, 1e-7, 1e-2)
    clock = pygame.time.Clock()
    lx, ly, x, y = 0, 0, 0, 0
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                sys.exit()
        
        clock.tick(100)
        screen.fill(0)
        fluid.FluidCubeStep()
        lx, ly = x, y 
        x, y = pygame.mouse.get_pos()
        fluid.FluidCubeAddDensity(x//SCALE, y//SCALE, 300)
        amtX = float(x - lx) 
        amtY = float(y - ly) 
        fluid.FluidCubeAddVelocity(x//SCALE, y//SCALE, amtX, amtY)
        fluid.fadeD()
        fluid.renderD()
        pygame.display.flip()

def main():
    run()

if __name__ == '__main__':
    main()