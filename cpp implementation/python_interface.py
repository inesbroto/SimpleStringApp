# click_gui.py
import pygame
import subprocess

WIDTH, HEIGHT = 600, 600  # Same size as your 2D plate grid
pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Click to Excite Plate")

running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

        elif event.type == pygame.MOUSEBUTTONDOWN:
            x, y = event.pos
            norm_x = x / WIDTH
            norm_y = y / HEIGHT
            print(f"Clicked at normalized coords: ({norm_x:.2f}, {norm_y:.2f})")

            # Call the C++ executable with coordinates
            #subprocess.run(["./synth", str(norm_x), str(norm_y)])
            subprocess.run(["./plate_interface",str(norm_x), str(norm_y)])

pygame.quit()
