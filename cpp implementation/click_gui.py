# click_gui.py
import pygame
import subprocess
import os

pygame.init()
WIDTH, HEIGHT = 800, 1100

screen = pygame.display.set_mode((WIDTH, HEIGHT))

# Load tapa image (must be in same directory)
image_path = "cajon_tapa.png"
if not os.path.isfile(image_path):
    raise FileNotFoundError(f"Image '{image_path}' not found.")

# Optional: Load a cajón tapa image instead of solid background
tapa_image = pygame.image.load("cajon_tapa.png").convert_alpha()
tapa_image = pygame.transform.scale(tapa_image, (WIDTH, HEIGHT))

#WIDTH, HEIGHT = tapa_image.get_size()

#screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Click to Excite Plate")

# Define light brown color (RGB)
LIGHT_BROWN = (205, 133, 63)
LIGHT_RED = (238, 75, 43)

# Mass-spring normalized coordinates
mass_spring_positions = [
    (0.4, 0.2), (0.2, 0.2), (0.3, 0.2),
    (0.6, 0.2), (0.8, 0.2), (0.7, 0.2),
    (0.4, 0.4), (0.2, 0.4), (0.3, 0.2),
    (0.6, 0.4), (0.8, 0.4), (0.7, 0.2)
]

running = True
while running:
    # Draw background
    screen.fill(LIGHT_BROWN)
    screen.blit(tapa_image, (0, 0))  # Use this if you have an image

    # Draw mass-spring system points
    for norm_x, norm_y in mass_spring_positions:
        x = int(norm_x * WIDTH)
        y = int(norm_y * HEIGHT)
        pygame.draw.circle(screen, LIGHT_RED, (x, y), 6)  # black dots

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

        elif event.type == pygame.MOUSEBUTTONDOWN:
            x, y = event.pos
            norm_x = x / WIDTH
            norm_y = y / HEIGHT
            print(f"Clicked at normalized coords: ({norm_x:.2f}, {norm_y:.2f})")

            subprocess.run(["./plate_interface", str(norm_y), str(norm_x)], capture_output=False)

    pygame.display.flip()

pygame.quit()
