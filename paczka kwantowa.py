import pygame
import numpy as np
import time

# Parametry fizyczne 
N = 400
L = 400e-9
dx = L / N
x = np.linspace(0, L, N)
dt = 0.06e-13
hbar = 1.0545718e-34
m_e = 9.10938356e-31

# Parametry pakietu 
sigma = 10e-9
x0 = L / 4
k0 = 5e7

# Funkcja falowa 
real_m = np.exp(-((x - x0) ** 2) / (2 * sigma ** 2)) * np.cos(k0 * x)
imag_m_half = np.exp(-((x - x0) ** 2) / (2 * sigma ** 2)) * np.sin(k0 * x)

V = np.zeros(N)

Y_AXIS_RANGES = {
    "V": (0.0, 1e-22),
    "Psi": (-10000, 10000),
    "Psi_sq": (0.0, 0.6*100000000)
}

pygame.init()
WIDTH, HEIGHT = 1200, 950
screen = pygame.display.set_mode((WIDTH, HEIGHT))
clock = pygame.time.Clock()
font = pygame.font.SysFont("Arial", 14)

marg = 50
button_height = 30
graph_height = (HEIGHT - 6 * marg - 2 * button_height) // 3
left_margin_with_labels = marg + 40


top_graph = pygame.Rect(left_margin_with_labels, marg, WIDTH - left_margin_with_labels - marg, graph_height)
middle_graph = pygame.Rect(left_margin_with_labels, top_graph.bottom + marg, WIDTH - left_margin_with_labels - marg, graph_height)
bottom_graph = pygame.Rect(left_margin_with_labels, middle_graph.bottom + marg, WIDTH - left_margin_with_labels - marg, graph_height)

button_y = bottom_graph.bottom + marg
start_stop_button = pygame.Rect(marg, button_y, 100, button_height)
reset_button = pygame.Rect(start_stop_button.right + 10, button_y, 100, button_height)
reset_wave_button = pygame.Rect(reset_button.right + 10, button_y, 100, button_height)
reset_potential_button = pygame.Rect(reset_wave_button.right + 10, button_y, 140, button_height)
mode_buttons = {
    "manual": pygame.Rect(reset_potential_button.right + 10, button_y, 100, button_height),
    "drop": pygame.Rect(reset_potential_button.right + 120, button_y, 100, button_height),
    "invdrop": pygame.Rect(reset_potential_button.right + 230, button_y, 120, button_height),
}

modes = ["manual", "drop", "invdrop"]
mode_names = {
    "manual": "Tryb: ręczny",
    "drop": "Tryb: spadek",
    "invdrop": "Tryb: odw. spadek"
}
current_mode = "manual"

start_time = time.time()
paused_time = 0.0
paused = False
running = True

start_wave = real_m.copy(), imag_m_half.copy()
start_potential = V.copy()

def draw_axes(rect, label, y_range):
    pygame.draw.rect(screen, (0, 0, 0), rect, 1)
    screen.blit(font.render(label, True, (0, 0, 0)), (rect.left + 5, rect.top + 5))
    for i in range(6):
        x_pos = rect.left + i * rect.width // 5
        pygame.draw.line(screen, (200, 200, 200), (x_pos, rect.top), (x_pos, rect.bottom), 1)
        if rect.bottom + 2 < HEIGHT:
            x_label_val = i * L / 5 * 1e9  
            x_label = font.render(f"{x_label_val:.0f} nm", True, (0, 0, 0))
            screen.blit(x_label, (x_pos - x_label.get_width() // 2, bottom_graph.bottom + 5))
    y_min, y_max = y_range
    for i in range(5):
        y_pos = rect.top + i * rect.height // 4
        pygame.draw.line(screen, (220, 220, 220), (rect.left, y_pos), (rect.right, y_pos), 1)
        val_y = y_max - i * (y_max - y_min) / 4
        label_text = f"{val_y:.2e}"
        label_surf = font.render(label_text, True, (0, 0, 0))
        screen.blit(label_surf, (rect.left - 38, y_pos - 7))

def calculate_scale_y_fixed(y_range, rect):
    y_min, y_max = y_range
    return rect.height / (y_max - y_min)

def draw_curve(rect, values, color, y_range):
    scale_y = calculate_scale_y_fixed(y_range, rect)
    y_min, y_max = y_range
    for i in range(N - 1):
        y1 = rect.bottom - (values[i] - y_min) * scale_y
        y2 = rect.bottom - (values[i + 1] - y_min) * scale_y
        if np.isfinite(y1) and np.isfinite(y2):
            x1 = int(rect.left + i * rect.width / N)
            x2 = int(rect.left + (i + 1) * rect.width / N)
            if rect.top <= y1 <= rect.bottom and rect.top <= y2 <= rect.bottom:
                pygame.draw.line(screen, color, (x1, int(y1)), (x2, int(y2)), 2)

def draw_button(rect, text, active):
    pygame.draw.rect(screen, (180, 180, 180) if active else (220, 220, 220), rect)
    pygame.draw.rect(screen, (0, 0, 0), rect, 2)
    label = font.render(text, True, (0, 0, 0))
    screen.blit(label, (rect.centerx - label.get_width() // 2, rect.centery - label.get_height() // 2))

def reset_wave():
    global real_m, imag_m_half, start_time, paused_time
    real_m = start_wave[0].copy()
    imag_m_half = start_wave[1].copy()
    start_time = time.time()
    paused_time = 0.0

def reset_potential():
    global V
    V = start_potential.copy()

def apply_mode(mode):
    global V
    if mode == "drop":
        V[:N//2] =  0.5*1.60218e-22
        V[N//2:] = 0
    elif mode == "invdrop":
        V[:N//2] = 0
        V[N//2:] = 0.5*1.60218e-22
    elif mode == "manual":
        
        pass

reset_wave()

while running:
    screen.fill((255, 255, 255))

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        if event.type == pygame.MOUSEBUTTONDOWN:
            mx, my = pygame.mouse.get_pos()
            if top_graph.collidepoint(mx, my) and current_mode == "manual":
                i = int((mx - top_graph.left) * N / top_graph.width)
                if 0 <= i < N:
                    width_barrier = 5
                    left = max(0, i - width_barrier)
                    right = min(N, i + width_barrier + 1)
                    V[left:right] = 0.5*1.60218e-22
            elif start_stop_button.collidepoint(mx, my):
                paused = not paused
                if paused:
                    paused_time += time.time() - start_time
                else:
                    start_time = time.time()
            elif reset_button.collidepoint(mx, my):
                reset_wave()
                reset_potential()
            elif reset_wave_button.collidepoint(mx, my):
                reset_wave()
            elif reset_potential_button.collidepoint(mx, my):
                reset_potential()
            for mode, rect in mode_buttons.items(): 
                if rect.collidepoint(mx, my):
                    current_mode = mode
                    apply_mode(mode)

    if not paused:
        A = hbar * dt / (2 * m_e * dx**2)
        B = dt / hbar

        lap_im = np.zeros_like(imag_m_half)
        lap_im[1:-1] = imag_m_half[2:] - 2 * imag_m_half[1:-1] + imag_m_half[:-2]
        real_m1 = real_m - A * lap_im + B * V * imag_m_half

        lap_real = np.zeros_like(real_m1)
        lap_real[1:-1] = real_m1[2:] - 2 * real_m1[1:-1] + real_m1[:-2]
        imag_m3_half = imag_m_half + A * lap_real - B * V * real_m1

        norm_val = np.sqrt(np.sum(real_m1**2 + imag_m3_half**2) * dx)
        if norm_val > 0:
            real_m1 /= norm_val
            imag_m3_half /= norm_val

        real_m = real_m1
        imag_m_half = imag_m3_half

    psi_sq = real_m**2 + imag_m_half**2

    draw_axes(top_graph, "Potencjał V(x) [J]", Y_AXIS_RANGES["V"])
    draw_axes(middle_graph, "Funkcja falowa: Re (niebieski), Im (czerwony) [1]", Y_AXIS_RANGES["Psi"])
    draw_axes(bottom_graph, "Gęstość prawdopodobieństwa |psi|^2 [1/m]", Y_AXIS_RANGES["Psi_sq"])



#Rysowanie pionowych kresek zależnie od trybu:
    barrier_threshold = 1e-30  

    if current_mode == "manual":
        barrier_indices = np.where(V > barrier_threshold)[0]
        if barrier_indices.size > 0:
            left_barrier_pos = barrier_indices[0]
            right_barrier_pos = barrier_indices[-1]
            if left_barrier_pos != right_barrier_pos:
                for graph_rect in [top_graph, middle_graph, bottom_graph]:
                    x_line_left = int(graph_rect.left + left_barrier_pos * graph_rect.width / N)
                    x_line_right = int(graph_rect.left + right_barrier_pos * graph_rect.width / N)
                    pygame.draw.line(screen, (211, 211, 211), (x_line_left, graph_rect.top), (x_line_left, graph_rect.bottom), 2)
                    pygame.draw.line(screen, (211, 211, 211), (x_line_right, graph_rect.top), (x_line_right, graph_rect.bottom), 2)
            else:
                for graph_rect in [top_graph, middle_graph, bottom_graph]:
                    x_line = int(graph_rect.left + left_barrier_pos * graph_rect.width / N)
                    pygame.draw.line(screen, (211, 211, 211), (x_line, graph_rect.top), (x_line, graph_rect.bottom), 2)
    elif current_mode == "drop" or current_mode == "invdrop":
        if abs(V[0] - V[N//2]) > barrier_threshold:
            for graph_rect in [top_graph, middle_graph, bottom_graph]:
                x_line = int(graph_rect.left + (N//2) * graph_rect.width / N)
                pygame.draw.line(screen, (211, 211, 211), (x_line, graph_rect.top), (x_line, graph_rect.bottom), 2)



    draw_curve(bottom_graph, psi_sq, (0, 150, 0), Y_AXIS_RANGES["Psi_sq"])
    draw_curve(top_graph, V, (255, 165, 0), Y_AXIS_RANGES["V"])
    draw_curve(middle_graph, real_m, (0, 0, 255), Y_AXIS_RANGES["Psi"])
    draw_curve(middle_graph, imag_m_half, (200, 0, 0), Y_AXIS_RANGES["Psi"])
    
    # Przyciski
    draw_button(start_stop_button, "Start/Stop", not paused)
    draw_button(reset_button, "Reset wszystko", False)
    draw_button(reset_wave_button, "Reset fali", False)
    draw_button(reset_potential_button, "Reset potencjału", False)
    for mode, rect in mode_buttons.items():
        draw_button(rect, mode_names[mode], current_mode == mode)

    # Czas symulacji
    if paused:
        display_time = paused_time
    else:
        display_time = time.time() - start_time + paused_time

    time_text = font.render(f"Czas symulacji: {display_time:.2f} s", True, (0, 0, 0))
    screen.blit(time_text, (marg, button_y + button_height + 5))

    # Norma funkcji falowej (unormowanie)
    norm = np.sqrt(np.sum(real_m**2 + imag_m_half**2) * dx)
    norm_text = font.render(f"Unormowanie: {norm:.5f}", True, (0, 0, 0))
    screen.blit(norm_text, (marg + 220, button_y + button_height + 5))

    pygame.display.flip()
    clock.tick(60)

pygame.quit()
