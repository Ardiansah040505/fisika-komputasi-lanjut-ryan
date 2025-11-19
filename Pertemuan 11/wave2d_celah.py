import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameter domain
nx, ny = 200, 200
dx = dy = 1.0
c = 1.0
dt = 0.5
steps = 500

# Parameter simulasi
normal=True  # Jika True gambar grafik tajam (normalisasi). False, tampil normal

# Stabilitas CFL
assert c * dt / dx < 1 / np.sqrt(2), "CFL condition tidak terpenuhi!"

# Inisialisasi grid
u = np.zeros((nx, ny))
u_prev = np.zeros((nx, ny))
u_next = np.zeros((nx, ny))

# Sumber Gaussian negatif di tengah
x0, y0 = nx // 2, ny // 2
sigma = 5.0
gaussian = -np.exp(-(((np.arange(nx)[:, None] - x0)**2 + (np.arange(ny)[None, :] - y0)**2) / (2 * sigma**2)))
u = -gaussian.copy()
u_prev = -gaussian.copy()

# Domain celah: penghalang vertikal di x = nx//4
mask = np.ones((nx, ny), dtype=bool)
barrier_x = nx // 4
gap_width = ny // 16
gap_start = ny // 2 - gap_width // 2
gap_end = ny // 2 + gap_width // 2

# Set penghalang kecuali celah
mask[barrier_x, :] = False
mask[barrier_x, gap_start:gap_end] = True

# Overlay celah untuk visualisasi
celah_overlay = np.zeros((nx, ny))
celah_overlay[barrier_x, :] = 1.0                      # warnai seluruh penghalang
celah_overlay[barrier_x, gap_start:gap_end] = 0.0      # celah tetap transparan

# Fungsi update
def update(frame):
    global u, u_prev, u_next
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            if not mask[i, j]:
                continue
            laplacian = (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1] - 4 * u[i, j]) / dx**2
            u_next[i, j] = 2 * u[i, j] - u_prev[i, j] + (c * dt)**2 * laplacian

    if normal:  # Nilai amplitudo dinormalisasi menjadi 1
        u_norm = u_next / np.max(np.abs(u_next))
        im_wave.set_array(u_norm)
    else:
        im_wave.set_array(u_next)
    im_celah.set_array(celah_overlay)
    u_prev, u = u, u_next.copy()
    return [im_wave, im_celah]

# Visualisasi
fig, ax = plt.subplots()
im_wave = ax.imshow(u, cmap='seismic', vmin=-1, vmax=1, animated=True)
cbar = fig.colorbar(im_wave, ax=ax)
if normal:
    cbar.set_label("Amplitudo (dinormalisasi)")
else:
    cbar.set_label("Amplitudo") 
im_celah = ax.imshow(celah_overlay, cmap='Blues', alpha=0.5, animated=True)

# Animasi
ani = animation.FuncAnimation(fig, update, frames=steps, interval=30, blit=True)
plt.title("Gelombang 2D dengan Celah")
# ani.save("simulasi_gelombang_celah.mp4", writer='ffmpeg', fps=30)
plt.show()