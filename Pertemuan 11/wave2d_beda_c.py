import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameter domain
nx, ny = 200, 200
dx = dy = 1.0
dt = 0.5
steps = 500

# Inisialisasi grid
u = np.zeros((nx, ny))
u_prev = np.zeros((nx, ny))
u_next = np.zeros((nx, ny))

# Sumber Gaussian negatif di tengah
x0, y0 = nx // 2, ny // 2
sigma = 5.0
gaussian = -np.exp(-(((np.arange(nx)[:, None] - x0)**2 + (np.arange(ny)[None, :] - y0)**2) / (2 * sigma**2)))
u = gaussian.copy()
u_prev = gaussian.copy()

# Peta kecepatan gelombang: c = 1.0 di bawah, c = 0.5 di atas
c_map = np.ones((nx, ny))
barrier_x = nx
gap_width = ny // 8
gap_start = 150
gap_end = 200
c_map[0:barrier_x, gap_start:gap_end] = 0.5  # material lambat di atas, sejajar celah

# Fungsi update
def update(frame):
    global u, u_prev, u_next
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            laplacian = (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1] - 4 * u[i, j]) / dx**2
            u_next[i, j] = 2 * u[i, j] - u_prev[i, j] + (c_map[i, j] * dt)**2 * laplacian

    # Normalisasi amplitudo
    u_norm = u_next / np.max(np.abs(u_next))
    im.set_array(u_norm)
    u_prev, u = u, u_next.copy()
    return [im]

# # Visualisasi
# fig, ax = plt.subplots()
# im = ax.imshow(u, cmap='seismic', vmin=-1, vmax=1, animated=True)
# ax.axvline(gap_start, color='red', linestyle='--', linewidth=2, label='Batas Material')
# ax.legend(loc='upper right')

# cbar = fig.colorbar(im, ax=ax)
# cbar.set_label("Amplitudo (dinormalisasi)")
# ani = animation.FuncAnimation(fig, update, frames=steps, interval=30, blit=True)
# plt.title("Gelombang 2D dengan Dua Jenis Material")
# plt.show()

# Visualisasi
fig, ax = plt.subplots()
im = ax.imshow(u, cmap='seismic', vmin=-1, vmax=1, animated=True)

# Garis vertikal di batas material (misalnya di kolom gap_start)
garis_batas = ax.axvline(gap_start, color='red', linestyle='--', linewidth=2, label='Batas Material')

# Tambahkan legend setelah garis dibuat
ax.legend(loc='upper right')

# Colorbar dan animasi
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Amplitudo (dinormalisasi)")
ani = animation.FuncAnimation(fig, update, frames=steps, interval=30, blit=False)
plt.title("Gelombang 2D dengan Dua Jenis Material")
# Simpan animasi sebagai video MP4
# ani.save("simulasi_gelombang_beda_material.mp4", writer='ffmpeg', fps=30)
plt.show()