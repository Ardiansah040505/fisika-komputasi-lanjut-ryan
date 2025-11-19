import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameter domain
nx, ny = 200, 200         # jumlah grid
dx = dy = 1.0             # resolusi spasial
c = 1.0                   # kecepatan gelombang
dt = 0.5                  # langkah waktu
steps = 500               # jumlah langkah waktu

# Stabilitas (CFL condition)
assert c * dt / dx < 1 / np.sqrt(2), "CFL condition tidak terpenuhi!"

# Inisialisasi grid
u = np.zeros((nx, ny))        # kondisi saat ini
u_prev = np.zeros((nx, ny))   # kondisi sebelumnya
u_next = np.zeros((nx, ny))   # kondisi berikutnya

# Sumber Gaussian negatif di tengah
x0, y0 = nx // 2, ny // 2
sigma = 5.0
gaussian = -np.exp(-(((np.arange(nx)[:, None] - x0)**2 + (np.arange(ny)[None, :] - y0)**2) / (2 * sigma**2)))
u = gaussian.copy()
u_prev = gaussian.copy()

# Fungsi update
def update(frame):
    global u, u_prev, u_next
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            laplacian = (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1] - 4 * u[i, j]) / dx**2
            u_next[i, j] = 2 * u[i, j] - u_prev[i, j] + (c * dt)**2 * laplacian

    im.set_array(u_next)
    u_prev, u = u, u_next.copy()
    return [im]

# Visualisasi animasi
fig, ax = plt.subplots()
im = ax.imshow(u, cmap='seismic', vmin=-1, vmax=1, animated=True)
ani = animation.FuncAnimation(fig, update, frames=steps, interval=30, blit=True)
plt.title("Simulasi Gelombang 2D")
ani.save("simulasi_gelombang_2D.mp4", writer='ffmpeg', fps=30)
plt.show()