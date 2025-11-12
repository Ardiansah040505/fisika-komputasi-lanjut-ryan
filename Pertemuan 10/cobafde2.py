import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# ============================
# PARAMETER SIMULASI
# ============================
L = 1.0           # Panjang domain [m]
N = 50            # Jumlah grid points
a = 0.2           # Konduktivitas termal
d = L / (N - 1)   # Spasi grid
h = 0.8 * (d**2) / (4 * a)  # Time step (memenuhi CFL)

# Grid coordinates
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
X, Y = np.meshgrid(x, y)

# ============================
# KONDISI AWAL & SYARAT BATAS
# ============================
u = np.zeros((N, N))  # Suhu awal: 0°C di seluruh domain

# Sumber panas di titik tengah
center_i, center_j = N // 2, N // 2
u[center_i, center_j] = 1.0  # Suhu 1°C di tengah

# Syarat batas Dirichlet (suhu tetap 0 di semua batas)
u[0, :] = 0    # Batas bawah
u[-1, :] = 0   # Batas atas
u[:, 0] = 0    # Batas kiri
u[:, -1] = 0   # Batas kanan

# ============================
# SIMULASI
# ============================
def update_temperature(u):
    """Update suhu menggunakan skema eksplisit"""
    u_new = u.copy()
    
    # Loop melalui titik interior (tidak termasuk batas)
    for i in range(1, N-1):
        for j in range(1, N-1):
            u_new[i, j] = u[i, j] + (a * h / d**2) * (
                u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1] - 4*u[i, j]
            )
    
    # Pertahankan syarat batas
    u_new[0, :] = 0
    u_new[-1, :] = 0
    u_new[:, 0] = 0
    u_new[:, -1] = 0
    
    return u_new

# ============================
# VISUALISASI
# ============================
fig = plt.figure(figsize=(12, 5))

# Plot 2D
ax1 = fig.add_subplot(121)
im = ax1.imshow(u, cmap='hot', origin='lower', extent=[0, L, 0, L], vmin=0, vmax=1)
plt.colorbar(im, label='Suhu (°C)')
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_title('Distribusi Suhu 2D')

# Plot 3D
ax2 = fig.add_subplot(122, projection='3d')
surf = ax2.plot_surface(X, Y, u, cmap='hot', linewidth=0, antialiased=True)
ax2.set_xlabel('x (m)')
ax2.set_ylabel('y (m)')
ax2.set_zlabel('Suhu (°C)')
ax2.set_title('Distribusi Suhu 3D')
ax2.set_zlim(0, 1)

plt.tight_layout()

# ============================
# ANIMASI
# ============================
def animate(frame):
    global u
    u = update_temperature(u)
    
    # Update plot 2D
    ax1.clear()
    im = ax1.imshow(u, cmap='hot', origin='lower', extent=[0, L, 0, L], vmin=0, vmax=1)
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_title(f'Distribusi Suhu 2D - Langkah waktu: {frame}')
    
    # Update plot 3D
    ax2.clear()
    surf = ax2.plot_surface(X, Y, u, cmap='hot', linewidth=0, antialiased=True)
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('y (m)')
    ax2.set_zlabel('Suhu (°C)')
    ax2.set_title(f'Distribusi Suhu 3D - Langkah waktu: {frame}')
    ax2.set_zlim(0, 1)
    
    return [im, surf]

# Buat animasi
ani = FuncAnimation(fig, animate, frames=100, interval=50, blit=False)

plt.show()

# ============================
# SIMPAN VIDEO (opsional)
# ============================
# Uncomment baris berikut untuk menyimpan video
ani.save('simulasi_panas_2d.mp4', writer='ffmpeg', fps=10)
print("Simulasi selesai!")
print(f"Parameter: L = {L}, N = {N}, a = {a}")
print(f"Spasi grid: d = {d:.4f}")
print(f"Time step: h = {h:.6f}")
print(f"CFL number: {a * h / d**2:.4f} (harus ≤ 0.25)")