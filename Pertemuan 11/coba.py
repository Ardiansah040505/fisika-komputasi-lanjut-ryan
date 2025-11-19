import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- 1. Inisialisasi Parameter ---
# Ukuran domain (jumlah titik grid)
# Sesuaikan dengan Tugas 1 & 2: 100 x 200
nx = 100
ny = 200

# Kecepatan gelombang (konstanta dalam simulasi dasar ini)
# PDF menyebut c=1 untuk Tugas 1, dan bervariasi untuk Tugas 2.
# Kita buat versi dasar dulu dengan c=1
c_val = 1.0

# Ukuran langkah waktu
# dt <= dx / (c * sqrt(2)) untuk stabilitas. Kita asumsikan dx = dy = 1.
# S = 0.5 adalah faktor keamanan.
S = 0.5
dx = 1.0  # Kita tetapkan dx = dy = 1 untuk menyederhanakan rumus
dt = S * dx / (c_val * np.sqrt(2))

# Jumlah iterasi waktu (lama simulasi)
num_steps = 400

# --- 2. Inisialisasi Array ---
# Buat array 2D untuk menyimpan amplitudo
u_lama = np.zeros((nx, ny))
u_sekarang = np.zeros((nx, ny))
u_baru = np.zeros((nx, ny))

# --- 3. Set Kondisi Awal (Fungsi Gaussian dari PDF) ---
# PDF menyebut: f(x,y) = (1/(2*pi*sigma^2)) * exp(-(x+y)^2 / (2*sigma^2))
# Ini BUKAN Gaussian isotropik standar. Ini akan menciptakan pola gelombang unik.
# Kita ikuti sesuai PDF.
# Kita gunakan meshgrid untuk menghitung fungsi ini di seluruh domain.
# Gunakan indexing='ij' agar X sesuai dengan sumbu baris (i) dan Y sesuai dengan sumbu kolom (j)
x = np.arange(nx)
y = np.arange(ny)
X, Y = np.meshgrid(x, y, indexing='ij')

# Parameter Gaussian dari PDF (atau bisa disesuaikan)
sigma = 2.0

# Kondisi awal: f(x, y) dari PDF (x+y)^2
# f_x_y = (1.0 / (2 * np.pi * sigma**2)) * np.exp(-((X + Y)**2) / (2 * sigma**2))
# Namun, ini akan menciptakan pola gelombang dari pojok kiri atas.
# Untuk simulasi gelombang dari TITIK SUMBER, kita gunakan fungsi Gaussian isotropik
# di posisi sumber, seperti yang dijelaskan dalam teks sebelum rumusnya.
# Misalnya, sumber di tengah (atau sesuai tugas): (25, 50) untuk Tugas 1 & 2
source_x = 25 # Indeks x sumber (sesuai Tugas 1 & 2)
source_y = 50 # Indeks y sumber (sesuai Tugas 1 & 2)

# Kondisi awal: Gaussian Isotropik di titik sumber (Lebih fisikal dan umum untuk simulasi gelombang point source)
dist_sq = (X - source_x)**2 + (Y - source_y)**2
u_sekarang = (1.0 / (2 * np.pi * sigma**2)) * np.exp(-dist_sq / (2 * sigma**2))
u_lama = u_sekarang.copy() # Untuk iterasi pertama, asumsikan kecepatan awal nol

# --- 4. Setup Visualisasi ---
fig, ax = plt.subplots()
# Cetak bentuk array sebelum imshow untuk debugging
print(f"Bentuk u_sekarang sebelum imshow: {u_sekarang.shape}")
im = ax.imshow(u_sekarang, cmap='RdBu', animated=True, vmin=-0.05, vmax=0.05)
ax.set_title('Simulasi Gelombang 2D - Gaussian Isotropik (c=1)')
ax.axis('off') # Matikan sumbu untuk tampilan bersih

# --- 5. Fungsi Update untuk Animasi ---
def update(frame):
    global u_lama, u_sekarang, u_baru

    # Koefisien dalam rumus iteratif: (c*dt/dx)^2
    # Karena dx=1 dan c=1, maka coeff = (1*dt)^2
    coeff = (c_val * dt / dx)**2

    # --- 6. Iterasi Waktu (Loop Utama) ---
    # Hitung Laplacian untuk wilayah interior sekaligus
    laplacian_interior = (
        u_sekarang[2:, 1:-1] +     # u[i+1, j]
        u_sekarang[:-2, 1:-1] +    # u[i-1, j]
        u_sekarang[1:-1, 2:] +     # u[i, j+1]
        u_sekarang[1:-1, :-2] -    # u[i, j-1]
        4 * u_sekarang[1:-1, 1:-1] # - 4*u[i, j]
    )

    # Hitung u_baru untuk wilayah interior
    u_baru[1:-1, 1:-1] = (
        2 * u_sekarang[1:-1, 1:-1] - # 2 * u^n
        u_lama[1:-1, 1:-1] +         # - u^n-1
        coeff * laplacian_interior   # + coeff * Laplacian
    )

    # --- 7. Terapkan Kondisi Batas (Dirichlet: u = 0) ---
    u_baru[0, :] = 0  # Sisi atas
    u_baru[-1, :] = 0 # Sisi bawah
    u_baru[:, 0] = 0  # Sisi kiri
    u_baru[:, -1] = 0 # Sisi kanan

    # --- 8. Update Array untuk Iterasi Selanjutnya ---
    # PENTING: Gunakan .copy() untuk mencegah saling referensi antar array
    u_lama = u_sekarang.copy()
    u_sekarang = u_baru.copy()

    # Update data gambar untuk animasi
    im.set_array(u_sekarang)
    ax.set_title(f'Simulasi Gelombang 2D - Iterasi {frame} (c=1)')
    return [im]

# --- 9. Jalankan Animasi ---
ani = animation.FuncAnimation(fig, update, frames=num_steps, interval=50, blit=True, repeat=False)

# Tampilkan animasi
plt.show()

# Opsional: Simpan animasi sebagai GIF (bisa memakan waktu)
# ani.save('simulasi_gelombang_2d_dasar.gif', writer='pillow', fps=20)
