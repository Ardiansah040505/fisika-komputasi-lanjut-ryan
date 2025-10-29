import numpy as np
import matplotlib.pyplot as plt

def solve_steady_state_heat(nx, ny, T_left, T_right, T_top, T_bottom, max_iter=1000, tolerance=1e-4):
    """
    Menyelesaikan persamaan Laplace 2D untuk konduksi panas steady-state
    menggunakan metode finite difference dengan iterasi Gauss-Seidel
    
    Parameters:
    -----------
    nx, ny : int
        Jumlah titik grid dalam arah x dan y
    T_left, T_right, T_top, T_bottom : float
        Suhu pada batas-batas (boundary conditions)
    max_iter : int
        Maksimum iterasi
    tolerance : float
        Toleransi konvergensi
    
    Returns:
    --------
    T : ndarray
        Matrix suhu 2D
    iterations : int
        Jumlah iterasi yang dilakukan
    """
    
    # Inisialisasi grid suhu
    T = np.zeros((ny, nx))
    
    # Set kondisi batas
    T[:, 0] = T_left     # Batas kiri
    T[:, -1] = T_right   # Batas kanan
    T[0, :] = T_bottom   # Batas bawah
    T[-1, :] = T_top     # Batas atas
    
    # Iterasi Gauss-Seidel
    for it in range(max_iter):
        T_old = T.copy()
        
        # Update nilai internal points
        for i in range(1, ny-1):
            for j in range(1, nx-1):
                T[i,j] = 0.25 * (T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1])
        
        # Check konvergensi
        if np.max(np.abs(T - T_old)) < tolerance:
            return T, it + 1
            
    return T, max_iter

# Parameter simulasi
nx = ny = 50  # Ukuran grid
T_left = 100   # Suhu pada batas kiri (°C)
T_right = 0    # Suhu pada batas kanan (°C)
T_top = 50     # Suhu pada batas atas (°C)
T_bottom = 25  # Suhu pada batas bawah (°C)

# Jalankan simulasi
T, iterations = solve_steady_state_heat(nx, ny, T_left, T_right, T_top, T_bottom)

# Visualisasi hasil
plt.figure(figsize=(10, 8))

# Plot distribusi suhu
plt.subplot(111)
im = plt.imshow(T, cmap='hot', interpolation='nearest')
plt.colorbar(im, label='Temperatur (°C)')
plt.title(f'Distribusi Suhu Steady-State\nKonvergen setelah {iterations} iterasi')
plt.xlabel('Posisi X')
plt.ylabel('Posisi Y')

# Tampilkan plot
plt.tight_layout()
plt.show()

# Print informasi tambahan
print(f"Simulasi selesai setelah {iterations} iterasi")
print(f"Suhu maksimum: {np.max(T):.2f}°C")
print(f"Suhu minimum: {np.min(T):.2f}°C")
print(f"Suhu rata-rata: {np.mean(T):.2f}°C")