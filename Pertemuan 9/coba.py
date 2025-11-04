import numpy as np
import matplotlib.pyplot as plt

#parameter
L = 10.0  # panjang domain dalam meter
N = 100  # jumlah grid dalam setiap arah
W = 5.0   # lebar domain dalam meter
dx = dy = L / (N-1)  # ukuran grid

#inisialisasi temperatur
T = np.zeros((N, N))


#syarat batas
T[:, 0] = 0.0  # kiri
T[:, -1] = 50.0  # kanan
T[0, :] = 100.0   # atas
T[-1, :] = 0.0   # bawah

#iterasi Gauss-Seidel
max_iter = 1000
tolerance = 1e-5
for iteration in range(max_iter):
    T_old = T.copy()

    for i in range(1, N-1):
        for j in range(1, N-1):
            T[i, j] = 0.25 * (T[i+1, j] + T[i-1, j] + T[i, j+1] + T[i, j-1])
    
    #cek konvergensi
    np.max(np.abs(T - T_old))
    error = np.max(np.abs(T - T_old ))
    if error < tolerance:
        print(f"Konvergen setelah {iteration + 1} iterasi.")
        break

#visualisasi hasil
plt.figure(figsize=(8, 4))
plt.imshow(T, cmap='hot', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Temperatur (°C)')
plt.xlabel('Panjang (m)')
plt.ylabel('Lebar (m)')
plt.title('Distribusi Temperatur Steady-State')
plt.grid()
plt.show()
