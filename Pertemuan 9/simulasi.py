import numpy as np
import matplotlib.pyplot as plt 

def simulasi_hantaran_panas(L, tol=1e-4, max_iter=10000):
    plt.figure(figsize=(6, 5))
    for it in range(1, 6):
        # Inisialisasi grid suhu
        T = np.zeros((L, L))
        
        # Set kondisi batas
        T[:, 0] = 100     # Batas kiri
        T[:, -1] = 0      # Batas kanan
        T[0, :] = 50      # Batas bawah
        T[-1, :] = 25     # Batas atas
        
        # Iterasi Gauss-Seidel
        for iteration in range(max_iter):
            T_old = T.copy()
            
            # Update nilai internal points
            for i in range(1, L-1):
                for j in range(1, L-1):
                    T[i,j] = 0.25 * (T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1])
            
            # Check konvergensi
            if np.max(np.abs(T - T_old)) < tol:
                break
        
        # Plot hasil setiap iterasi tertentu
        if it in [1, 2, 3, 4, 5]:
            plt.subplot(3, 2, it)
            plt.imshow(T, cmap='hot', origin='lower')
            plt.colorbar(label='Suhu (°C)')
            plt.title(f'Iterasi ke-{iteration+1}')
            plt.xlabel('X')
            plt.ylabel('Y')
plt.tight_layout()
plt.show()
    