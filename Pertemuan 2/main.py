import numpy as np
import matplotlib.pyplot as plt

# Buat data x dari 0 sampai 2π dengan 100 titik
x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)
z = np.cos(x)
yz = y * z

# Buat plot
plt.figure(figsize=(8, 4))
plt.plot(x, y, label='sin(x)', color='red', linewidth=2, linestyle ="dashed")
plt.plot(x, z, label='cos(x)', color='blue', linewidth=2, linestyle ="dashed")
plt.plot(x, yz, label='cos(x)', color='green', linewidth=2, linestyle ="dashed")
# Tambahkan label dan judul
plt.title('Grafik Fungsi Sinus')
plt.xlabel('x (radian)')
plt.ylabel('sin(x) dan cos(x)')
plt.grid(True)
plt.legend()

# Tampilkan plot
plt.show()