
# Animasi gerak jatuh bebas dengan pantulan tanah
# Dibuat oleh: Yoyok Adisetio Laksono
# Departemen Fisika, FMIPA, Universitas Negeri Malang
# 2025
# Program ini adalah materi di kuliah Fisika Komputasi Lanjut

import matplotlib.pyplot as plt
import matplotlib.animation as animation
#from tqdm import tqdm

#######################################################################################################
# Definisi fungsi

# Generator posisi bola dari waktu ke waktu
def get_pos(t=0):           # Fungsi generator dengan t adalah waktu awal = 0
    x, y, vx, vy = x0, y0, vx0, vy0
    while x < XMAX:
        t += dt             # waktu t bertambah tiap langkah dt
        x += vx0 * dt       # menambah jarak gerak horizontal konstan
        y += vy * dt        # menambah jarak gerak vertikal
        vy -= g * dt        # kecepatan vertikal akibat percepatan gravitasi
        if y < R:           # jika menyentuh tanah
            y = R           # set y setinggi R agar bola tidak nampak tenggelam
            vy = -vy * cor  # pantulan dengan mengubah polaritas dan dikurangi akibat energi hilang
        yield x, y          # hasilkan posisi sekarang sebagai tuple (x,y) untuk animasi selanjutnya

# Inisialisasi tampilan awal animasi
def init():
    ax.set_xlim(0, XMAX)          # mengatur batas aksis x
    ax.set_ylim(0, y0)            # mengatur batas aksis y
    ax.set_xlabel("$x$ (m)")      # teks keterangan aksis x, $x$ untuk kode latex, x tampil italic 
    ax.set_ylabel("$y$ (m)")      # teks keterangan aksis x
    line.set_data(xdata, ydata)   # buat garis dari xdata, ydata
    ball.set_center((x0, y0))     # pengaturan lokasi tengan bola sebagai koordinat acuan bola
    height_text.set_text(f"Ketinggian: {y0:.2f} m")  # menulis variabel y0 ke obyek height_text  
    return line, ball, height_text  # mengirim hasil pengubahan ke pemanggil fungsi init

# Fungsi untuk memperbarui animasi setiap frame
def animate(pos):
    x, y = pos                    # ambil x dan y dari parameter pos
    xdata.append(x)               # tambahkan x ke xdata
    ydata.append(y)               # tambahkan y ke ydata
    line.set_data(xdata, ydata)   # masukkan ke line
    ball.set_center((x, y))       # letakkan bola di posisi x, y
    height_text.set_text(f"Ketinggian: {y:.2f} m")
    return line, ball, height_text

#######################################################################################################
# Parameter fisika dan simulasi
g = 9.81        # Percepatan gravitasi (m/s^2)
XMAX = 5       # Batas horizontal lintasan bola
cor = 0.65      # Koefisien restitusi (pantulan)
dt = 0.005      # Langkah waktu untuk animasi

# Posisi dan kecepatan awal
x0, y0    = 0, 4
vx0, vy0  = 1, 0
R         = 0.08  # Jari-jari bola

#######################################################################################################
# Atur figure dan objek visual
fig, ax = plt.subplots()         # fig adalah gambar keseluruhan, ax adalah gambar grafik di dalamnya
ax.set_aspect("equal")           # Mengatur ukuran sumbu x dan y sama
(line,) = ax.plot([], [], lw=2)  # garis lintasan dengan data x dan y kosong
ball = plt.Circle((x0, y0), R)   # bola dengan radius R ==> belum ada di axes, masih definisi
ax.add_patch(ball)               # tambahkan bola ke axes
xdata, ydata = [], []            # data lintasan masih kosong
# teks ketinggian, f di depan string adalah f-string dimana string bisa menampilkan nilai variabel 
# dengan menulis {variabel}
height_text = ax.text(XMAX * 0.5, y0 * 0.8, f"Ketinggian: {y0:.2f} m") 

#######################################################################################################
# Jalankan animasi
interval = 2000 * dt
ani = animation.FuncAnimation(
                fig=fig,            # obyek gambar yang akan dianimasikan
                func=animate,       # fungsi animasi
                frames=get_pos,     # penghitungan 
                blit=True,
                interval=interval,
                repeat=False,
                init_func=init,
                cache_frame_data=False
                )
plt.show()