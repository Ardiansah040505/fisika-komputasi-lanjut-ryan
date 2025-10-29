import serial
import time
import csv
import threading
import matplotlib.pyplot as plt
from collections import deque
from datetime import datetime

# === Konfigurasi ===
SERIAL_PORT = '/dev/ttyUSB0'   # Ganti sesuai port Anda (cek dengan `ls /dev/tty*`)
BAUD_RATE = 9600               # Sesuaikan dengan modul Geiger Anda
CONVERSION_FACTOR = 0.0057     # Contoh: CPM → µSv/h (untuk tabung SBM-20)
LOG_FILE = 'geiger_log.csv'

# === Variabel global ===
counts = 0
cpm_history = deque(maxlen=60)  # Simpan 60 menit terakhir
time_history = deque(maxlen=60)
stop_event = threading.Event()

# === Fungsi pembaca serial ===
def read_geiger():
    global counts
    try:
        with serial.Serial(SERIAL_PORT, BAUD_RATE, timeout=1) as ser:
            print(f"Terhubung ke {SERIAL_PORT} pada {BAUD_RATE} baud...")
            while not stop_event.is_set():
                if ser.in_waiting > 0:
                    # Anggap setiap byte = 1 deteksi
                    ser.read()  # Baca dan abaikan isi, hanya hitung event
                    counts += 1
                time.sleep(0.01)  # Ringan di CPU
    except serial.SerialException as e:
        print(f"Error serial: {e}")
        stop_event.set()

# === Fungsi utama penghitung & logging ===
def main_loop():
    global counts
    # Buat header CSV
    with open(LOG_FILE, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Timestamp', 'CPM', 'uSv_h'])

    last_minute = time.time()
    print("Mulai monitoring radiasi...\nTekan Ctrl+C untuk berhenti.\n")

    plt.ion()
    fig, ax = plt.subplots(figsize=(10, 4))
    line, = ax.plot([], [], 'b-o')
    ax.set_xlabel('Waktu (menit terakhir)')
    ax.set_ylabel('CPM')
    ax.set_title('Geiger Counter - Real-time CPM')
    ax.grid(True)

    try:
        while not stop_event.is_set():
            now = time.time()
            if now - last_minute >= 60:
                cpm = counts  # karena 60 detik
                uSv_h = cpm * CONVERSION_FACTOR
                timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

                # Simpan ke CSV
                with open(LOG_FILE, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([timestamp, cpm, uSv_h])

                # Update riwayat untuk grafik
                cpm_history.append(cpm)
                time_history.append(len(cpm_history))

                # Tampilkan di terminal
                print(f"[{timestamp}] CPM: {cpm} | Estimasi: {uSv_h:.3f} µSv/h")

                # Update grafik real-time
                line.set_xdata(list(time_history))
                line.set_ydata(list(cpm_history))
                ax.relim()
                ax.autoscale_view()
                plt.draw()
                plt.pause(0.01)

                # Reset hitungan
                counts = 0
                last_minute = now

            time.sleep(0.1)

    except KeyboardInterrupt:
        print("\nMenghentikan monitoring...")
    finally:
        stop_event.set()
        plt.ioff()
        plt.show()  # Biarkan grafik tetap terbuka

# === Jalankan ===
if __name__ == '__main__':
    reader_thread = threading.Thread(target=read_geiger, daemon=True)
    reader_thread.start()
    main_loop()
    reader_thread.join(timeout=1)

