import tkinter as tk
from tkinter import ttk
import threading
import time
import csv
import serial
import serial.tools.list_ports
from datetime import datetime

# === Konfigurasi ===
BAUD_RATE = 9600
CONVERSION_FACTOR = 0.0057  # Sesuaikan: SBM-20 ≈ 0.0057, J305 ≈ 0.00812
CHECK_INTERVAL = 2  # detik

class GeigerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Geiger Counter - Windows")
        self.root.geometry("520x320")
        self.root.resizable(False, False)

        self.status_var = tk.StringVar(value="Menunggu perangkat Geiger...")
        self.cpm_var = tk.StringVar(value="CPM: --")
        self.dose_var = tk.StringVar(value="Dosis: -- µSv/h")

        tk.Label(root, textvariable=self.status_var, font=("Arial", 12)).pack(pady=10)
        tk.Label(root, textvariable=self.cpm_var, font=("Arial", 16)).pack()
        tk.Label(root, textvariable=self.dose_var, font=("Arial", 16)).pack()

        self.log_text = tk.Text(root, height=8, state='disabled', font=("Consolas", 10))
        self.log_text.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

        # File log
        self.log_file = open("geiger_log.csv", "w", newline="", encoding="utf-8")
        self.csv_writer = csv.writer(self.log_file)
        self.csv_writer.writerow(["Timestamp", "CPM", "uSv_h"])
        self.log_file.flush()

        # State
        self.running = False
        self.counts = 0
        self.start_time = time.time()
        self.detected_port = None

        # Mulai pemindaian
        self.log("Memindai port COM setiap 2 detik...")
        self.scan_thread = threading.Thread(target=self.scan_loop, daemon=True)
        self.scan_thread.start()

    def log(self, msg):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.config(state='normal')
        self.log_text.insert(tk.END, f"[{timestamp}] {msg}\n")
        self.log_text.see(tk.END)
        self.log_text.config(state='disabled')

    def find_geiger_port(self):
        """Cari port COM yang aktif (asumsi: Geiger kirim data serial)."""
        ports = serial.tools.list_ports.comports()
        for port in ports:
            # Filter: abaikan perangkat Bluetooth, printer, dll
            if "Arduino" in port.description or "CH340" in port.description or "CP210" in port.description or "USB Serial" in port.description:
                return port.device
            # Atau: coba semua COM port (opsi lebih agresif)
        # Jika tidak yakin, kembalikan semua COM port aktif
        active_ports = [p.device for p in ports if 'COM' in p.device]
        return active_ports[0] if active_ports else None

    def scan_loop(self):
        while not self.running:
            port = self.find_geiger_port()
            if port and port != self.detected_port:
                self.detected_port = port
                self.status_var.set(f"Terhubung ke {port}")
                self.log(f"Perangkat terdeteksi di {port}!")
                self.running = True
                self.read_thread = threading.Thread(target=self.read_serial, daemon=True)
                self.read_thread.start()
                break
            time.sleep(CHECK_INTERVAL)

    def read_serial(self):
        try:
            with serial.Serial(self.detected_port, BAUD_RATE, timeout=1) as ser:
                self.log("Mulai membaca data radiasi...")
                while self.running:
                    if ser.in_waiting > 0:
                        ser.read()  # Anggap setiap byte = 1 click
                        self.counts += 1
                    time.sleep(0.01)
        except Exception as e:
            self.log(f"Error: {e}")
            self.reset()

    def update_display(self):
        if not self.running:
            return

        elapsed = time.time() - self.start_time
        if elapsed >= 60:
            cpm = self.counts
            uSv_h = cpm * CONVERSION_FACTOR

            # Update GUI
            self.cpm_var.set(f"CPM: {cpm}")
            self.dose_var.set(f"Dosis: {uSv_h:.3f} µSv/h")

            # Simpan ke CSV
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            self.csv_writer.writerow([timestamp, cpm, uSv_h])
            self.log_file.flush()

            self.log(f"CPM: {cpm} | {uSv_h:.3f} µSv/h")

            # Reset hitungan
            self.counts = 0
            self.start_time = time.time()

        self.root.after(1000, self.update_display)

    def reset(self):
        self.running = False
        self.detected_port = None
        self.status_var.set("Menunggu perangkat Geiger...")
        self.cpm_var.set("CPM: --")
        self.dose_var.set("Dosis: -- µSv/h")
        self.scan_thread = threading.Thread(target=self.scan_loop, daemon=True)
        self.scan_thread.start()

    def on_close(self):
        self.running = False
        self.log_file.close()
        self.root.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    app = GeigerApp(root)
    root.protocol("WM_DELETE_WINDOW", app.on_close)
    root.after(1000, app.update_display)
    root.mainloop()