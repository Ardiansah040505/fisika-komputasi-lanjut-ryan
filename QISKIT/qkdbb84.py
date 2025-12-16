"""
bb84_pyqt5_full.py
GUI BB84 simulator (PyQt5) — robust against Qiskit/Aer availability.
- If qiskit + qiskit-aer present: uses AerSimulator or execute (faster/more realistic).
- If qiskit Aer missing: falls back to a simple classical BB84 emulator so the GUI still runs.

Reference: original script uploaded by user (main.py). :contentReference[oaicite:1]{index=1}
"""

import sys
import traceback
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from functools import partial

# PyQt5 imports
from PyQt5.QtCore import QObject, QThread, pyqtSignal, Qt
from PyQt5.QtWidgets import (
    QApplication, QComboBox, QDoubleSpinBox, QGridLayout, QLabel,
    QMainWindow, QPushButton, QSpinBox, QTextEdit, QWidget, QMessageBox
)

AER_BACKEND = AerSimulator()
# Try to import Qiskit in a resilient way
USE_QISKIT = False
USE_EXECUTE = False
USE_AER_SIM = False
_qiskit_import_exc = None

try:
    # try modern-style imports first (Qiskit >= 1.0)
    from qiskit import QuantumCircuit, transpile
    from qiskit.providers.aer import AerSimulator
    AER_BACKEND = AerSimulator()
    USE_QISKIT = True
    USE_AER_SIM = True
except Exception as e_modern:
    _qiskit_import_exc = e_modern
    try:
        # try legacy-style (older qiskit that exposes Aer & execute)
        from qiskit import QuantumCircuit, Aer, execute, transpile
        AER_BACKEND = Aer.get_backend("qasm_simulator")
        USE_QISKIT = True
        USE_EXECUTE = True
    except Exception as e_legacy:
        # Qiskit not usable; fallback will be used
        _qiskit_import_exc = e_legacy
        USE_QISKIT = False
        AER_BACKEND = None

# Default number of qubits if user doesn't change in GUI
DEFAULT_NUM_QUBITS = 32


# -------------------------
# Core BB84 logic functions
# -------------------------
def generate_bits(n: int) -> np.ndarray:
    return np.random.randint(2, size=n)


def encode_message_qiskit(bits: np.ndarray, basis: np.ndarray, num_qubits: int) -> list:
    """
    Build list of single-qubit QuantumCircuit objects representing encoded qubits.
    Used only when Qiskit is available.
    """
    message = []
    for i in range(num_qubits):
        qc = QuantumCircuit(1, 1)
        if basis[i] == 0:  # Z-basis
            if bits[i] == 1:
                qc.x(0)
        else:  # X-basis
            if bits[i] == 0:
                qc.h(0)
            else:
                qc.x(0)
                qc.h(0)
        qc.barrier()
        message.append(qc)
    return message


def simulate_quantum_channel_qiskit(message: list, error_rate: float) -> list:
    """
    Apply bit-flip (X) to each circuit with probability error_rate.
    Works when representing qubits as QuantumCircuit objects.
    """
    noisy = []
    for qc in message:
        qc2 = qc.copy()
        if np.random.random() < error_rate:
            qc2.x(0)
        noisy.append(qc2)
    return noisy


def measure_message_qiskit(message: list, basis: np.ndarray, num_qubits: int) -> list:
    """
    Measure list of single-qubit circuits using Qiskit (either execute or AerSimulator.run).
    Returns list of ints (0/1).
    """
    if USE_EXECUTE:
        measurements = []
        for q in range(num_qubits):
            qc = message[q].copy()
            if basis[q] == 1:
                qc.h(0)
            qc.measure(0, 0)
            result = execute(qc, AER_BACKEND, shots=1, memory=True).result()
            measurements.append(int(result.get_memory()[0]))
        return measurements

    if USE_AER_SIM:
        # Build circuits and run in batch
        circuits = []
        for q in range(num_qubits):
            qc = message[q].copy()
            if basis[q] == 1:
                qc.h(0)
            qc.measure(0, 0)
            circuits.append(qc)
        tc = transpile(circuits, backend=AER_BACKEND)
        job = AER_BACKEND.run(tc, shots=1, memory=True)
        result = job.result()
        measurements = [int(result.get_memory(i)[0]) for i in range(len(circuits))]
        return measurements

    raise RuntimeError("Qiskit backend not available for qiskit-based measurement.")


# -------------------------
# Fallback pure-Python BB84 emulator
# -------------------------
def encode_message_fallback(bits: np.ndarray, basis: np.ndarray, num_qubits: int):
    """
    Fallback encoding: we represent each qubit as (alice_basis, bit).
    This is not a full quantum simulation but reproduces BB84 behavior:
    - If Bob measures in same basis -> he recovers Alice's bit (unless channel error flips it).
    - If Bob measures in different basis -> result is random (0/1).
    """
    # return list of tuples for each qubit
    return [(int(basis[i]), int(bits[i])) for i in range(num_qubits)]


def simulate_quantum_channel_fallback(message: list, error_rate: float) -> list:
    """
    Simulate an X error: for simplicity, we flip the stored logical bit with probability error_rate.
    This approximates a noisy channel effect on the resulting measured bits.
    """
    noisy = []
    for (ab, bit) in message:
        bit2 = bit
        if np.random.random() < error_rate:
            bit2 = 1 - bit2
        noisy.append((ab, bit2))
    return noisy


def measure_message_fallback(message: list, basis: np.ndarray, num_qubits: int) -> list:
    """
    Measure the fallback message representation:
    - If bob_basis == alice_basis: return stored bit
    - Else: return random 0/1
    This matches ideal BB84 probabilities.
    """
    results = []
    for i in range(num_qubits):
        ab, bit = message[i]
        if basis[i] == ab:
            results.append(int(bit))
        else:
            results.append(int(np.random.randint(2)))
    return results


# -------------------------
# Common utility functions
# -------------------------
def remove_garbage(a_basis: np.ndarray, b_basis: np.ndarray, bits):
    """
    Keep only bits where bases match. bits may be list or numpy array.
    """
    # bits can be numpy array or list; assume indexable
    kept = [bits[i] for i in range(len(a_basis)) if a_basis[i] == b_basis[i]]
    return kept


def check_keys(key1: list, key2: list):
    equal = (key1 == key2)
    return equal, key1, key2


# -------------------------
# Worker for running simulations in background thread
# -------------------------
class SimulatorWorker(QObject):
    progress = pyqtSignal(str)
    finished = pyqtSignal()
    error = pyqtSignal(str)

    def __init__(self, num_qubits: int, error_rate: float, iterations: int, mode: str):
        super().__init__()
        self.num_qubits = int(num_qubits)
        self.error_rate = float(error_rate)
        self.iterations = int(iterations)
        self.mode = mode
        self._is_running = True

    def stop(self):
        self._is_running = False

    def run(self):
        try:
            self.progress.emit(f"Using Qiskit: {USE_QISKIT}, execute: {USE_EXECUTE}, AerSim: {USE_AER_SIM}")
            if USE_QISKIT:
                self.progress.emit("Running using Qiskit backend.")
            else:
                self.progress.emit("Qiskit/Aer not available — using fallback classical emulator.")

            for it in range(self.iterations):
                if not self._is_running:
                    self.progress.emit("Simulation stopped by user.")
                    break

                # Alice generates bits & basis
                alice_bits = generate_bits(self.num_qubits)
                alice_basis = generate_bits(self.num_qubits)

                if USE_QISKIT:
                    # build quantum circuits
                    message = encode_message_qiskit(alice_bits, alice_basis, self.num_qubits)
                    message = simulate_quantum_channel_qiskit(message, self.error_rate)
                    if self.mode == "eavesdrop":
                        # simple eavesdrop: measure in random bases and re-prepare (approx)
                        # For simplicity, emulate eavesdrop by random measurement disturbance:
                        for i in range(len(message)):
                            if np.random.random() < 0.5:
                                # measure in Z by applying measure - we just flip sometimes
                                if np.random.random() < 0.1:
                                    message[i].x(0)
                            else:
                                # measure in X: apply H then measure; approximate disturbance
                                if np.random.random() < 0.1:
                                    message[i].x(0)
                    bob_basis = generate_bits(self.num_qubits)
                    bob_results = measure_message_qiskit(message, bob_basis, self.num_qubits)

                else:
                    # fallback path (no Qiskit)
                    message = encode_message_fallback(alice_bits, alice_basis, self.num_qubits)
                    message = simulate_quantum_channel_fallback(message, self.error_rate)
                    if self.mode == "eavesdrop":
                        # eavesdropper introduces additional disturbance: flip some bits randomly
                        noisy = []
                        for (ab, bit) in message:
                            if np.random.random() < 0.1:
                                bit = 1 - bit
                            noisy.append((ab, bit))
                        message = noisy
                    bob_basis = generate_bits(self.num_qubits)
                    bob_results = measure_message_fallback(message, bob_basis, self.num_qubits)

                alice_key = remove_garbage(alice_basis, bob_basis, alice_bits)
                bob_key = remove_garbage(alice_basis, bob_basis, bob_results)
                eq, akey, bkey = check_keys(alice_key, bob_key)
                header = f"Iteration {it+1}/{self.iterations} | Qubits: {self.num_qubits} | Error: {self.error_rate} | Mode: {self.mode}"
                self.progress.emit(header)
                self.progress.emit(f"Alice key ({len(akey)}): {''.join(map(str, akey))}")
                self.progress.emit(f"Bob   key ({len(bkey)}): {''.join(map(str, bkey))}")
                self.progress.emit("Status: " + ("MATCH" if eq else "DIFFER (noise/eavesdrop detected)"))
                self.progress.emit("-" * 60)

            self.finished.emit()
        except Exception as ex:
            tb = traceback.format_exc()
            self.error.emit(tb)
            self.finished.emit()


# -------------------------
# PyQt5 Main Window
# -------------------------
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("BB84 Simulator - PyQt5 (robust)")
        self._worker_thread = None
        self._worker = None
        self._build_ui()

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        layout = QGridLayout()
        central.setLayout(layout)

        # Number of qubits
        layout.addWidget(QLabel("Number of qubits:"), 0, 0)
        self.num_qubits_spin = QSpinBox()
        self.num_qubits_spin.setRange(2, 1024)
        self.num_qubits_spin.setValue(DEFAULT_NUM_QUBITS)
        self.num_qubits_spin.setSingleStep(2)
        layout.addWidget(self.num_qubits_spin, 0, 1)

        # Error rate
        layout.addWidget(QLabel("Error rate (0.0 - 1.0):"), 1, 0)
        self.error_spin = QDoubleSpinBox()
        self.error_spin.setRange(0.0, 1.0)
        self.error_spin.setSingleStep(0.01)
        self.error_spin.setDecimals(3)
        self.error_spin.setValue(0.075)
        layout.addWidget(self.error_spin, 1, 1)

        # Iterations
        layout.addWidget(QLabel("Iterations:"), 2, 0)
        self.iter_spin = QSpinBox()
        self.iter_spin.setRange(1, 100)
        self.iter_spin.setValue(5)
        layout.addWidget(self.iter_spin, 2, 1)

        # Mode
        layout.addWidget(QLabel("Mode:"), 3, 0)
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["normal", "eavesdrop"])
        layout.addWidget(self.mode_combo, 3, 1)

        # Buttons
        self.run_btn = QPushButton("Run Simulation")
        self.run_btn.clicked.connect(self.start_simulation)
        layout.addWidget(self.run_btn, 4, 0)

        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self.stop_simulation)
        layout.addWidget(self.stop_btn, 4, 1)

        # Output text area
        self.output = QTextEdit()
        self.output.setReadOnly(True)
        layout.addWidget(self.output, 5, 0, 1, 2)

    def append_output(self, text: str):
        self.output.append(text)
        # scroll to end
        self.output.verticalScrollBar().setValue(self.output.verticalScrollBar().maximum())

    def start_simulation(self):
        if self._worker_thread is not None:
            QMessageBox.warning(self, "Simulation running", "A simulation is already running.")
            return

        num_qubits = int(self.num_qubits_spin.value())
        error_rate = float(self.error_spin.value())
        iterations = int(self.iter_spin.value())
        mode = self.mode_combo.currentText()

        self.output.clear()
        self.append_output("Starting simulation...")
        # show qiskit status
        if USE_QISKIT:
            self.append_output(f"Qiskit detected. execute={USE_EXECUTE}, AerSim={USE_AER_SIM}")
        else:
            self.append_output("Qiskit/Aer not available — using fallback emulator (classical).")

        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)

        # Create worker + thread
        self._worker = SimulatorWorker(num_qubits=num_qubits, error_rate=error_rate, iterations=iterations, mode=mode)
        self._worker_thread = QThread()
        self._worker.moveToThread(self._worker_thread)
        self._worker_thread.started.connect(self._worker.run)
        self._worker.progress.connect(self.append_output)
        self._worker.finished.connect(self._on_finished)
        self._worker.error.connect(self._on_error)

        # ensure thread is quit when worker finished; keep reference to thread to wait later
        self._worker.finished.connect(self._worker_thread.quit)

        # cleanup connections when thread finishes
        self._worker.finished.connect(self._worker.deleteLater)
        self._worker_thread.finished.connect(self._worker_thread.deleteLater)

        self._worker_thread.start()

    def stop_simulation(self):
        """
        Request worker to stop, then ask thread to quit and wait a bit for it to finish.
        """
        if not self._worker_thread:
            return

        self.append_output("Stop requested...")
        try:
            if self._worker:
                self._worker.stop()  # set flag; worker checks flag between iterations
        except Exception:
            pass

        try:
            # ask thread to quit (if it's running an event loop)
            if self._worker_thread.isRunning():
                self._worker_thread.quit()
            # wait for thread to finish (short timeout to keep UI responsive)
            # but if backend.run is blocking, this will wait until it finishes
            self._worker_thread.wait(2000)  # wait up to 2 seconds
        except Exception:
            pass

        # reset UI regardless; finished handler will finalize if thread terminates later
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _on_finished(self):
        """
        Called when worker emits finished. Ensure thread has stopped before cleaning up.
        """
        self.append_output("Simulation finished (worker signalled).")
        try:
            if self._worker_thread and self._worker_thread.isRunning():
                # give it a chance to quit cleanly
                self._worker_thread.quit()
                self._worker_thread.wait(2000)
        except Exception:
            pass

        # reset state
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self._worker = None
        self._worker_thread = None

    def _on_error(self, tb: str):
        QMessageBox.critical(self, "Error in simulation", "An error occurred:\n" + tb)
        self.append_output("Error occurred (see dialog).")
        # ensure thread shutdown
        try:
            if self._worker_thread and self._worker_thread.isRunning():
                self._worker_thread.quit()
                self._worker_thread.wait(2000)
        except Exception:
            pass
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def closeEvent(self, event):
        """
        Ensure worker/thread are stopped when window is closed.
        """
        if self._worker_thread and self._worker_thread.isRunning():
            self.append_output("Window closing: stopping simulation...")
            try:
                if self._worker:
                    self._worker.stop()
            except Exception:
                pass
            try:
                self._worker_thread.quit()
                # wait a bit for the thread to finish; if it doesn't, proceed to close anyway
                self._worker_thread.wait(2000)
            except Exception:
                pass
        # proceed with normal close
        event.accept()


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.resize(800, 600)
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
