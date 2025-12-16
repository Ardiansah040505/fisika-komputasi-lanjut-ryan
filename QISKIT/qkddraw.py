#!/usr/bin/env python3
"""
qkddraw_fixed.py

Generate and save PNG diagrams for BB84 stages:
 - Alice encoding
 - (conceptual) Eve intercept-resend (simulated outcomes -> reprepare)
 - Bob measurement

Requirements:
    pip install qiskit matplotlib numpy
"""

import os
from typing import List
import numpy as np
from qiskit import QuantumCircuit
import matplotlib.pyplot as plt

# --- Config ---
NUM_QUBITS = 8   # small so diagrams remain readable
OUTPUT_DIR = "qkd_diagrams"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# --- Utility generators ---
def random_bits(n: int) -> List[int]:
    return list(np.random.randint(2, size=n))


# --- Circuit builders ---
def build_alice_encoding(bits: List[int], basis: List[int]) -> QuantumCircuit:
    """
    Build a QuantumCircuit showing Alice's encoding of bits using basis.
    basis: 0 -> Z basis (|0>,|1>), 1 -> X basis (|+>,|->)
    """
    n = len(bits)
    qc = QuantumCircuit(n, n, name="Alice: encode")
    for i in range(n):
        if basis[i] == 0:  # Z basis
            if bits[i] == 1:
                qc.x(i)
        else:  # X basis
            if bits[i] == 0:
                qc.h(i)            # |+>
            else:
                qc.x(i)
                qc.h(i)            # |->
    qc.barrier()
    return qc


def build_eve_intercept_resend(bits: List[int], alice_basis: List[int], eve_basis: List[int]) -> QuantumCircuit:
    """
    Build a conceptual intercept-resend circuit for Eve.
    Implementation:
      - Draw Alice encoding first (visual continuity)
      - Apply H where Eve measures in X-basis, then measure into classical bits
      - Simulate measurement outcomes in Python:
          if eve_basis == alice_basis -> measured = alice_bit
          else -> measured = random bit (50/50)
      - Reset qubits and reprepare them according to measured outcomes:
          If Eve measured in Z-basis: prepare |1> if measured==1 (apply X)
          If Eve measured in X-basis: prepare |+> if measured==0 (apply H),
                                     prepare |-> if measured==1 (apply X then H)
    This produces a clear diagram without using conditional gates (.c_if).
    """
    n = len(bits)
    # Create circuit with n qubits and n classical bits (we'll use first n classical for Eve's measurement)
    qc = QuantumCircuit(n, n, name="Eve: intercept-resend")

    # Visual: Alice encoding (same gates as Alice)
    for i in range(n):
        if alice_basis[i] == 0:
            if bits[i] == 1:
                qc.x(i)
        else:
            if bits[i] == 0:
                qc.h(i)
            else:
                qc.x(i)
                qc.h(i)
    qc.barrier()

    # Apply H where Eve measures in X-basis (so measuring in Z simulates X-basis)
    for i in range(n):
        if eve_basis[i] == 1:
            qc.h(i)
    # Measure into classical bits 0..n-1
    qc.measure(range(n), range(n))
    qc.barrier()

    # Simulate measurement outcomes in Python so we can draw the reprepare gates deterministically:
    measured = []
    for i in range(n):
        if eve_basis[i] == alice_basis[i]:
            # measuring in same basis yields the original bit deterministically
            measured_bit = bits[i]
        else:
            # measurement in different basis is random
            measured_bit = int(np.random.randint(2))
        measured.append(measured_bit)

    # Reset qubits (visual representation of preparing fresh qubits for resend)
    for i in range(n):
        qc.reset(i)

    # Reprepare qubits according to measured results (simulate resend)
    for i in range(n):
        if eve_basis[i] == 0:
            # Z-basis resend: prepare |1> if measured == 1
            if measured[i] == 1:
                qc.x(i)
        else:
            # X-basis resend:
            # measured == 0 -> prepare |+>  (apply H)
            # measured == 1 -> prepare |->  (apply X then H)
            if measured[i] == 0:
                qc.h(i)
            else:
                qc.x(i)
                qc.h(i)
    qc.barrier()
    return qc


def build_bob_measurement(basis: List[int]) -> QuantumCircuit:
    """
    Build Bob's measurement stage as a circuit:
    Apply H where Bob measures in X basis, then measure to classical bits.
    """
    n = len(basis)
    qc = QuantumCircuit(n, n, name="Bob: measure")
    for i in range(n):
        if basis[i] == 1:
            qc.h(i)
    qc.measure(range(n), range(n))
    return qc


# --- Drawing helpers ---
def save_circuit_png(qc: QuantumCircuit, filename: str, dpi: int = 200, fold: int = 100):
    """
    Draw circuit using matplotlib backend and save to PNG.
    """
    fig = qc.draw(output="mpl", fold=fold)
    outpath = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


# --- Example main ---
def main():
    # generate random data (or you may supply your own)
    bits = random_bits(NUM_QUBITS)
    alice_basis = random_bits(NUM_QUBITS)
    bob_basis = random_bits(NUM_QUBITS)
    eve_basis = random_bits(NUM_QUBITS)

    print("Alice bits   :", bits)
    print("Alice basis  :", alice_basis)
    print("Bob basis    :", bob_basis)
    print("Eve basis    :", eve_basis)

    # Alice's encoding diagram
    qc_a = build_alice_encoding(bits, alice_basis)
    save_circuit_png(qc_a, "bb84_alice.png", fold=100)

    # Eve intercept-resend diagram (conceptual & deterministic for drawing)
    qc_e = build_eve_intercept_resend(bits, alice_basis, eve_basis)
    save_circuit_png(qc_e, "bb84_eve.png", fold=100)

    # Bob measurement diagram
    qc_b = build_bob_measurement(bob_basis)
    save_circuit_png(qc_b, "bb84_bob.png", fold=100)

    print("All diagrams saved in directory:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
