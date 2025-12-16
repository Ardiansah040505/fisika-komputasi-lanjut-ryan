from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from numpy.random import randint
import numpy as np
from art import *

NUM_QUBITS = 32


def generate_bits(n: int) -> np.ndarray:
    return randint(2, size=n)


def encode_message(bits: np.ndarray, basis: np.ndarray) -> list:
    message = []
    for i in range(NUM_QUBITS):
        qc = QuantumCircuit(1, 1)
        if basis[i] == 0:
            if bits[i] == 1:
                qc.x(0)
        else:
            if bits[i] == 0:
                qc.h(0)
            else:
                qc.x(0)
                qc.h(0)
        qc.barrier()
        message.append(qc)
    return message


def simulate_quantum_channel(message: list, error_rate: float) -> list:
    noisy_message = []
    for qc in message:
        qc_copy = qc.copy()
        if np.random.random() < error_rate:
            qc_copy.x(0)
        noisy_message.append(qc_copy)
    return noisy_message


def measure_message(message: list, basis: np.ndarray) -> list:
    """
    Penggantian utama: jangan pakai Aer.get_backend atau execute.
    Gunakan AerSimulator().run(circuit, shots=1, memory=True).
    """
    sim = AerSimulator()
    measurements = []

    for q in range(NUM_QUBITS):
        qc = message[q].copy()

        if basis[q] == 1:  # measure in X-basis
            qc.h(0)

        qc.measure(0, 0)
        job = sim.run(qc, shots=1, memory=True)
        result = job.result()
        measurements.append(int(result.get_memory()[0]))

    return measurements


def remove_garbage(a_basis: np.ndarray, b_basis: np.ndarray, bits: np.ndarray) -> list:
    return [bits[q] for q in range(NUM_QUBITS) if a_basis[q] == b_basis[q]]


def check_keys(key1: list, key2: list) -> None:
    print("\nAlice's key: ", key1)
    print("Bob's key:   ", key2)
    if key1 == key2:
        print("Keys are the same and secure.")
    else:
        print("Error: keys are different.")


def quantum_communication(error_rate: float, num_iterations: int):
    for i in range(num_iterations):
        alice_bits = generate_bits(NUM_QUBITS)
        alice_basis = generate_bits(NUM_QUBITS)
        message = encode_message(alice_bits, alice_basis)
        message = simulate_quantum_channel(message, error_rate)
        bob_basis = generate_bits(NUM_QUBITS)
        bob_results = measure_message(message, bob_basis)
        alice_key = remove_garbage(alice_basis, bob_basis, alice_bits)
        bob_key = remove_garbage(alice_basis, bob_basis, bob_results)
        check_keys(alice_key, bob_key)


def simulate_eavesdropping(message: list, eavesdropper_basis: np.ndarray, error_rate: float) -> list:
    """
    Simulate intercept-resend. Setiap qubit di-measure, lalu dibuat ulang sesuai hasil
    agar Bob menerima 'resend' qubit.
    """
    sim = AerSimulator()
    intercepted = []

    for i in range(len(message)):
        qc = message[i].copy()

        if np.random.random() < 0.5:  # Z-basis
            qc.measure(0, 0)
            job = sim.run(qc, shots=1, memory=True)
            m = job.result().get_memory()[0]
            resend = QuantumCircuit(1, 1)
            if int(m) == 1:
                resend.x(0)
        else:  # X-basis
            qc.h(0)
            qc.measure(0, 0)
            job = sim.run(qc, shots=1, memory=True)
            m = job.result().get_memory()[0]
            resend = QuantumCircuit(1, 1)
            if int(m) == 0:
                resend.h(0)  # |+>
            else:
                resend.x(0)
                resend.h(0)  # |->
        if np.random.random() < error_rate:
            resend.x(0)
        intercepted.append(resend)

    return intercepted


def main():
    print("-" * 50 + "\n")
    print(text2art("BB84", font="small"))
    print("-" * 50)
    print("CTG Assignment CSF02 2023\n")
    print(
        "This script simulates QKD using the BB84 protocol.\n"
        "It includes functions for generating random bits,\n"
        "encoding messages with qubits,simulating a quantum\n"
        "channel with errors, measuring quantum messages,\nand removing bits measured in different bases.\n"
        "Users can choose between simulating communication\nwith no error, low error rate, high error rate,\n"
        "or simulating an eavesdropping attempt to observe\nthe effects on the integrity of the \n"
        "data transferred."
    )

    print("-" * 50 + "\n\n")

    num_iterations = 5
    while True:
        print("Choose an option:")
        print("[1] Simulate with 0 error")
        print("[2] Simulate with low error rate")
        print("[3] Simulate with high error rate")
        print("[4] Simulate eavesdropping attempt")
        print("[5] Exit")
        choice = input("Enter your choice (1/2/3/4/5): ")

        if choice == "1":
            quantum_communication(0.0, num_iterations)
            print()
        elif choice == "2":
            quantum_communication(0.075, num_iterations)
            print()
        elif choice == "3":
            quantum_communication(0.2, num_iterations)
            print()
        elif choice == "4":
            for _ in range(num_iterations):
                alice_bits = generate_bits(NUM_QUBITS)
                alice_basis = generate_bits(NUM_QUBITS)
                message = encode_message(alice_bits, alice_basis)
                eavesdropper_basis = generate_bits(NUM_QUBITS)
                message = simulate_eavesdropping(message, eavesdropper_basis, 0.1)
                bob_basis = generate_bits(NUM_QUBITS)
                bob_results = measure_message(message, bob_basis)
                alice_key = remove_garbage(alice_basis, bob_basis, alice_bits)
                bob_key = remove_garbage(alice_basis, bob_basis, bob_results)
                check_keys(alice_key, bob_key)
            print()
        elif choice == "5":
            break
        else:
            print("Invalid choice. Please enter a valid option.")


if __name__ == "__main__":
    main()
