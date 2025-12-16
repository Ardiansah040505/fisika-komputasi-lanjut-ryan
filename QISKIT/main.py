from qiskit import QuantumCircuit
from qiskit.primitives import StatevectorSampler
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt

qc = QuantumCircuit(3)
qc.h(0)
qc.cx(0,1)
qc.cx(1,2)
qc.measure_all()

from qiskit import QuantumCircuit, Aer, transpile
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt

qc = QuantumCircuit(3)
qc.h(0)
qc.cx(0,1)
qc.cx(1,2)
qc.measure_all()

# Gunakan Aer simulator untuk menjalankan sirkuit dan mendapatkan counts
sim = Aer.get_backend('aer_simulator')
tcirc = transpile(qc, sim)
job = sim.run(tcirc, shots=1024)
result = job.result()
counts = result.get_counts()

plot_histogram(counts)
from qiskit import QuantumCircuit, Aer, transpile
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt

qc = QuantumCircuit(3)
qc.h(0)
qc.cx(0,1)
qc.cx(1,2)
qc.measure_all()

# Gunakan Aer simulator untuk menjalankan sirkuit dan mendapatkan counts
sim = Aer.get_backend('aer_simulator')
tcirc = transpile(qc, sim)
job = sim.run(tcirc, shots=1024)
result = job.result()
counts = result.get_counts()

plot_histogram(counts)
print(counts)
plt.show()