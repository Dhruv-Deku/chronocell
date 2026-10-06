"""
Quantum lab (experimental): ChronoCell problems written in the forms quantum computers take, solved on a
quantum SIMULATOR on this computer and compared with classical methods on the same input.

sim       statevector simulator (little-endian like Qiskit: qubit 0 is the lowest bit of a basis index),
          circuits with OpenQASM 2.0 export, QAOA, sampling, an approximate hardware-noise model.
qubo      QUBO / Ising problems; exact (enumeration), simulated annealing and simulated quantum annealing
          (path-integral Monte Carlo, a classical algorithm that imitates a quantum annealer).
problems  TAD boundaries, variant set, drug combination and gene group as QUBOs.
lattice   lattice folding of a short chain to match contacts (a higher-order diagonal Hamiltonian).
chem      minimal-basis (STO-3G) quantum chemistry from scratch, Jordan-Wigner qubits, VQE.
kernels   quantum feature-map kernels (QSVM), amplitude encoding and the swap test.
walk      continuous-time quantum walk on a contact graph, against the classical random walk.

Nothing here runs on quantum hardware, and nothing claims a speed-up: every result is labelled "simulated
quantum" and is shown next to the classical answer. No module imports Streamlit.
"""

LABEL = "simulated quantum (statevector simulator on this computer, not quantum hardware)"
