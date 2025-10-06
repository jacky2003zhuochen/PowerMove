import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, transpile
from qiskit_aer.aerprovider import AerSimulator as Aer
from qiskit_aer.noise import NoiseModel, depolarizing_error
from scipy.optimize import minimize
from qiskit.quantum_info import SparsePauliOp
from  qiskit.quantum_info import Statevector


def create_noise_model(p1, p2, T1, T2, gate_time): 
    noise_model = NoiseModel()
    

    single_qubit_error = depolarizing_error(p1, 1)
    noise_model.add_all_qubit_quantum_error(single_qubit_error, ["rx", "ry", "rz", "u3"])
    

    two_qubit_error = depolarizing_error(p2, 2)
    noise_model.add_all_qubit_quantum_error(two_qubit_error, ["cz"])
    

    thermal_error = depolarizing_error(1 - np.exp(-gate_time / T2), 2)
    noise_model.add_all_qubit_quantum_error(thermal_error, ["cz"])

    return noise_model

def create_vqa_circuit(theta):

    n = len(theta)

    qc = QuantumCircuit(n)
    

    for i in range(n):
        qc.rx(theta[i], i)
    

    for i in range(n - 1):
        qc.cz(i, i + 1)
    qc.save_statevector()
        
    return qc

def get_hamiltonian(n):
    pauli_list = []
    coeffs = []
    for i in range(n - 1):
        pauli_string = ["I"] * n
        pauli_string[i] = "Z"
        pauli_string[i + 1] = "Z"
        pauli_list.append("".join(pauli_string))
        coeffs.append(-1.0)
    return SparsePauliOp(pauli_list, coeffs)

def vqa_energy(result, hamiltonian):
    
    psi = Statevector(result).data

    expectation = psi @ hamiltonian.to_matrix() @ psi.conj().T
    
    return np.real(expectation)


def vqa_optimization(noise_model, backend, theta_init, hamiltonian, shots=500):
    energy_history = []
    last_result = [0.00]

    def objective(theta):
        result_index = 0
        r = 1
        for i in range(r):
            qc = create_vqa_circuit(theta)
            transpiled = transpile(qc, backend)
            
            job = backend.run(transpiled, shots=shots, noise_model=noise_model)

            result = vqa_energy(job.result().get_statevector(transpiled), hamiltonian)
            result_index += result

        last_result[0] = result_index / r
        
        return result_index / r
    
    def callback(X):

        energy_history.append(last_result[0])
        print(last_result[0])


    
    result = minimize(objective, theta_init, method='Nelder-Mead', \
                      options={'maxiter': 301}, callback=callback)
    

    print(result)
    num_iterations = result.nit
    return result.fun, num_iterations, energy_history


def vqa_energy_sim(n_qubits):

    noise_model0 = create_noise_model(0, 0.005, None, 1.5 * 10 ** 6, 200)
    noise_model1 = create_noise_model(0, 0.0075, None, 1.5 * 10 ** 6, 400)
    noise_model2 = create_noise_model(0, 0, None, 1.5 * 10 ** 6, 0)

    backend0 = Aer(noise_model=noise_model0,
                    basis_gates=["rx", "ry", "rz", "u3", "cz"], method="statevector")
    backend1 = Aer(noise_model=noise_model1,
                    basis_gates=["rx", "ry", "rz", "u3", "cz"], method="statevector")
    backend2 = Aer(noise_model=noise_model2,
                    basis_gates=["rx", "ry", "rz", "u3", "cz"], method="statevector")

    if n_qubits == 15:
        theta_init = [2.05324197, 1.73213234, 0.57041856, 0.52622191, 0.05020943, 0.03605259, 1.60804159, 3.13068123, 1.52368694, 2.58169054, 0.15731502, 2.91539992, 0.77207426, 2.03051643, 0.01661375]
    elif n_qubits == 10:
        theta_init = [1.14568191, 1.70545605, 0.20067731, 0.84539167, 0.58839702, 1.97085317, 0.20749827, 0.35483746, 2.40901398, 1.81153558]
    hamiltonian = get_hamiltonian(len(theta_init))

    noiseless_result, noiseless_iters, energy_history0 = vqa_optimization(noise_model0, backend0, theta_init, hamiltonian)
    noisy_result, noisy_iters, energy_history1 = vqa_optimization(noise_model1, backend1, theta_init, hamiltonian)
    noisy_result, noisy_iters, energy_history2 = vqa_optimization(noise_model2, backend2, theta_init, hamiltonian)

    labels = ["Noiseless", "Noisy"]
    values = [noiseless_result, noisy_result]
    iterations = [noiseless_iters, noisy_iters]

    plt.figure(figsize=(8,5))
    plt.bar(labels, values, color=["blue", "red"])
    plt.ylabel("Final Expectation Value")
    plt.title("VQA Convergence with Noise Effects")

    for i, v in enumerate(values):
        plt.text(i, v + 0.02, f"{v:.3f}", ha='center', fontsize=12)

    # with open("data/vqa_energy_" + str(n_qubits) + "qubits.txt", 'w') as file:
    #     file.write(str(energy_history0) + '\n')
    #     file.write(str(energy_history1) + '\n')
    #     file.write(str(energy_history2) + '\n')
    print(theta_init)

for n in [10]:
    vqa_energy_sim(n)
