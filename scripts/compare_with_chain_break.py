import sys
import os

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

from PowerMove import *
from Enola import *
import time
import random
import math
import sys
import os
from Construct_Circuit import *
from PowerMove import *
from qiskit.qasm3 import loads

def pauli_string_to_qasm(pauli_str_list, qreg_name="q"):
    qasm_lines = []
    n = len(pauli_str_list[0])
    header = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n"
    reg = f"qreg q[{n}];\n"
    qasm_lines.append(header)
    qasm_lines.append(reg)
    for pauli_str in pauli_str_list:
        for idx, op in enumerate(pauli_str):
            if op == 'I':
                continue
            elif op in ['X', 'Y', 'Z']:
                qasm_lines.append(f"{op.lower()} {qreg_name}[{idx}];")
            else:
                raise ValueError(f"Invalid character '{op}' in Pauli string.")
    return "\n".join(qasm_lines)

def find_threshold_key(d, threshold=0.7):
    sorted_items = sorted(d.items()) 
    
    total = sum(v for _, v in sorted_items)
    target = threshold * total
    
    cumulative = 0
    for key, value in sorted_items:
        if cumulative > target:
            return key
        cumulative += value

    return sorted_items[-1][0] if sorted_items else None


# N_Qubit_List = [5, 10, 20, 30 ,50, 100]
N_Qubit_List = [6, 10, 20, 30, 40, 50, 60, 80, 100]
P_List = [0.1, 0.2, 0.3, 0.4, 0.5]
# P_List = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
I_List = range(10)

index = random.choice(I_List)
type = 'rand'
# type = 'regular'
d = 1
method_list = ['base', 'move_split', 'change_dest', 'change_dest+move_split', "break_chains", "break_chains+change_dest", "break_chains+change_dest+move_split"]
# method_list = ["break_chains", "break_chains+change_dest", "break_chains+change_dest+move_split"]
for method in method_list:
    print(method)
    for P in P_List:
    # for degree in range(1,11):
    # for P in [0.5]:
        mvqc_transfer_duration_list = []
        mvqc_move_duration_list = [] 
        mvqc_cir_fidelity_list = [] 
        mvqc_cir_fidelity_1q_gate_list = [] 
        mvqc_cir_fidelity_2q_gate_list = [] 
        mvqc_cir_fidelity_2q_gate_for_idle_list = [] 
        mvqc_cir_fidelity_atom_transfer_list = [] 
        mvqc_cir_fidelity_coherence_list = []
        mvqc_nstage_list = []
        mvqc_runtime = []
        chain_length_list = []
        loop_num_list = []
        threshold_length_list = []
        q_list = []

        for n in N_Qubit_List:
            Row = math.ceil(math.sqrt(n))

            path = f"benchmarks/qaoa/{type}/q{n}_p{P}/i{index}.txt"
            # path = f"benchmarks/qaoa/{type}/q{n}_regular{degree}/i{index}.txt"
            # path = f"benchmarks/qsim/{type}/q{n}_10_p{P}/i{index}.txt"
            try:
                with open(path, "r") as fid:
                    gates = eval(fid.read())
                    q_list.append(n)
            except FileNotFoundError:
                continue
            with open(path, "r") as fid:
                gates = eval(fid.read())
            # print(gates)
            # circuit = pauli_string_to_qasm(gates)
            # print(circuit)
            # cir = QuantumCircuit.from_qasm_str(circuit)
            # print(cir)
            # mvqc_start_time = time.time()
            mvqc_transfer_duration, mvqc_move_duration, mvqc_cir_fidelity, mvqc_cir_fidelity_1q_gate, mvqc_cir_fidelity_2q_gate, mvqc_cir_fidelity_2q_gate_for_idle, mvqc_cir_fidelity_atom_transfer, mvqc_cir_fidelity_coherence, mvqc_nstage, count, loop_num = mvqc([gates], Row, n, False, d, 1, method)
            
            # sorted(count.items())
            # print("loop num", loop_num)
            # print("chains length", count)
            count = dict(sorted(count.items()))
            threshold_length = find_threshold_key(count, 0.7)
            # print(threshold_length)
            threshold_length_list.append(threshold_length)
            mvqc_transfer_duration_list.append(mvqc_transfer_duration)
            mvqc_move_duration_list.append(mvqc_move_duration)
            mvqc_cir_fidelity_list.append(mvqc_cir_fidelity)
            mvqc_cir_fidelity_1q_gate_list.append(mvqc_cir_fidelity_1q_gate)
            mvqc_cir_fidelity_2q_gate_list.append(mvqc_cir_fidelity_2q_gate)
            mvqc_cir_fidelity_2q_gate_for_idle_list.append(mvqc_cir_fidelity_2q_gate_for_idle)
            mvqc_cir_fidelity_atom_transfer_list.append(mvqc_cir_fidelity_atom_transfer)
            mvqc_cir_fidelity_coherence_list.append(mvqc_cir_fidelity_coherence)   
            mvqc_nstage_list.append(mvqc_nstage)
            chain_length_list.append(count)
            loop_num_list.append(loop_num)

            # enola_start_time = time.time()
            # enolar_transfer_duration, enola_move_duration, enola_cir_fidelity, enola_cir_fidelity_1q_gate, enola_cir_fidelity_2q_gate, enola_cir_fidelity_2q_gate_for_idle, enola_cir_fidelity_atom_transfer, enola_cir_fidelity_coherence, enola_nstage = enola([gates], Row, n, d)
            # enola_runtime.append(time.time() - enola_start_time)
            
            # enola_transfer_duration_list.append(enolar_transfer_duration)
            # enola_move_duration_list.append(enola_move_duration)
            # enola_cir_fidelity_list.append(enola_cir_fidelity)
            # enola_cir_fidelity_1q_gate_list.append(enola_cir_fidelity_1q_gate)
            # enola_cir_fidelity_2q_gate_list.append(enola_cir_fidelity_2q_gate)
            # enola_cir_fidelity_2q_gate_for_idle_list.append(enola_cir_fidelity_2q_gate_for_idle)
            # enola_cir_fidelity_atom_transfer_list.append(enola_cir_fidelity_atom_transfer)
            # enola_cir_fidelity_coherence_list.append(enola_cir_fidelity_coherence)  
            # enola_nstage_list.append(enola_nstage)
        with open(f"data/qaoa_{type}_no_storage{P}_compare_{method}.txt", 'w') as file:
        # with open(f"data/qaoa_{type}{degree}_no_storage_compare_{method}.txt", 'w') as file:
            # file.write(str(N_Qubit_List) + '\n')
            # file.write(str(q_list) + '\n')        
            file.write(str([x + y for x, y in zip(mvqc_transfer_duration_list, mvqc_move_duration_list)]) + '\n') 
            file.write(str(loop_num_list) + '\n') 
            file.write(str(chain_length_list) + '\n') 
            file.write(str(threshold_length_list) + '\n')
            file.write(str(mvqc_transfer_duration_list) + '\n') 
            file.write(str(mvqc_move_duration_list) + '\n') 
            file.write(str(mvqc_cir_fidelity_list) + '\n') 
            file.write(str(mvqc_cir_fidelity_1q_gate_list) + '\n') 
            file.write(str(mvqc_cir_fidelity_2q_gate_list) + '\n') 
            file.write(str(mvqc_cir_fidelity_2q_gate_for_idle_list) + '\n') 
            file.write(str(mvqc_cir_fidelity_atom_transfer_list) + '\n') 
            file.write(str(mvqc_cir_fidelity_coherence_list) + '\n')
            file.write(str(mvqc_nstage_list) + '\n')
            file.write(str(mvqc_runtime) + '\n')

            # file.write(str(enola_transfer_duration_list) + '\n') 
            # file.write(str(enola_move_duration_list) + '\n') 
            # file.write(str(enola_cir_fidelity_list) + '\n') 
            # file.write(str(enola_cir_fidelity_1q_gate_list) + '\n') 
            # file.write(str(enola_cir_fidelity_2q_gate_list) + '\n') 
            # file.write(str(enola_cir_fidelity_2q_gate_for_idle_list) + '\n') 
            # file.write(str(enola_cir_fidelity_atom_transfer_list) + '\n') 
            # file.write(str(enola_cir_fidelity_coherence_list) + '\n')
            # file.write(str(enola_nstage_list) + '\n')
            # file.write(str(enola_runtime) + '\n')
