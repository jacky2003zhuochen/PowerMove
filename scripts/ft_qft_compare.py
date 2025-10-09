import sys
import os

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

from Construct_Circuit import *
from PowerMove import *
from Enola import *
import random
import math
from qiskit.qasm3 import loads

# n = 18
Distance_List = [5, 10, 15, 20, 25]
I_List = range(10)

# mvqc_transfer_duration_list = []
# mvqc_move_duration_list = [] 
# mvqc_cir_fidelity_list = [] 
# mvqc_cir_fidelity_1q_gate_list = [] 
# mvqc_cir_fidelity_2q_gate_list = [] 
# mvqc_cir_fidelity_2q_gate_for_idle_list = [] 
# mvqc_cir_fidelity_atom_transfer_list = [] 
# mvqc_cir_fidelity_coherence_list = []
# mvqc_nstage_list = []

# enola_transfer_duration_list = []
# enola_move_duration_list = [] 
# enola_cir_fidelity_list = [] 
# enola_cir_fidelity_1q_gate_list = [] 
# enola_cir_fidelity_2q_gate_list = [] 
# enola_cir_fidelity_2q_gate_for_idle_list = [] 
# enola_cir_fidelity_atom_transfer_list = [] 
# enola_cir_fidelity_coherence_list = []
# enola_nstage_list = []

# no_storage_transfer_duration_list = [] 
# no_storage_move_duration_list = [] 
# no_storage_cir_fidelity_list = [] 
# no_storage_cir_fidelity_1q_gate_list = [] 
# no_storage_cir_fidelity_2q_gate_list = [] 
# no_storage_cir_fidelity_2q_gate_for_idle_list = [] 
# no_storage_cir_fidelity_atom_transfer_list = [] 
# no_storage_cir_fidelity_coherence_list = []
# no_storage_nstage_list = []

method_list = ['base', "move_split", "change_dest", "change_dest+move_split"]
for method in method_list:
    no_storage_transfer_duration_list = [] 
    no_storage_move_duration_list = [] 
    no_storage_cir_fidelity_list = [] 
    no_storage_cir_fidelity_1q_gate_list = [] 
    no_storage_cir_fidelity_2q_gate_list = [] 
    no_storage_cir_fidelity_2q_gate_for_idle_list = [] 
    no_storage_cir_fidelity_atom_transfer_list = [] 
    no_storage_cir_fidelity_coherence_list = []
    no_storage_nstage_list = []
    chain_length_list = []
    loop_num_list = []

    # for d in Distance_List:
    d=1
    for n in [6, 15, 30, 40, 50, 70, 80]:
    # for n in [6, 10, 20, 30, 40, 50, 60, 80, 100]:
    # for n in [10, 20, 30, 50]:
        index = random.choice(I_List)
        Row = math.ceil(math.sqrt(n))

        # circ = QuantumCircuit.from_qasm_file(f"benchmarks/vqe/vqe_n{n}.qasm")
        with open(f"benchmarks/qft/qftn{n}.qasm", "r") as f:
        # with open(f"benchmarks/bv/new_bv_n{n}.qasm", "r") as f:
        # with open(f"benchmarks/rca/rca_n{n}.qasm", "r") as f:
        # with open(f"benchmarks/vqe/vqe_n{n}.qasm", "r") as f:
            qasm_code = f.read()
        circ = loads(qasm_code)

        test_circuit = transpile(circ, basis_gates=["u1", "u2", "u3", "cz", "id"],  optimization_level=2)

        cz_blocks = get_cz_blocks(test_circuit)

        # mvqc_transfer_duration, mvqc_move_duration, mvqc_cir_fidelity, mvqc_cir_fidelity_1q_gate, mvqc_cir_fidelity_2q_gate, mvqc_cir_fidelity_2q_gate_for_idle, mvqc_cir_fidelity_atom_transfer, mvqc_cir_fidelity_coherence, mvqc_nstage = mvqc(cz_blocks, Row, n, True, d, 1)
        # mvqc_transfer_duration_list.append(mvqc_transfer_duration)
        # mvqc_move_duration_list.append(mvqc_move_duration)
        # mvqc_cir_fidelity_list.append(mvqc_cir_fidelity)
        # mvqc_cir_fidelity_1q_gate_list.append(mvqc_cir_fidelity_1q_gate)
        # mvqc_cir_fidelity_2q_gate_list.append(mvqc_cir_fidelity_2q_gate)
        # mvqc_cir_fidelity_2q_gate_for_idle_list.append(mvqc_cir_fidelity_2q_gate_for_idle)
        # mvqc_cir_fidelity_atom_transfer_list.append(mvqc_cir_fidelity_atom_transfer)
        # mvqc_cir_fidelity_coherence_list.append(mvqc_cir_fidelity_coherence)   
        # mvqc_nstage_list.append(mvqc_nstage)

        # enola_transfer_duration, enola_move_duration, enola_cir_fidelity, enola_cir_fidelity_1q_gate, enola_cir_fidelity_2q_gate, enola_cir_fidelity_2q_gate_for_idle, enola_cir_fidelity_atom_transfer, enola_cir_fidelity_coherence, enola_nstage = enola(cz_blocks, Row, n, d)
        # enola_transfer_duration_list.append(enola_transfer_duration)
        # enola_move_duration_list.append(enola_move_duration)
        # enola_cir_fidelity_list.append(enola_cir_fidelity)
        # enola_cir_fidelity_1q_gate_list.append(enola_cir_fidelity_1q_gate)
        # enola_cir_fidelity_2q_gate_list.append(enola_cir_fidelity_2q_gate)
        # enola_cir_fidelity_2q_gate_for_idle_list.append(enola_cir_fidelity_2q_gate_for_idle)
        # enola_cir_fidelity_atom_transfer_list.append(enola_cir_fidelity_atom_transfer)
        # enola_cir_fidelity_coherence_list.append(enola_cir_fidelity_coherence)  
        # enola_nstage_list.append(enola_nstage)

        no_storage_transfer_duration, no_storage_move_duration, no_storage_cir_fidelity, no_storage_cir_fidelity_1q_gate, no_storage_cir_fidelity_2q_gate, no_storage_cir_fidelity_2q_gate_for_idle, no_storage_cir_fidelity_atom_transfer, no_storage_cir_fidelity_coherence, no_storage_nstage, count, loop_num = mvqc(cz_blocks, Row, n, False, d, 1, method)

        sorted(count.items())
        print(loop_num)
        print(count)
        no_storage_transfer_duration_list.append(no_storage_transfer_duration)
        no_storage_move_duration_list.append(no_storage_move_duration)
        no_storage_cir_fidelity_list.append(no_storage_cir_fidelity)
        no_storage_cir_fidelity_1q_gate_list.append(no_storage_cir_fidelity_1q_gate)
        no_storage_cir_fidelity_2q_gate_list.append(no_storage_cir_fidelity_2q_gate)
        no_storage_cir_fidelity_2q_gate_for_idle_list.append(no_storage_cir_fidelity_2q_gate_for_idle)
        no_storage_cir_fidelity_atom_transfer_list.append(no_storage_cir_fidelity_atom_transfer)
        no_storage_cir_fidelity_coherence_list.append(no_storage_cir_fidelity_coherence)   
        no_storage_nstage_list.append(no_storage_nstage)
        chain_length_list.append(count)
        loop_num_list.append(loop_num)

    with open(f"data/new_compare_qft_{method}.txt", 'w') as file:
        # file.write(str(Distance_List) + '\n')
        # file.write(str(mvqc_transfer_duration_list) + '\n') 
        # file.write(str(mvqc_move_duration_list) + '\n') 
        # file.write(str([x+y for x,y in zip(mvqc_transfer_duration_list,mvqc_move_duration_list)]) + '\n') 
        # file.write(str(mvqc_cir_fidelity_list) + '\n') 
        # file.write(str(mvqc_cir_fidelity_1q_gate_list) + '\n') 
        # file.write(str(mvqc_cir_fidelity_2q_gate_list) + '\n') 
        # file.write(str(mvqc_cir_fidelity_2q_gate_for_idle_list) + '\n') 
        # file.write(str(mvqc_cir_fidelity_atom_transfer_list) + '\n') 
        # file.write(str(mvqc_cir_fidelity_coherence_list) + '\n')
        # file.write(str(mvqc_nstage_list) + '\n')

        # file.write(str(enola_transfer_duration_list) + '\n') 
        # file.write(str(enola_move_duration_list) + '\n') 
        # file.write(str(enola_cir_fidelity_list) + '\n') 
        # file.write(str(enola_cir_fidelity_1q_gate_list) + '\n') 
        # file.write(str(enola_cir_fidelity_2q_gate_list) + '\n') 
        # file.write(str(enola_cir_fidelity_2q_gate_for_idle_list) + '\n') 
        # file.write(str(enola_cir_fidelity_atom_transfer_list) + '\n') 
        # file.write(str(enola_cir_fidelity_coherence_list) + '\n')
        # file.write(str(enola_nstage_list) + '\n')
        
        file.write(str([x + y for x, y in zip(no_storage_transfer_duration_list, no_storage_move_duration_list)]) + '\n') 
        file.write(str(loop_num_list) + '\n')  
        file.write(str(chain_length_list) + '\n')  
        file.write(str(no_storage_transfer_duration_list) + '\n')  
        file.write(str(no_storage_move_duration_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_1q_gate_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_2q_gate_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_2q_gate_for_idle_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_atom_transfer_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_coherence_list) + '\n') 
        file.write(str(no_storage_nstage_list) + '\n') 

