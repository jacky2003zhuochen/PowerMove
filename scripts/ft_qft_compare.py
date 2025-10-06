import sys
import os

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

from Construct_Circuit import *
from PowerMove import *
from Enola import *
import random
import math

n = 18
Distance_List = [5, 10, 15, 20, 25]
I_List = range(10)

mvqc_transfer_duration_list = []
mvqc_move_duration_list = [] 
mvqc_cir_fidelity_list = [] 
mvqc_cir_fidelity_1q_gate_list = [] 
mvqc_cir_fidelity_2q_gate_list = [] 
mvqc_cir_fidelity_2q_gate_for_idle_list = [] 
mvqc_cir_fidelity_atom_transfer_list = [] 
mvqc_cir_fidelity_coherence_list = []
mvqc_nstage_list = []

enola_transfer_duration_list = []
enola_move_duration_list = [] 
enola_cir_fidelity_list = [] 
enola_cir_fidelity_1q_gate_list = [] 
enola_cir_fidelity_2q_gate_list = [] 
enola_cir_fidelity_2q_gate_for_idle_list = [] 
enola_cir_fidelity_atom_transfer_list = [] 
enola_cir_fidelity_coherence_list = []
enola_nstage_list = []

no_storage_transfer_duration_list = [] 
no_storage_move_duration_list = [] 
no_storage_cir_fidelity_list = [] 
no_storage_cir_fidelity_1q_gate_list = [] 
no_storage_cir_fidelity_2q_gate_list = [] 
no_storage_cir_fidelity_2q_gate_for_idle_list = [] 
no_storage_cir_fidelity_atom_transfer_list = [] 
no_storage_cir_fidelity_coherence_list = []
no_storage_nstage_list = []

for d in Distance_List:
    index = random.choice(I_List)
    Row = math.ceil(math.sqrt(n))

    circ = QuantumCircuit.from_qasm_file(f"benchmarks/qft/qft_n{n}.qasm")

    test_circuit = transpile(circ, basis_gates=["u1", "u2", "u3", "cz", "id"],  optimization_level=2)

    cz_blocks = get_cz_blocks(test_circuit)

    mvqc_transfer_duration, mvqc_move_duration, mvqc_cir_fidelity, mvqc_cir_fidelity_1q_gate, mvqc_cir_fidelity_2q_gate, mvqc_cir_fidelity_2q_gate_for_idle, mvqc_cir_fidelity_atom_transfer, mvqc_cir_fidelity_coherence, mvqc_nstage = mvqc(cz_blocks, Row, n, True, d, 1)
    mvqc_transfer_duration_list.append(mvqc_transfer_duration)
    mvqc_move_duration_list.append(mvqc_move_duration)
    mvqc_cir_fidelity_list.append(mvqc_cir_fidelity)
    mvqc_cir_fidelity_1q_gate_list.append(mvqc_cir_fidelity_1q_gate)
    mvqc_cir_fidelity_2q_gate_list.append(mvqc_cir_fidelity_2q_gate)
    mvqc_cir_fidelity_2q_gate_for_idle_list.append(mvqc_cir_fidelity_2q_gate_for_idle)
    mvqc_cir_fidelity_atom_transfer_list.append(mvqc_cir_fidelity_atom_transfer)
    mvqc_cir_fidelity_coherence_list.append(mvqc_cir_fidelity_coherence)   
    mvqc_nstage_list.append(mvqc_nstage)

    enola_transfer_duration, enola_move_duration, enola_cir_fidelity, enola_cir_fidelity_1q_gate, enola_cir_fidelity_2q_gate, enola_cir_fidelity_2q_gate_for_idle, enola_cir_fidelity_atom_transfer, enola_cir_fidelity_coherence, enola_nstage = enola(cz_blocks, Row, n, d)
    enola_transfer_duration_list.append(enola_transfer_duration)
    enola_move_duration_list.append(enola_move_duration)
    enola_cir_fidelity_list.append(enola_cir_fidelity)
    enola_cir_fidelity_1q_gate_list.append(enola_cir_fidelity_1q_gate)
    enola_cir_fidelity_2q_gate_list.append(enola_cir_fidelity_2q_gate)
    enola_cir_fidelity_2q_gate_for_idle_list.append(enola_cir_fidelity_2q_gate_for_idle)
    enola_cir_fidelity_atom_transfer_list.append(enola_cir_fidelity_atom_transfer)
    enola_cir_fidelity_coherence_list.append(enola_cir_fidelity_coherence)  
    enola_nstage_list.append(enola_nstage)

    no_storage_transfer_duration, no_storage_move_duration, no_storage_cir_fidelity, no_storage_cir_fidelity_1q_gate, no_storage_cir_fidelity_2q_gate, no_storage_cir_fidelity_2q_gate_for_idle, no_storage_cir_fidelity_atom_transfer, no_storage_cir_fidelity_coherence, no_storage_nstage = mvqc(cz_blocks, Row, n, False, d, 1)
    no_storage_transfer_duration_list.append(no_storage_transfer_duration)
    no_storage_move_duration_list.append(no_storage_move_duration)
    no_storage_cir_fidelity_list.append(no_storage_cir_fidelity)
    no_storage_cir_fidelity_1q_gate_list.append(no_storage_cir_fidelity_1q_gate)
    no_storage_cir_fidelity_2q_gate_list.append(no_storage_cir_fidelity_2q_gate)
    no_storage_cir_fidelity_2q_gate_for_idle_list.append(no_storage_cir_fidelity_2q_gate_for_idle)
    no_storage_cir_fidelity_atom_transfer_list.append(no_storage_cir_fidelity_atom_transfer)
    no_storage_cir_fidelity_coherence_list.append(no_storage_cir_fidelity_coherence)   
    no_storage_nstage_list.append(no_storage_nstage)

with open("data/fault_tolerant_compare_qft.txt", 'w') as file:
    file.write(str(Distance_List) + '\n')
    file.write(str(mvqc_transfer_duration_list) + '\n') 
    file.write(str(mvqc_move_duration_list) + '\n') 
    file.write(str(mvqc_cir_fidelity_list) + '\n') 
    file.write(str(mvqc_cir_fidelity_1q_gate_list) + '\n') 
    file.write(str(mvqc_cir_fidelity_2q_gate_list) + '\n') 
    file.write(str(mvqc_cir_fidelity_2q_gate_for_idle_list) + '\n') 
    file.write(str(mvqc_cir_fidelity_atom_transfer_list) + '\n') 
    file.write(str(mvqc_cir_fidelity_coherence_list) + '\n')
    file.write(str(mvqc_nstage_list) + '\n')

    file.write(str(enola_transfer_duration_list) + '\n') 
    file.write(str(enola_move_duration_list) + '\n') 
    file.write(str(enola_cir_fidelity_list) + '\n') 
    file.write(str(enola_cir_fidelity_1q_gate_list) + '\n') 
    file.write(str(enola_cir_fidelity_2q_gate_list) + '\n') 
    file.write(str(enola_cir_fidelity_2q_gate_for_idle_list) + '\n') 
    file.write(str(enola_cir_fidelity_atom_transfer_list) + '\n') 
    file.write(str(enola_cir_fidelity_coherence_list) + '\n')
    file.write(str(enola_nstage_list) + '\n')

    file.write(str(no_storage_transfer_duration_list) + '\n')  
    file.write(str(no_storage_move_duration_list) + '\n')  
    file.write(str(no_storage_cir_fidelity_list) + '\n')  
    file.write(str(no_storage_cir_fidelity_1q_gate_list) + '\n')  
    file.write(str(no_storage_cir_fidelity_2q_gate_list) + '\n')  
    file.write(str(no_storage_cir_fidelity_2q_gate_for_idle_list) + '\n')  
    file.write(str(no_storage_cir_fidelity_atom_transfer_list) + '\n')  
    file.write(str(no_storage_cir_fidelity_coherence_list) + '\n') 
    file.write(str(no_storage_nstage_list) + '\n') 

