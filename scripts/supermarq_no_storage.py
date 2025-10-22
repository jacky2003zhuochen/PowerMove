import sys
import os
from pathlib import Path
import supermarq

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

from Construct_Circuit import *
from PowerMove import *
from Enola import *
import random
import math
from qiskit.qasm3 import loads

# n = 18
# Distance_List = [5, 10, 15, 20, 25]
Distance_List = [5, 10, 15, 20, 25]
I_List = range(10)

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

P_List = [0.1, 0.2, 0.3, 0.4, 0.5]
cost_para_list = {'ghz':0.8, 'cat':0.8, 'ising':0.7, 'wstate':0.5, 'qft':0.9}

bench_list = ['hs','mb']
method_list = ['base', 'move_split', 'change_dest', 'change_dest+move_split', "break_chains", "break_chains+change_dest", "break_chains+change_dest+move_split"]
# method_list = ['move_split', 'change_dest+move_split', "break_chains+change_dest+move_split"]
# method_list = ['break_chains+change_dest+move_split']
# method_list = ['move_split_double',  'change_dest+move_split_double', "break_chains+change_dest+move_split_double"]
# method_list = ['break_chains', "break_chains+change_dest", "break_chains+change_dest+move_split"]
para2 = 0
for method in method_list:
    print(method)

    # for d in Distance_List:
    d=1

    # for thre in [0.1,0.5,1,2,5,10,20]:
    # for thre in [0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1]:
    # for thre in [60,90]:
    # for thre in [-10,-9,-8,-7,-6,-5,-4,-3,-1]:
    # for thre in [-9,-8,-7,-5,-2,0,2,4]:

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
    threshold_length_list = []
    loop_num_list = []
    q_list = []
    for n in [10, 30, 50, 70, 100, 150, 200]:
        print(n)
        index = random.choice(I_List)
        Row = math.ceil(math.sqrt(n))

        hs = supermarq.benchmarks.hamiltonian_simulation.HamiltonianSimulation(num_qubits=n)
        hs_circuit = hs.qiskit_circuit()
        
        test_circuit_hs = transpile(hs_circuit, basis_gates=["u1", "u2", "u3", "cz", "id"],  optimization_level=2)
        cz_blocks_hs = get_cz_blocks(test_circuit_hs)

        no_storage_transfer_duration, no_storage_move_duration, no_storage_cir_fidelity, no_storage_cir_fidelity_1q_gate, no_storage_cir_fidelity_2q_gate, no_storage_cir_fidelity_2q_gate_for_idle, no_storage_cir_fidelity_atom_transfer, no_storage_cir_fidelity_coherence, no_storage_nstage, count, loop_num, split_succ, split_fail = mvqc(cz_blocks_hs, Row, n, False, d, 1, method, cost_para=100, para1=0, para2=0)

        sorted(count.items())
        # print(loop_num)
        # print(count)
        count = dict(sorted(count.items()))
        threshold_length = find_threshold_key(count, 0.7)
        threshold_length_list.append(threshold_length)
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
    # print(split_succ, split_fail)

    with open(f"data/compare_hs_{method}.txt", 'w') as file:
    # with open(f"data/compare_{benchm}_{method}_sort_rand.txt", 'w') as file:
    # with open(f"data/compare_{benchm}_{method}.txt", 'w') as file:
        file.write(str([x + y for x, y in zip(no_storage_transfer_duration_list, no_storage_move_duration_list)]) + '\n') 
        file.write(str(loop_num_list) + '\n')  
        file.write(str(chain_length_list) + '\n')  
        file.write(str(threshold_length_list) + '\n')
        file.write(str(no_storage_transfer_duration_list) + '\n')  
        file.write(str(no_storage_move_duration_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_1q_gate_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_2q_gate_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_2q_gate_for_idle_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_atom_transfer_list) + '\n')  
        file.write(str(no_storage_cir_fidelity_coherence_list) + '\n') 
        file.write(str(no_storage_nstage_list) + '\n') 

        # no_storage_transfer_duration_list = [] 
        # no_storage_move_duration_list = [] 
        # no_storage_cir_fidelity_list = [] 
        # no_storage_cir_fidelity_1q_gate_list = [] 
        # no_storage_cir_fidelity_2q_gate_list = [] 
        # no_storage_cir_fidelity_2q_gate_for_idle_list = [] 
        # no_storage_cir_fidelity_atom_transfer_list = [] 
        # no_storage_cir_fidelity_coherence_list = []
        # no_storage_nstage_list = []
        # chain_length_list = []
        # threshold_length_list = []
        # loop_num_list = []
        # q_list = []

        # for n in [10, 30, 50]:
        #     print(n)
        #     index = random.choice(I_List)
        #     Row = math.ceil(math.sqrt(n))
        #     mb = supermarq.benchmarks.mermin_bell.MerminBell(num_qubits=n)
        #     mb_circuit = mb.qiskit_circuit()
        #     test_circuit_mb = transpile(mb_circuit, basis_gates=["u1", "u2", "u3", "cz", "id"],  optimization_level=2)
        #     cz_blocks_mb = get_cz_blocks(test_circuit_mb)

        #     no_storage_transfer_duration, no_storage_move_duration, no_storage_cir_fidelity, no_storage_cir_fidelity_1q_gate, no_storage_cir_fidelity_2q_gate, no_storage_cir_fidelity_2q_gate_for_idle, no_storage_cir_fidelity_atom_transfer, no_storage_cir_fidelity_coherence, no_storage_nstage, count, loop_num, split_succ, split_fail = mvqc(cz_blocks_mb, Row, n, False, d, 1, method, cost_para=thre, para1=0, para2=0)
        #     sorted(count.items())
        #     # print(loop_num)
        #     # print(count)
        #     count = dict(sorted(count.items()))
        #     threshold_length = find_threshold_key(count, 0.7)
        #     threshold_length_list.append(threshold_length)
        #     no_storage_transfer_duration_list.append(no_storage_transfer_duration)
        #     no_storage_move_duration_list.append(no_storage_move_duration)
        #     no_storage_cir_fidelity_list.append(no_storage_cir_fidelity)
        #     no_storage_cir_fidelity_1q_gate_list.append(no_storage_cir_fidelity_1q_gate)
        #     no_storage_cir_fidelity_2q_gate_list.append(no_storage_cir_fidelity_2q_gate)
        #     no_storage_cir_fidelity_2q_gate_for_idle_list.append(no_storage_cir_fidelity_2q_gate_for_idle)
        #     no_storage_cir_fidelity_atom_transfer_list.append(no_storage_cir_fidelity_atom_transfer)
        #     no_storage_cir_fidelity_coherence_list.append(no_storage_cir_fidelity_coherence)   
        #     no_storage_nstage_list.append(no_storage_nstage)
        #     chain_length_list.append(count)
        #     loop_num_list.append(loop_num)
            
        # with open(f"data/compare_mb_{method}_sim_{thre}.txt", 'w') as file:
        # # with open(f"data/compare_{benchm}_{method}_sort_rand.txt", 'w') as file:
        # # with open(f"data/compare_{benchm}_{method}.txt", 'w') as file:
        #     file.write(str([x + y for x, y in zip(no_storage_transfer_duration_list, no_storage_move_duration_list)]) + '\n') 
        #     file.write(str(loop_num_list) + '\n')  
        #     file.write(str(chain_length_list) + '\n')  
        #     file.write(str(threshold_length_list) + '\n')
        #     file.write(str(no_storage_transfer_duration_list) + '\n')  
        #     file.write(str(no_storage_move_duration_list) + '\n')  
        #     file.write(str(no_storage_cir_fidelity_list) + '\n')  
        #     file.write(str(no_storage_cir_fidelity_1q_gate_list) + '\n')  
        #     file.write(str(no_storage_cir_fidelity_2q_gate_list) + '\n')  
        #     file.write(str(no_storage_cir_fidelity_2q_gate_for_idle_list) + '\n')  
        #     file.write(str(no_storage_cir_fidelity_atom_transfer_list) + '\n')  
        #     file.write(str(no_storage_cir_fidelity_coherence_list) + '\n') 
        #     file.write(str(no_storage_nstage_list) + '\n') 
