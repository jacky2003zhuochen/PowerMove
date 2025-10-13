import sys
import os
from pathlib import Path

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



# bench_list = [ 'ghz', 'knn', 'multiplier', 'qft', ]
# bench_list = ['adder', 'bv', 'cc', 'dnn', 'ghz', 'knn', 'multiplier', 'qft', 'cat', 'ising', 'qugan', 'square_root', 'swap_test', 'wstate', ]#'bwt', #'vqe_uccsd']
bench_list = ['ghz','qft', 'cat', 'ising', 'wstate', ]#'bwt', #'vqe_uccsd']
# bench_list = ['bwt',] #'vqe_uccsd']
N_Bench_list = {'adder':[10,28,64,118,433], 'bv':[14,19,30,70,140,280], 'bwt':[21,37,57,97,177], 'cc':[12,32,64,151,301], 'dnn':[8,16,33,51], 'ghz':[23,40,78,127,255], 'knn':[25,31,41,67,129,341], 'multiplier':[15,45,75,350,400], 'qft':[4,18,29,63,160,320], 'cat':[4,22,35,65,130,260], 'ising':[10,26,34,42,66,98,420], 'qugan':[39, 71, 111, 395], 'square_root':[18,45,60], 'swap_test':[25,41,83,115,361], 'vqe_uccsd':[4,6,8,28], 'wstate':[3,27,36,76,118,380]}
# method_list = ['base', "move_split", "change_dest", "change_dest+move_split"]
# method_list = ['base', 'move_split', 'change_dest', 'change_dest+move_split', "break_chains", "break_chains+change_dest", "break_chains+change_dest+move_split"]
# method_list = ['base', 'move_split', 'change_dest', 'change_dest+move_split', "break_chains", "break_chains+change_dest", "break_chains+change_dest+move_split"]
method_list = ['move_split_double',  'change_dest+move_split_double', "break_chains+change_dest+move_split_double"]
# method_list = ["break_chains", "break_chains+change_dest", "break_chains+change_dest+move_split"]
for method in method_list:
    print(method)

    # for d in Distance_List:
    d=1
    # for n in [6, 15, 30, 40, 50, 70, 80]:
    # for n in [6, 10, 20, 30, 40, 50, 60, 80, 100, 150, 250]:
    # for n in [10, 26, 34, 42, 66, 98, 420]:#Ising
    # for n in [4, 22, 35, 65, 130, 260]:#cat
    for benchm in bench_list:
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
        print(benchm)
        for n in N_Bench_list[benchm]:
        # for n in [10, 20, 30, 50]:
            index = random.choice(I_List)
            Row = math.ceil(math.sqrt(n))


            try:
                # circ = QuantumCircuit.from_qasm_file(f"benchmarks\QASMBench-master\ising_n{n}\ising_n{n}_transpiled.qasm")
                circ = QuantumCircuit.from_qasm_file(f"benchmarks\\QASMBench-master\\{benchm}\\{benchm}_n{n}\\{benchm}_n{n}_transpiled.qasm")
                # circ_path = Path("benchmarks") / "QASMBench-master" / benchm / f"{benchm}_n{n}" / f"{benchm}_n{n}_transpiled.qasm"
                # circ = QuantumCircuit.from_qasm_file(circ_path)
            except FileNotFoundError:
                # circ = QuantumCircuit.from_qasm_file(f"benchmarks\QASMBench-master\ising_n{n}\ising_n{n}.qasm")
                # continue
                circ = QuantumCircuit.from_qasm_file(f"benchmarks\\QASMBench-master\\{benchm}\\{benchm}_n{n}\\{benchm}_n{n}.qasm")

            # circ = QuantumCircuit.from_qasm_file(f"benchmarks\QASMBench-master\medium\qec9xz_n17\qec9xz_n17.qasm")
            # circ = QuantumCircuit.from_qasm_file(f"benchmarks/vqe/vqe_n{n}.qasm")
            # with open(f"benchmarks/qft/qftn{n}.qasm", "r") as f:
            # with open(f"benchmarks/bv/new_bv_n{n}.qasm", "r") as f:
            # with open(f"benchmarks/rca/rca_n{n}.qasm", "r") as f:
            # with open(f"benchmarks/vqe/vqe_n{n}.qasm", "r") as f:
                # qasm_code = f.read()
            # circ = loads(qasm_code)

            test_circuit = transpile(circ, basis_gates=["u1", "u2", "u3", "cz", "id"],  optimization_level=2)
            cz_blocks = get_cz_blocks(test_circuit)
            # cz_blocks = get_cz_blocks(circ)

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
            try:
                no_storage_transfer_duration, no_storage_move_duration, no_storage_cir_fidelity, no_storage_cir_fidelity_1q_gate, no_storage_cir_fidelity_2q_gate, no_storage_cir_fidelity_2q_gate_for_idle, no_storage_cir_fidelity_atom_transfer, no_storage_cir_fidelity_coherence, no_storage_nstage, count, loop_num, split_succ, split_fail = mvqc(cz_blocks, Row, n, False, d, 1, method)
            except:
                print("split", benchm)
                break

            sorted(count.items())
            # print(loop_num)
            # print(count)
            print(split_succ, split_fail)
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

        with open(f"data/compare_{benchm}_{method}.txt", 'w') as file:
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

