import sys
import os

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

from PowerMove import *
import random
import math

N_Qubit_List = [5, 10, 20, 30 ,50]
NAOD_list = [1, 2, 3, 4]
P_List = [0.1, 0.2, 0.3, 0.4, 0.5]
I_List = range(10)


p = 0.5
index = random.choice(I_List)
type = 'rand'

with open("data/qaoa_random_storage_compare_multi_aod.txt", 'w') as file:
    file.write(str(NAOD_list) + '\n')
    for NAOD in NAOD_list:
        mvqc_transfer_duration_list = []
        mvqc_move_duration_list = [] 
        mvqc_cir_fidelity_list = [] 
        mvqc_cir_fidelity_1q_gate_list = [] 
        mvqc_cir_fidelity_2q_gate_list = [] 
        mvqc_cir_fidelity_2q_gate_for_idle_list = [] 
        mvqc_cir_fidelity_atom_transfer_list = [] 
        mvqc_cir_fidelity_coherence_list = []
        mvqc_nstage_list = []
        for n in N_Qubit_List:
            Row = math.ceil(math.sqrt(n))

            path = f"benchmarks/qaoa/{type}/q{n}_p{p}/i{index}.txt"


            with open(path, "r") as fid:
                gates = eval(fid.read())

            mvqc_transfer_duration, mvqc_move_duration, mvqc_cir_fidelity, mvqc_cir_fidelity_1q_gate, mvqc_cir_fidelity_2q_gate, mvqc_cir_fidelity_2q_gate_for_idle, mvqc_cir_fidelity_atom_transfer, mvqc_cir_fidelity_coherence, mvqc_nstage = mvqc([gates], Row, n, True, 1, NAOD)
            mvqc_transfer_duration_list.append(mvqc_transfer_duration)
            mvqc_move_duration_list.append(mvqc_move_duration)
            mvqc_cir_fidelity_list.append(mvqc_cir_fidelity)
            mvqc_cir_fidelity_1q_gate_list.append(mvqc_cir_fidelity_1q_gate)
            mvqc_cir_fidelity_2q_gate_list.append(mvqc_cir_fidelity_2q_gate)
            mvqc_cir_fidelity_2q_gate_for_idle_list.append(mvqc_cir_fidelity_2q_gate_for_idle)
            mvqc_cir_fidelity_atom_transfer_list.append(mvqc_cir_fidelity_atom_transfer)
            mvqc_cir_fidelity_coherence_list.append(mvqc_cir_fidelity_coherence)   
            mvqc_nstage_list.append(mvqc_nstage)

        file.write(str(N_Qubit_List) + '\n')
        file.write(str(mvqc_transfer_duration_list) + '\n') 
        file.write(str(mvqc_move_duration_list) + '\n') 
        file.write(str(mvqc_cir_fidelity_list) + '\n') 
        file.write(str(mvqc_cir_fidelity_1q_gate_list) + '\n') 
        file.write(str(mvqc_cir_fidelity_2q_gate_list) + '\n') 
        file.write(str(mvqc_cir_fidelity_2q_gate_for_idle_list) + '\n') 
        file.write(str(mvqc_cir_fidelity_atom_transfer_list) + '\n') 
        file.write(str(mvqc_cir_fidelity_coherence_list) + '\n')
        file.write(str(mvqc_nstage_list) + '\n')
