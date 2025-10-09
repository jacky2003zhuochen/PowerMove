import sys
import os

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

from PowerMove import *
from Enola import *
import random
import math
import time

N_Qubit_List = [6, 10, 20, 30, 40, 50, 60, 80, 100]
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
# mvqc_runtime = []

# enola_transfer_duration_list = []
# enola_move_duration_list = [] 
# enola_cir_fidelity_list = [] 
# enola_cir_fidelity_1q_gate_list = [] 
# enola_cir_fidelity_2q_gate_list = [] 
# enola_cir_fidelity_2q_gate_for_idle_list = [] 
# enola_cir_fidelity_atom_transfer_list = [] 
# enola_cir_fidelity_coherence_list = []
# enola_nstage_list = []
# enola_runtime = []

index = random.choice(I_List)
type = 'regular'
d = 1
for degree in range(1,11):
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
    q_list = []
    for n in N_Qubit_List:
        Row = math.ceil(math.sqrt(n))

        path = f"benchmarks/qaoa/{type}/q{n}_regular{degree}/i{index}.txt"

        try:
            with open(path, "r") as fid:
                gates = eval(fid.read())
                q_list.append(n)
        except FileNotFoundError:
            continue
        
        mvqc_start_time = time.time()
        mvqc_transfer_duration, mvqc_move_duration, mvqc_cir_fidelity, mvqc_cir_fidelity_1q_gate, mvqc_cir_fidelity_2q_gate, mvqc_cir_fidelity_2q_gate_for_idle, mvqc_cir_fidelity_atom_transfer, mvqc_cir_fidelity_coherence, mvqc_nstage = mvqc([gates], Row, n, False, d, 1, method)
        mvqc_runtime.append(time.time() - mvqc_start_time)
        
        mvqc_transfer_duration_list.append(mvqc_transfer_duration)
        mvqc_move_duration_list.append(mvqc_move_duration)
        mvqc_cir_fidelity_list.append(mvqc_cir_fidelity)
        mvqc_cir_fidelity_1q_gate_list.append(mvqc_cir_fidelity_1q_gate)
        mvqc_cir_fidelity_2q_gate_list.append(mvqc_cir_fidelity_2q_gate)
        mvqc_cir_fidelity_2q_gate_for_idle_list.append(mvqc_cir_fidelity_2q_gate_for_idle)
        mvqc_cir_fidelity_atom_transfer_list.append(mvqc_cir_fidelity_atom_transfer)
        mvqc_cir_fidelity_coherence_list.append(mvqc_cir_fidelity_coherence)   
        mvqc_nstage_list.append(mvqc_nstage)

        # enola_start_time = time.time()
        # enola_transfer_duration, enola_move_duration, enola_cir_fidelity, enola_cir_fidelity_1q_gate, enola_cir_fidelity_2q_gate, enola_cir_fidelity_2q_gate_for_idle, enola_cir_fidelity_atom_transfer, enola_cir_fidelity_coherence, enola_nstage = enola([gates], Row, n, d)
        # enola_runtime.append(time.time() - enola_start_time)
        
        # enola_transfer_duration_list.append(enola_transfer_duration)
        # enola_move_duration_list.append(enola_move_duration)
        # enola_cir_fidelity_list.append(enola_cir_fidelity)
        # enola_cir_fidelity_1q_gate_list.append(enola_cir_fidelity_1q_gate)
        # enola_cir_fidelity_2q_gate_list.append(enola_cir_fidelity_2q_gate)
        # enola_cir_fidelity_2q_gate_for_idle_list.append(enola_cir_fidelity_2q_gate_for_idle)
        # enola_cir_fidelity_atom_transfer_list.append(enola_cir_fidelity_atom_transfer)
        # enola_cir_fidelity_coherence_list.append(enola_cir_fidelity_coherence)  
        # enola_nstage_list.append(enola_nstage)

    with open(f"data/qaoa_regular{degree}_no_storage_compare_change_dest+move_split.txt", 'w') as file:
        # file.write(str(N_Qubit_List) + '\n')
        file.write(str(q_list) + '\n')        
        file.write(str([x + y for x, y in zip(mvqc_transfer_duration_list, mvqc_move_duration_list)]) + '\n') 
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
