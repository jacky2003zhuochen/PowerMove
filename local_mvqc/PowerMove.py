from scheduler.gate_scheduler import gate_scheduling
from placer.placer import place_qubit
from Coll_Moves_Scheduler import *
from Continuous_Router import *
from Stage_Scheduler import *
import matplotlib.pyplot as plt
import networkx as nx
import math
from collections import deque

N_Block = 1
X_SEP = 10
Y_SEP = 10
SITE_SEP = 2
STORAGE_SEP = 3
Infinity = math.inf
Storage_Y_SEP = 10
Fidelity_2Q_Gate = 0.995
Fidelity_1Q_Gate = 0.995
Fidelity_Atom_Transfer = 0.999
Coherence_Time = 1.5e6 # ms
MUS_PER_FRM = 8

def update(pos, qubit_map, moves, storage_in_move):
    for m in moves:
        if qubit_map.has_edge(m[0], m[1]):
            qubit_map.remove_edge(m[0], m[1])
        # empty_space[pos[m[0]]].remove(m[0])
        # empty_space[pos[m[1]]].append(m[0])
        qubit_map.nodes[m[0]]['pos'] = qubit_map.nodes[m[1]]['pos']
        pos[m[0]] = pos[m[1]]

    for sq in storage_in_move.keys():
        # empty_space[pos[sq]].remove(sq)
        # empty_space[storage_in_move[sq]].append(sq)
        qubit_map.nodes[sq]['pos'] = storage_in_move[sq]
        pos[sq] = storage_in_move[sq]

    return pos, qubit_map

def sabre_initial_layout(cz_blocks, N, dimension):
    from qiskit import QuantumCircuit, transpile  # type: ignore
    from qiskit.transpiler import CouplingMap  # type: ignore

    X = dimension[0]
    Y = dimension[1]

    # build coupling graph
    pos_dict = {}
    coupling_graph = []
    for i in range(X):
        for j in range(Y):
            pos_dict[i * Y + j] = (i, j)
            if i:
                coupling_graph.append([i * Y + j, (i - 1) * Y + j])
                coupling_graph.append([(i - 1) * Y + j, i * Y + j])
            if j:
                coupling_graph.append([i * Y + j, i * Y + j - 1])
                coupling_graph.append([i * Y + j - 1, i * Y + j])

    coupling_graph = CouplingMap(couplinglist=coupling_graph)

    # build a virtual circuit
    circ = QuantumCircuit(N)
    for block in cz_blocks:
        for cz_gate in block:
            circ.cz(cz_gate[0], cz_gate[1])

    tcirc = transpile(
        circ,
        coupling_map=coupling_graph,
        layout_method="sabre",
        routing_method="sabre",
        seed_transpiler=42,
    )
    qiskit_layout = tcirc.layout

    if qiskit_layout is None:
        raise RuntimeError("Sabre algorithm failed to find a layout")

    layout = qiskit_layout.initial_layout

    qubit_mapping = {}
    idx = 0
    for lbit in layout.get_virtual_bits():
        if idx >= N:
            break
        qubit_mapping[idx] = pos_dict[layout[lbit]]
        idx += 1
    return qubit_mapping

def update_counter_auto(counter: dict, arr: np.ndarray):
    unique, counts = np.unique(arr, return_counts=True)
    for u, c in zip(unique, counts):
        if u not in counter:
            counter[u] = 0
        counter[u] += c
    return counter

def mvqc(cz_blocks, Row, n, storage_flag, d, num_aod, method, cost_para, para1, cost_para2, para2, thre, h, a, location_size=2, iter_num=2):
    N_Block = d
    # Row *= 2
    # print(Row**2)
    Row = Row*a
    Row1 = Row
    Row2 = 4*Row
    Column1 = math.ceil(Row/(a**2))
    Column1 -= h
    Column2 = math.ceil(Row/(2*a**2))
    # Column2 = math.ceil(n/Row2)

    if not storage_flag:
        list_gates = []
        for gates in cz_blocks:
            list_gates += storage_gate_scheduling(gates, storage_flag, Row*Row)

    else:
        L = Column1*Row1
        # print("L", L)
        list_gates = []
        for gates in cz_blocks:
            list_gates += storage_gate_scheduling(gates, storage_flag, L)

    # list_gates.sort(key=len, reverse=True)  # 可选：按每个 stage 中 gate 的数量排序，便于后续处理

    if storage_flag and h != 0:
        list_gates = sort_stages(list_gates)

    # for gates in list_gates:
        # print(len(gates), "num of cz")
    # print("cz")
    qubit_mapping = list(sabre_initial_layout(cz_blocks, n, (Row2, Column2)).values())
    # qubit_mapping = place_qubit((Row, Row), n, list_gates, True)

    cir_fidelity_2q_gate = 1
    cir_fidelity_2q_gate_for_idle = 1 
    cir_fidelity_atom_transfer = 1 
    cir_fidelity_1q_gate = 1
    cir_fidelity_coherence = 1
    fidelity_2q_gate_for_idle = 1 - (1-Fidelity_2Q_Gate)/2

    num_movement_stage = 0
    
    cir_qubit_idle_time = [0] * n
    list_movement_duration = []
    list_transfer_duration = []
    count = defaultdict(int)
    loop_num = 0

    empty_space = {}
    for i in range(Row2):
        for j in range(Column2):
            empty_space[(i, -1-j)] = []
    for i in range(Row1):
        for j in range(Column1):
            empty_space[(i, j)] = []

    # draw cz_graph layout
    qubit_map = nx.Graph()
    # max_y = 0
    # for qk in qubit_mapping:
    #     if qk[1] > max_y:
    #         max_y = qk[1]

    for i in range(len(qubit_mapping)):
        if not storage_flag:
            qubit_pos = qubit_mapping[i]
            qubit_map.add_node(i)
            qubit_map.nodes[i]['pos'] = qubit_pos
            empty_space[qubit_pos].append(i)
        else:
            qubit_pos = qubit_mapping[i]
            qubit_map.add_node(i)
            qubit_map.nodes[i]['pos'] = (qubit_pos[0], -1 - qubit_pos[1]) # storage位置的y坐标取负数，方便区分
            empty_space[qubit_map.nodes[i]['pos']].append(i)    
    # print("qubit mapping", qubit_mapping)
    initial_space = {k: v[:] for k, v in empty_space.items()}
    for gates in cz_blocks:
        for gate in gates:
            qubit_map.add_edge(gate[0], gate[1])

    pos = nx.get_node_attributes(qubit_map, 'pos')
    location_index = {}
    target_location_index = {}
    for q in range(n):
        location_index[q] = 0
        target_location_index[q] = 0

    # nx.draw(qubit_map, pos = pos, node_size = 50)
    # nx.draw_networkx_labels(qubit_map, pos, labels={i: str(i) for i in qubit_map.nodes()}, font_color='black')

    s_index = 0
    if storage_flag:
        qubits_not_in_storage = set()
    else:
        qubits_not_in_storage = set(range(n))

    if "SLM" in method:
        location_size = int(method[-1])
        # location_size = int(method[method.find("slm")+1])
    # storage_occ = {}

    # for i in range(Row):
    #     storage_occ[i] = -2
 
    # print(qubit_map.nodes(), pos)
    #[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] 
    #{0: (2, 3), 1: (3, 1), 2: (0, 2), 3: (2, 1), 4: (1, 2), 5: (3, 2), 6: (2, 0), 7: (1, 3), 8: (3, 4), 9: (4, 1), 10: (0, 3), 11: (2, 2), 12: (4, 0), 13: (2, 4), 14: (3, 3), 15: (1, 1), 16: (4, 3), 17: (3, 0), 18: (0, 1), 19: (1, 4)}
    list_gates = deque(list_gates)
    # print(list_gates, "list_gates")
    split_succ, split_fail = 0, 0
    pick_drop_orig_num, pick_drop_split_num = 0, 0
    aod_sum, AOD_moveTime_sum, AOD_move_num_sum, max_length = 0, 0, 0, 0
    AOD_moveTime_span = (10000000,0)
    counter = {i: 0 for i in range(1, 3)} 
    Max_AOD_size = 0
    while len(list_gates)>0:
        mg = list_gates.popleft()
        # print(mg)

        move_in_qubits = []
        move_out_qubits = []

        if storage_flag:
            interaction_qubits = set()
            for m in mg:
                interaction_qubits.add(m[0])
                interaction_qubits.add(m[1])
            # print(mg)
            for q in interaction_qubits:
                if q not in qubits_not_in_storage:
                    move_out_qubits.append(q)
                    qubits_not_in_storage.add(q)

            move_in_qubits = [q for q in qubits_not_in_storage if q not in interaction_qubits]
            qubits_not_in_storage -= set(move_in_qubits)

        

        location_index = target_location_index.copy()
        initial_space = {k: v[:] for k, v in empty_space.items()}

        move_group, empty_space, rq_moved_pos, qmg, storage_in_move, target_location_index = continuous_router(mg, pos, empty_space, move_out_qubits, move_in_qubits, storage_flag, Row, location_index, target_location_index, location_size, method, 0, a)

        move_distance = {}
        # for move in move_group:
        #     if move[1][1] >= 0 and pos[move[0]][1] >= 0:
        #         move_distance[move] = abs(move[1][0] - pos[move[0]][0]) * X_SEP * d + abs(move[1][1] - pos[move[0]][1]) * Y_SEP * d 
        #     elif move[1][1] >= 0:
        #         move_distance[move] = abs(move[1][0] - pos[move[0]][0]) * X_SEP * d  + abs(move[1][1] + 2) * Y_SEP * d  + abs(-2 - pos[move[0]][1]) * Storage_Y_SEP * d 
        #     elif pos[move[0]][1] >= 0:
        #         move_distance[move] = abs(move[1][0] - pos[move[0]][0]) * X_SEP * d  + abs(pos[move[0]][1] + 2) * Y_SEP * d  + abs(-2 - move[1][1] ) * Storage_Y_SEP * d 
        #     else:
        #         move_distance[move] = abs(move[1][0] - pos[move[0]][0]) * X_SEP * d  + abs(move[1][1] - pos[move[0]][1]) * Storage_Y_SEP * d 

        # revised move distance calculation with SiTE_SEP
        for move in move_group:
            source = move[1]
            dest = move[2]
            if dest[1] >= 0 and source[1] >= 0:
                move_distance[move] = abs((dest[0] - source[0]) * X_SEP + (dest[2] - source[2]) * SITE_SEP) * d + abs(dest[1] - source[1]) * Y_SEP * d 
            elif dest[1] >= 0:
                move_distance[move] = abs((dest[0]*4 - source[0]) * X_SEP/4 + (dest[2] - source[2]) * SITE_SEP) * d + abs(dest[1]) * Y_SEP * d + Storage_Y_SEP * d + abs(source[1]+1) * SOTORAGE_SEP * d 
            elif source[1] >= 0:
                move_distance[move] = abs((dest[0] - source[0]*4) * X_SEP/4 + (dest[2] - source[2]) * SITE_SEP) * d + abs(source[1]) * Y_SEP * d + Storage_Y_SEP * d + abs(dest[1]+1) * SOTORAGE_SEP * d 
            else:
                move_distance[move] = abs((dest[0] - source[0]) * X_SEP/4 + (dest[2] - source[2]) * SITE_SEP) * d + abs(dest[1] - source[1]) * SOTORAGE_SEP * d 

        # for m in qmg:
        #     if empty_space[pos[m[1]]].find(m[0]) == -1:
        #         empty_space[pos[m[1]]].append(m[0])
        # coll moves
        empty_space, parallel_move_groups, num_movement_stage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration, list_movement_duration, target_location_index, change_dest, move_in_loop, count, loop_num, aod_num, AOD_moveTime, AOD_move_num, span, m_length = coll_moves_scheduler(empty_space, initial_space, n, Row, move_distance, move_group, num_aod, move_in_qubits, move_out_qubits,
        qubits_not_in_storage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration,
        list_movement_duration, num_movement_stage, location_index, target_location_index, location_size, 
        method, count, loop_num, split_fail, split_succ,cost_para, para1, cost_para2, para2, thre, h, a, iter_num=iter_num)

        # print(aod_num, "aod_num")
        aod_sum += aod_num
        AOD_moveTime_sum += AOD_moveTime
        AOD_move_num_sum += np.sum(np.array(AOD_move_num))
        update_counter_auto(counter, AOD_move_num)
        max_length = max(max_length, m_length)
        AOD_moveTime_span = (min(AOD_moveTime_span[0],span[0]), max(AOD_moveTime_span[1],span[1]))
        # pick_drop_orig_num += pick_drop_orig
        # pick_drop_split_num += pick_drop_split

        ###################################################################################################################
        # trivial task reduction
        new_mg = []

        if len(list(set(move_in_loop))) != 0:
            # print("move in loop", move_in_loop)
            for m in move_in_loop:
                q1 = m[0]
                q_list = empty_space[(m[2][0], m[2][1])][:]
                q_list.remove(q1)
                if len(q_list)>0:
                    q2 = q_list[0]
                    new_mg.append((q1,q2))
                    if (q1,q2) in mg:
                        mg.remove((q1,q2))
                    else:
                        mg.remove((q2,q1))
                else:
                    q_list = empty_space[(m[1][0], m[1][1])][:]
                    new_mg.append((q_list[0],q_list[1]))
                    if (q_list[0],q_list[1]) in mg:
                        mg.remove((q_list[0],q_list[1]))
                    else:
                        # if q_list == [76, 38]:
                        #     print(initial_space[(m[1][0], m[1][1])])
                        #     print(mg)
                        mg.remove((q_list[1],q_list[0]))
            empty_space = {k: v[:] for k, v in initial_space.items()}
            target_location_index = location_index.copy()
            list_gates.appendleft(new_mg)
            list_gates.appendleft(mg)
            if storage_flag:
                qubits_not_in_storage.update(move_in_qubits)
                qubits_not_in_storage -= set(move_out_qubits)
            # print("task split")
            continue

        cir_fidelity_2q_gate *= pow(Fidelity_2Q_Gate, len(mg))
        if not storage_flag:
            cir_fidelity_2q_gate_for_idle *= pow(fidelity_2q_gate_for_idle, (n - 2*len(mg)) * (N_Block ** 2))
        else:
            cir_fidelity_2q_gate_for_idle *= pow(fidelity_2q_gate_for_idle, (len(qubits_not_in_storage) - 2*len(mg)) * (N_Block ** 2))
            if (len(qubits_not_in_storage) - 2*len(mg)) != 0:
                print("error: idle qubits in storage also suffer fidelity decay")

        # print("one stage finished")

        # parallel_move_groups.sort(reverse = True, key = len)
        # print("move steps", len(parallel_move_groups))
        # for selective_move_group in parallel_move_groups:

        for rq in rq_moved_pos.keys():
            # empty_space[pos[rq]].remove(rq)
            pos[rq] = rq_moved_pos[rq]
            # empty_space[pos[rq]].append(rq)
            qubit_map.nodes[rq]['pos'] = pos[rq]

        # print(change_dest)

        pos, qubit_map = update(pos, qubit_map, qmg, storage_in_move)
        for cq in change_dest.keys():
            pos[cq] = change_dest[cq]
            qubit_map.nodes[cq]['pos'] = pos[cq]
        # print("empty space")
        # for p in empty_space.keys():
        #     if len(empty_space[p]):
        #         print(p, empty_space[p])
        
        s_index += 1

        # index = 0
        # for p in empty_space.keys():
        #     index += len(empty_space[p])
        # print("index", index)
        # pos_redundant_graph = nx.Graph()
        # for p in empty_space.keys():
        #     pos_qubits = list(empty_space[p])
        #     for q in pos_qubits:
        #         for m in qmg:
        #             if q == m[0]:
        #                 pos_qubits.remove(q)
        #                 pos_qubits.remove(m[1])
        #                 break
        #             elif q == m[1]:
        #                 pos_qubits.remove(q)
        #                 pos_qubits.remove(m[0])
        #                 break
        #     if len(pos_qubits):
        #         print(pos_qubits)
        #         pos_redundant_graph.add_node(str(pos_qubits))
        #         pos_redundant_graph.nodes[str(pos_qubits)]['pos'] = p

        # print(qubit_map.nodes(), pos)
        # s_index += 1

    # print("aod_num", aod_sum)
    for t in cir_qubit_idle_time:
        cir_fidelity_coherence = cir_fidelity_coherence * ((1 - t/Coherence_Time) ** (N_Block ** 2))
    cir_fidelity = cir_fidelity_1q_gate * cir_fidelity_2q_gate * cir_fidelity_2q_gate_for_idle \
                        * cir_fidelity_atom_transfer * cir_fidelity_coherence
    # print("cir_qubit_idle_time", cir_qubit_idle_time)
    # print(np.mean(np.array(cir_qubit_idle_time)), np.sum(np.array(cir_qubit_idle_time)))
    # print("cir_fidelity_1q_gate", cir_fidelity_1q_gate)
    # print("cir_fidelity_2q_gate", cir_fidelity_2q_gate)
    # print("cir_fidelity_2q_gate_for_idle", cir_fidelity_2q_gate_for_idle)
    # print("cir_fidelity_atom_transfer", cir_fidelity_atom_transfer)
    # print("cir_fidelity_coherence", cir_fidelity_coherence)
    # print("coherence_time", Coherence_Time)
    # print(s_index, "stages")
    return sum(list_transfer_duration), sum(list_movement_duration), cir_fidelity, cir_fidelity_1q_gate, cir_fidelity_2q_gate, cir_fidelity_2q_gate_for_idle, cir_fidelity_atom_transfer, cir_fidelity_coherence, num_movement_stage, count, loop_num, aod_sum, AOD_moveTime_sum/aod_sum, AOD_move_num_sum/aod_sum, AOD_moveTime_span, max_length, counter