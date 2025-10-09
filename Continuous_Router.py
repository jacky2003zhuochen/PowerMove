import copy 

def continuous_router(mg, pos, empty_space, move_out_qubits, move_in_qubits, storage_flag, Row, location_index, target_location_index, location_size, method):

    def get_pos(q):
        return pos[q][1]
    
    qmg = []
    static_pos = {}
    # static_num = {}
    moved_qubits = []
    redundant_qubits_to_be_moved = []
    init_space = copy.deepcopy(empty_space)
    for move in mg:
        q0 = move[0]
        q1 = move[1]
        pos_q0 = pos[q0]
        pos_q1 = pos[q1]
        
        if q0 in move_out_qubits and q1 in move_out_qubits:
            # print("count redundant0")
            qmg.append(move)
            moved_qubits.append(q0)
            redundant_qubits_to_be_moved.append(q1)
            empty_space[pos[q0]].remove(q0)
            empty_space[pos[q1]].remove(q1)
            
            # print('in move redundant', redundant_qubits_to_be_moved)
            
            # storage_occ[pos[q0][0]] += 1
            # storage_occ[pos[q1][0]] += 1
        elif q0 in move_out_qubits:
            # print("count redundant1")
            if pos_q1 not in static_pos.keys():
                # static_num[pos_q1] = static_num.get(pos_q1, 0) + 1
                # if static_num[pos_q1] >= location_size - 1:
                #     static_pos[pos_q1] = q1
                static_pos[pos_q1] = q1
                qmg.append(move)
                moved_qubits.append(q0)
                empty_space[pos[q0]].remove(q0)
            else:
                qmg.append(move)
                moved_qubits.append(q0)
                redundant_qubits_to_be_moved.append(q1)
                empty_space[pos[q0]].remove(q0)
                empty_space[pos[q1]].remove(q1)
                # print('in move redundant', redundant_qubits_to_be_moved)
            
            # storage_occ[pos[q0][0]] += 1                               
        elif q1 in move_out_qubits:
            # print("count redundant2")
            if pos_q0 not in static_pos.keys():
                # static_num[pos_q0] = static_num.get(pos_q0, 0) + 1
                # if static_num[pos_q0] >= location_size - 1:
                #     static_pos[pos_q0] = q0
                static_pos[pos_q1] = q0
                qmg.append((q1, q0))
                moved_qubits.append(q1)
                empty_space[pos[q1]].remove(q1)
            else:
                qmg.append((q1, q0))
                moved_qubits.append(q1)
                redundant_qubits_to_be_moved.append(q0)
                empty_space[pos[q0]].remove(q0)
                empty_space[pos[q1]].remove(q1)
                # print('in move redundant', redundant_qubits_to_be_moved)
            
            # storage_occ[pos[q1][0]] += 1                  
        else:
            # print("count redundant3")
            if pos_q1 not in static_pos.keys():
                # static_num[pos_q1] = static_num.get(pos_q1, 0) + 1
                # if static_num[pos_q1] >= location_size - 1:
                #     static_pos[pos_q1] = q1
                static_pos[pos_q1] = q1    
                qmg.append(move)
                moved_qubits.append(q0)
                empty_space[pos[q0]].remove(q0)
            elif pos_q0 not in static_pos.keys():
                # static_num[pos_q0] = static_num.get(pos_q0, 0) + 1
                # if static_num[pos_q0] >= location_size - 1:
                #     static_pos[pos_q0] = q0
                static_pos[pos_q0] = q0    
                qmg.append((q1, q0))
                moved_qubits.append(q1)
                empty_space[pos[q1]].remove(q1)

            else:
                qmg.append(move)
                moved_qubits.append(q0)
                redundant_qubits_to_be_moved.append(q1)
                # print(q0,q1)
                empty_space[pos[q0]].remove(q0)
                empty_space[pos[q1]].remove(q1)
                # print('in move redundant', redundant_qubits_to_be_moved)

    storage_in_move = {}
    # move in directly with the minimum distance
    if storage_flag:  
        move_in_qubits.sort(reverse = True, key = get_pos)
        for q in move_in_qubits:
            empty_space[pos[q]].remove(q)
            for y in range(-2, -8 * Row - 3, -1):
                if len(empty_space[(pos[q][0], y)]) == 0:
                    storage_in_move[q] = (pos[q][0], y)

    # print(redundant_qubits_to_be_moved)
    for p in empty_space.keys():
        placed_qubits = empty_space[p]
        for pq in placed_qubits:
            if pq not in move_in_qubits and pos[pq][1] >= 0:
                if p not in static_pos.keys():
                    static_pos[p] = pq
                if pq != static_pos[p] and pq not in moved_qubits and pq not in redundant_qubits_to_be_moved:
                    empty_space[pos[pq]].remove(pq)
                    redundant_qubits_to_be_moved.append(pq)
                    # print('not in move redundant', redundant_qubits_to_be_moved) # 将static位置上非static的qubit也加入到需要移动的队列中

    # print("redundant qubits to be moved", redundant_qubits_to_be_moved)

    rq_moved_pos = {}
    for rq in redundant_qubits_to_be_moved:
        pos_rq = pos[rq]
        pos_x = pos_rq[0]
        pos_y = pos_rq[1]
        pos_find_flag = False
        for r in range(20 * Row):
            for i in range(min(r + 1, Row)):
                j = r - i
                for a in [-1, 1]:
                    for b in [-1, 1]:
                        npos_x = pos_x + a * i
                        npos_y = pos_y + b * j
                        # print("npos_x,", npos_x, "npos_y,", npos_y)
                        # if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < Row:
                            # print("npos_x,", npos_x, "npos_y,", npos_y, "empty space", empty_space[(npos_x, npos_y)])
                        if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < Row and len(empty_space[(npos_x, npos_y)]) == 0:
                            # pos[rq] = (npos_x, npos_y)
                            empty_space[(npos_x, npos_y)].append(rq)
                            # qubit_map.nodes[rq]['pos'] = pos[rq]
                            rq_moved_pos[rq] = (npos_x, npos_y)
                            pos_find_flag = True
                            # if False:
                            if ((method == 'change_dest') or (method == 'change_dest+move_split')) and (len(init_space[(npos_x, npos_y)])==1):
                                target_location_index[rq] = (location_index[init_space[(npos_x,npos_y)][0]]+1)%2
                            else:
                                target_location_index[rq] = 0     
                            break
                    if pos_find_flag:
                        break

                if pos_find_flag:
                    break
            if pos_find_flag:
                break
        # print(pos_find_flag)

            # storage_in_move[q] = (pos[q][0], storage_occ[pos[q][0]])
            # storage_occ[pos[q][0]] -= 1

    move_group = []
    # formulate qmg (q0, q1) into move_group (q, pos)
    for qm in qmg:
        if pos[qm[1]] in static_pos.keys() and qm[1] == static_pos[pos[qm[1]]]:
            target_location_index[qm[0]] = (location_index[qm[1]] + 1) % location_size
            move_group.append((
                qm[0], 
                (pos[qm[0]][0], pos[qm[0]][1], location_index[qm[0]]), 
                (pos[qm[1]][0], pos[qm[1]][1], target_location_index[qm[0]])
            ))
            empty_space[pos[qm[1]]].append(qm[0])
        else:
            target_location_index[qm[0]] = (location_index[qm[1]] + 1) % location_size
            move_group.append((
                qm[0], 
                (pos[qm[0]][0], pos[qm[0]][1], location_index[qm[0]]), 
                (rq_moved_pos[qm[1]][0], rq_moved_pos[qm[1]][1], target_location_index[qm[0]])
            ))
            empty_space[rq_moved_pos[qm[1]]].append(qm[0])
    
    for rq in redundant_qubits_to_be_moved:
        move_group.append((
            rq, 
            (pos[rq][0], pos[rq][1], location_index[rq]), 
            (rq_moved_pos[rq][0], rq_moved_pos[rq][1], target_location_index[rq])
        ))

    for sq in storage_in_move.keys():
        target_location_index[sq] = 0
        move_group.append((
            sq, 
            (pos[sq][0], pos[sq][1], location_index[sq]), 
            (storage_in_move[sq][0], storage_in_move[sq][1], target_location_index[sq])
        ))
        empty_space[storage_in_move[sq]].append(sq)
    for m in move_group:
        if m[0] in empty_space[(m[1][0], m[1][1])]:
            empty_space[(m[1][0], m[1][1])].remove(m[0])
        if m[0] not in empty_space[(m[2][0], m[2][1])]:
            empty_space[(m[2][0], m[2][1])].append(m[0])
    return move_group, empty_space, rq_moved_pos, qmg, storage_in_move, target_location_index