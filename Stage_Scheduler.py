import math
Infinity = math.inf
def storage_gate_scheduling(gates, storage_flag):
    # print(gates)
    colored_gates = []
    for g in gates:
        new_color = True
        i = 0
        # colored_gates.sort(key = len)
        # print("color gates", colored_gates)
        for cgs in colored_gates:
            conflict_flag = False
            for cg in cgs:
                if g[0] == cg[0]:
                    conflict_flag = True
                    break
                elif g[0] == cg[1]:
                    conflict_flag = True
                    break
                elif g[1] == cg[0]:
                    conflict_flag = True
                    break
                elif g[1] == cg[1]:
                    conflict_flag = True
                    break
            if not conflict_flag:
                new_color = False
                colored_gates[i].append(g)
                break
            i += 1
        if new_color:
            colored_gates.append([g])
    if storage_flag:
        colored_gates.sort(key = len)
        lg = []
        color_num = len(colored_gates)
        for i in range(color_num):
            if i == 0:
                mg = colored_gates[0]
                lg.append(mg)
                colored_gates.remove(mg)
            else:
                min_diff = Infinity
                pre_mg = lg[-1]
                target_mg = []
                for mg in colored_gates:
                    pre_interaction_qubits = []
                    for m in pre_mg:
                        pre_interaction_qubits.append(m[0])
                        pre_interaction_qubits.append(m[1]) 
                    interaction_qubits = []
                    for m in mg:
                        interaction_qubits.append(m[0])
                        interaction_qubits.append(m[1])     
                    cur_diff = 0
                    for q in interaction_qubits:
                        if q not in pre_interaction_qubits:
                            cur_diff += 1

                    # for q in pre_interaction_qubits:
                    #     if q not in interaction_qubits:
                    #         cur_diff += 0.1   
                    if cur_diff < min_diff:
                        min_diff = cur_diff
                        target_mg = mg
                    print("stage", i, "cur_diff", cur_diff, pre_interaction_qubits, interaction_qubits)
                lg.append(target_mg)
                colored_gates.remove(target_mg) 
        return lg
    else:
            
        return colored_gates