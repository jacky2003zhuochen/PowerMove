import math
Infinity = math.inf
# def storage_gate_scheduling(gates, storage_flag, L):
#     # print(gates)
#     colored_gates = []
#     for g in gates:
#         new_color = True
#         i = 0
#         # colored_gates.sort(key = len)
#         # print("color gates", colored_gates)
#         for cgs in colored_gates:
#             conflict_flag = False
#             for cg in cgs:
#                 if g[0] == cg[0]:
#                     conflict_flag = True
#                     break
#                 elif g[0] == cg[1]:
#                     conflict_flag = True
#                     break
#                 elif g[1] == cg[0]:
#                     conflict_flag = True
#                     break
#                 elif g[1] == cg[1]:
#                     conflict_flag = True
#                     break
#             if not conflict_flag:
#                 if len(cgs) < L:
#                     new_color = False
#                     colored_gates[i].append(g)
#                     break
#             i += 1
#         if new_color:
#             colored_gates.append([g])
#     if storage_flag:
#         colored_gates.sort(key = len)
#         lg = []
#         color_num = len(colored_gates)
#         for i in range(color_num):
#             if i == 0:
#                 mg = colored_gates[0]
#                 lg.append(mg)
#                 colored_gates.remove(mg)
#             else:
#                 min_diff = Infinity
#                 pre_mg = lg[-1]
#                 target_mg = []
#                 for mg in colored_gates:
#                     pre_interaction_qubits = []
#                     for m in pre_mg:
#                         pre_interaction_qubits.append(m[0])
#                         pre_interaction_qubits.append(m[1]) 
#                     interaction_qubits = []
#                     for m in mg:
#                         interaction_qubits.append(m[0])
#                         interaction_qubits.append(m[1])     
#                     cur_diff = 0
#                     for q in interaction_qubits:
#                         if q not in pre_interaction_qubits:
#                             cur_diff += 1

#                     # for q in pre_interaction_qubits:
#                     #     if q not in interaction_qubits:
#                     #         cur_diff += 0.1   
#                     if cur_diff < min_diff:
#                         min_diff = cur_diff
#                         target_mg = mg
#                     # print("stage", i, "cur_diff", cur_diff, pre_interaction_qubits, interaction_qubits)
#                 lg.append(target_mg)
#                 colored_gates.remove(target_mg) 
#         return lg
#     else:
            
#         return colored_gates

import networkx as nx

def storage_gate_scheduling(gates, storage_flag, L):

    """
    基于容量限制图染色的 CZ Gate Stage Scheduling 算法
    
    参数:
        cz_list: list of tuples, 例如 [(0, 1), (1, 2), (0, 2), ...]
        L: int, 每个 stage 允许的最大 CZ 数量 (容量限制)
    返回:
        stages: list of lists, 包含所有划分好 stage 的 cz gates
    """
    if not gates:
        return []

    n_gates = len(gates)
    G = nx.Graph()
    
    # 1. 构建冲突图 (Conflict Graph)
    # 我们使用 index 作为节点，以防列表中存在完全相同的 cz gate 导致节点覆盖
    G.add_nodes_from(range(n_gates))
    
    for i in range(n_gates):
        for j in range(i + 1, n_gates):
            # 如果两个 CZ 门共享了至少一个 qubit，则连边表示冲突
            if set(gates[i]) & set(gates[j]):
                G.add_edge(i, j)
                
    # 2. 节点排序：按度数（冲突数量）降序排列 (Welsh-Powell 启发式)
    # 优先为冲突最多的门分配 Stage，更容易得到较少的总 Stage 数
    nodes_sorted = sorted(G.nodes, key=lambda x: G.degree[x], reverse=True)
    
    # stages_indices 用于记录每个 stage 包含的 gate 索引
    stages_indices = []
    
    # 3. 贪心染色
    for node in nodes_sorted:
        placed = False
        
        # 尝试放入已有的 stage
        for stage in stages_indices:
            # 约束 1: 检查容量限制 L
            if len(stage) >= L:
                continue
                
            # 约束 2: 检查冲突 (当前 stage 中是否有任何 gate 与 node 有边相连)
            has_conflict = any(G.has_edge(node, existing_node) for existing_node in stage)
            
            if not has_conflict:
                stage.append(node)
                placed = True
                break
                
        # 如果所有现存的 stage 都放不下，开辟一个新的 stage
        if not placed:
            stages_indices.append([node])
            
    # 4. 将 index 映射回实际的 gate tuple
    stages = [[gates[idx] for idx in stage] for stage in stages_indices]
    
    return stages

def sort_stages(colored_gates):

    colored_gates.sort(key = len)
    lg = []
    color_num = len(colored_gates)
    # print('color num', color_num)
    for i in range(color_num):
        if i == 0:
            mg = colored_gates[0]
            lg.append(mg)
            colored_gates.remove(mg)
        else:
            min_diff = -1
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
                    if q in pre_interaction_qubits:
                        cur_diff += 1
                # print(interaction_qubits)
                # print(pre_interaction_qubits)
                # print(cur_diff)
                # for q in pre_interaction_qubits:
                #     if q not in interaction_qubits:
                #         cur_diff += 0.1   
                if cur_diff > min_diff:
                    min_diff = cur_diff
                    target_mg = mg
                # print("stage", i, "cur_diff", cur_diff, pre_interaction_qubits, interaction_qubits)
            lg.append(target_mg)
            colored_gates.remove(target_mg) 
    
    # print('lg', len(lg))
    return lg

# ==========================================
# 测试用例
# ==========================================
# if __name__ == "__main__":
#     # 假设的 CZ 门列表
#     cz_gates_input = [(0, 1), (1, 2), (2, 3), (0, 3), (0, 2), (1, 3), (4, 5)]
    
#     # 设置每个 stage 最多允许并行执行的 CZ 门数量
#     L_limit = 2
    
#     scheduled_stages = storage_gate_scheduling(cz_gates_input, L_limit)
    
#     print(f"最大并行度 L = {L_limit}")
#     print("-" * 30)
#     for i, stage in enumerate(scheduled_stages):
#         print(f"Stage {i}: {stage}")