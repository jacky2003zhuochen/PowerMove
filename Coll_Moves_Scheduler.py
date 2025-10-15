import networkx as nx
import copy
from collections import defaultdict
import numpy as np

def weighted_cosine(v, L, alpha=0.5):
    sims = []
    for w in L:
        cos = np.dot(v, w) / (np.linalg.norm(v) * np.linalg.norm(w))
        sims.append((np.linalg.norm(v) * np.linalg.norm(w))**alpha * cos)
    return np.mean(sims)

def mag_angle_similarity(v, L, alpha=0.5, eps=1e-8):
    sims = []
    for w in L:
        nv, nw = np.linalg.norm(v), np.linalg.norm(w)
        if nv < eps or nw < eps:
            return 0.0
        cos_sim = np.dot(v, w) / (nv * nw)
        length_diff = abs(nv - nw) / max(nv, nw, eps)
        sims.append(cos_sim * np.exp(-alpha * length_diff))
    return np.mean(sims)

def magnitude_weighted_cosine(v, L, eps=1e-8):
    # v, w = np.array(v), np.array(w)
    sims = []
    for w in L:
        nv, nw = np.linalg.norm(v), np.linalg.norm(w)
        if nv < eps or nw < eps:
            return 0.0
        cos_sim = np.dot(v, w) / (nv * nw)
        mag_ratio = min(nv, nw) / max(nv, nw)
        sims.append(cos_sim * mag_ratio)
    return np.mean(sims)

def dual_threshold_similarity(v, w, angle_thresh=0.94, mag_thresh=0.7):
    v, w = np.array(v), np.array(w)
    nv, nw = np.linalg.norm(v), np.linalg.norm(w)
    if nv == 0 or nw == 0:
        return False
    cos_angle = np.dot(v, w) / (nv * nw)
    mag_ratio = min(nv, nw) / max(nv, nw)
    return (cos_angle >= angle_thresh) and (mag_ratio >= mag_thresh), (cos_angle, mag_ratio)

def euclidean1_similarity(v, w):
    return abs(w[0]-v[0])+abs(w[1]-v[1])

def euclidean2_similarity(v, w):
    return np.linalg.norm(w - v)

def normalized_euclidean1_similarity(v, w):
    return (abs(w[0]-v[0])+abs(w[1]-v[1]))/(np.linalg.norm(v, 1) + np.linalg.norm(w, 1))

def normalized_euclidean2_similarity(v, w):
    return np.linalg.norm(w - v)/(np.linalg.norm(v) + np.linalg.norm(w))

def length_ratio(v,w):
    return min(np.linalg.norm(v), np.linalg.norm(w))/max(np.linalg.norm(v),np.linalg.norm(w))

# 判断“加入v是否破坏整体相似度”
# def compare_similar(m, L, thre, method=mag_angle_similarity):
# def compare_similar(m, L, method=magnitude_weighted_cosine):
# def compare_similar(m, L, method=normalized_euclidean1_similarity):
# def compare_similar(m, L, method=normalized_euclidean2_similarity):
def compare_similar(m, L, method=length_ratio):
# def compare_similar(m, L, alpha, method=weighted_cosine):
    from itertools import combinations

    def avg_pairwise(L, method):
        pairs = combinations(L, 2)
        return np.mean([method(a, b) for a,b in pairs])

    a = np.array(((m[2][0]-m[1][0])*X_SEP+(m[2][2]-m[1][2])*SITE_SEP,(m[2][1]-m[1][1])*Y_SEP))
    # vec_list = []
    # for w in L:
    #     vec_list.append(np.array(((w[2][0]-w[1][0])*X_SEP+(w[2][2]-w[1][2])*SITE_SEP,(w[2][1]-w[1][1])*Y_SEP)))

    # before = avg_pairwise(vec_list, method)
    # after  = avg_pairwise(np.vstack([vec_list, a]), method)
    # if after < 0:
        # print(m, L)
    return np.mean([method(a, np.array(((w[2][0]-w[1][0])*X_SEP+(w[2][2]-w[1][2])*SITE_SEP,(w[2][1]-w[1][1])*Y_SEP))) for w in L])
def mean_sim_numpy(v, L):
    v = np.array(v)
    L = np.array(L)
    dot_products = L @ v
    norm_v = np.linalg.norm(v)
    norm_L = np.linalg.norm(L, axis=1)
    sims = dot_products / (norm_v * norm_L)
    return np.mean(sims)

def count_paths_and_loops_deg2(G: nx.DiGraph, count=None):
    if count is None:
        count = defaultdict(int)
    elif not isinstance(count, defaultdict):
        count = defaultdict(int, count)

    visited = set()

    def traverse_chain(start):
        path = [start]
        cur = start
        while True:
            succs = list(G.successors(cur))
            if len(succs) != 1:  # 无后继或分叉 -> 结束
                break
            nxt = succs[0]
            if nxt in path:  # 成环
                cycle_len = len(path) - path.index(nxt)
                count[cycle_len] += 1
                visited.update(path)
                return
            if nxt in visited:
                break
            path.append(nxt)
            cur = nxt
        # 非环链
        count[len(path) - 1] += 1
        visited.update(path)

    # 先从所有 indegree = 0 的节点出发（非环起点）
    for node in G.nodes:
        if node not in visited and G.in_degree(node) == 0:
            traverse_chain(node)

    # 处理剩下的环
    for node in G.nodes:
        if node not in visited:
            traverse_chain(node)

    return count

def find_threshold_key(d, threshold=0.7):
    sorted_items = sorted(d.items()) 
    
    total = sum(v for _, v in sorted_items)
    target = threshold * total
    
    cumulative = 0
    for key, value in sorted_items:
        if cumulative > target:
            return key
        cumulative += value

    return sorted_items[-1][0]

Fidelity_Atom_Transfer = 0.999
MUS_PER_FRM = 8
X_SEP = 19
Y_SEP = 15
SITE_SEP = 1.5


def check_conflict(pre_move, move, dim):
    src0 = move[1][dim]
    src1 = pre_move[1][dim]
    dst0 = move[2][dim]
    dst1 = pre_move[2][dim]

    if dim == 1:
        distance_src = src0 - src1
        distance_dst = dst0 - dst1
    else:
        distance_src = src0 - src1 + 0.1 * (move[1][2] - pre_move[1][2])
        distance_dst = dst0 - dst1 + 0.1 * (move[2][2] - pre_move[2][2])

    if distance_src == 0:
        dir_src = 0
    else:
        dir_src = distance_src // abs(distance_src)

    if distance_dst == 0:
        dir_dst = 0
    else:
        dir_dst = distance_dst // abs(distance_dst)

    return dir_dst != dir_src

def find_transfer_loc(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index):
    pos_x = move[1][0]
    pos_y = move[1][1]
    pos_find_flag = False
    for r in range(20 * Row):
        for i in range(min(r + 1, Row)):
            j = r - i
            # print("i, j", i, j)
            for a in [-1, 1]:
                for b in [-1, 1]:
                    npos_x = pos_x + a * i
                    npos_y = pos_y + b * j
                    if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < Row and len(empty_space[(npos_x, npos_y)]) < location_size:
                        if len(initial_space[(npos_x, npos_y)]) < location_size:
                            index_set = set(range(location_size))
                            for q in initial_space[(npos_x, npos_y)]:
                                index_set.discard(location_index[q])
                            for q in empty_space[(npos_x, npos_y)]:
                                index_set.discard(target_location_index[q])
                            for m in extra_move:
                                if m[2][0] == npos_x and m[2][1] == npos_y:
                                    index_set.discard(m[2][2])
                            if not index_set:
                                continue
                            loc = (npos_x, npos_y, index_set.pop())
                            pos_find_flag = True
                            break
                if pos_find_flag:
                    break
            if pos_find_flag:
                break
        if pos_find_flag:
            break
    return loc

def find_transfer_loc2(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index):
    pos_x = move[1][0]
    pos_y = move[1][1]
    dest_x = move[2][0]
    dest_y = move[2][1]
    pos_find_flag = False
    a = 1 if dest_x - pos_x > 0 or (dest_x - pos_x == 0 and pos_x < Row/2) else -1
    b = 1 if dest_y - pos_y > 0 or (dest_y - pos_y == 0 and pos_y < Row/2) else -1
    pos_x = move[1][0]
    pos_y = move[1][1]
    pos_find_flag = False
    for r in range(20 * Row):
        for i in range(min(r + 1, Row)):
            j = r - i
            npos_x = pos_x + a * i
            npos_y = pos_y + b * j
            if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < Row and len(empty_space[(npos_x, npos_y)]) < location_size:
                if len(initial_space[(npos_x, npos_y)]) < location_size:
                    index_set = set(range(location_size))
                    for q in initial_space[(npos_x, npos_y)]:
                        index_set.discard(location_index[q])
                    for q in empty_space[(npos_x, npos_y)]:
                        index_set.discard(target_location_index[q])
                    for m in extra_move:
                        if m[2][0] == npos_x and m[2][1] == npos_y:
                            index_set.discard(m[2][2])
                    if not index_set:
                        continue
                    loc = (npos_x, npos_y, index_set.pop())
                    pos_find_flag = True
                    break
        if pos_find_flag:
            break
    return loc

def find_transfer_loc_2qubit(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index):
    pos_x = move[1][0]
    pos_y = move[1][1]
    pos_x2 = move[2][0]
    pos_y2 = move[2][1]
    # make the distance for the two qubits as small as possible
    a = 1 if pos_x2 - pos_x > 0 or (pos_x2 - pos_x == 0 and pos_x < Row/2) else -1
    b = 1 if pos_y2 - pos_y > 0 or (pos_y2 - pos_y == 0 and pos_y < Row/2) else -1
    for r in range(20 * Row):
        for i in range(min(r + 1, Row)):
            j = r - i
            # for a in [-1, 1]:
            #     for b in [-1, 1]:
            # print("i, j", i, j)
            npos_x = pos_x + a * i
            npos_y = pos_y + b * j
            if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < Row and len(empty_space[(npos_x, npos_y)]) < location_size-1:
                if len(initial_space[(npos_x, npos_y)]) < location_size-1:
                    index_set = set(range(location_size))
                    for q in initial_space[(npos_x, npos_y)]:
                        index_set.discard(location_index[q])
                    for q in empty_space[(npos_x, npos_y)]:
                        index_set.discard(target_location_index[q])
                    for m in extra_move:
                        if m[2][0] == npos_x and m[2][1] == npos_y:
                            index_set.discard(m[2][2])
                    if len(index_set)<2:
                        continue
                    loc = (npos_x, npos_y, index_set.pop())
                    loc2 = (npos_x, npos_y, index_set.pop())
                    return loc, loc2
    for r in range(20 * Row):
        for i in range(min(r + 1, Row)):
            j = r - i
            for a in [-1, 1]:
                for b in [-1, 1]:
            # print("i, j", i, j)
                    npos_x = pos_x + a * i
                    npos_y = pos_y + b * j
                    if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < Row and len(empty_space[(npos_x, npos_y)]) < location_size-1:
                        if len(initial_space[(npos_x, npos_y)]) < location_size-1:
                            index_set = set(range(location_size))
                            for q in initial_space[(npos_x, npos_y)]:
                                index_set.discard(location_index[q])
                            for q in empty_space[(npos_x, npos_y)]:
                                index_set.discard(target_location_index[q])
                            for m in extra_move:
                                if m[2][0] == npos_x and m[2][1] == npos_y:
                                    index_set.discard(m[2][2])
                            if len(index_set)<2:
                                continue
                            loc = (npos_x, npos_y, index_set.pop())
                            loc2 = (npos_x, npos_y, index_set.pop())
                            return loc, loc2            
    
    return move[1], move[2]

def find_transfer_loc_1qubit(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index):
    pos_x = move[1][0]
    pos_y = move[1][1]
    pos_find_flag = False
    for r in range(20 * Row):
        for i in range(min(r + 1, Row)):
            j = r - i
            # print("i, j", i, j)
            for a in [-1, 1]:
                for b in [-1, 1]:
                    npos_x = pos_x + a * i
                    npos_y = pos_y + b * j
                    if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < Row and len(empty_space[(npos_x, npos_y)]) == 0:
                        if len(initial_space[(npos_x, npos_y)]) < location_size:
                            index_set = set(range(location_size))
                            for q in initial_space[(npos_x, npos_y)]:
                                index_set.discard(location_index[q])
                            for m in extra_move:
                                if m[2][0] == npos_x and m[2][1] == npos_y:
                                    index_set.discard(m[2][2])
                            if not index_set:
                                continue
                            loc = (npos_x, npos_y, index_set.pop())
                            pos_find_flag = True
                            break
                if pos_find_flag:
                    break
            if pos_find_flag:
                break
        if pos_find_flag:
            break
    return loc

def check_movelist_conflict(move_list, new_move):
    for move in move_list:
        if check_conflict(move, new_move, 0) or check_conflict(move, new_move, 1):
            return True
    return False

def split_move(empty_space, initial_space, move, parallel_move_groups, extra_move, Row, location_size, location_index, target_location_index, release_index, iter_num):
    if iter_num==0:
        return None
    pos_x = move[1][0]
    pos_y = move[1][1]
    dest_x = move[2][0]
    dest_y = move[2][1]
    comp_group = {}
    # a = 1 if dest_x - pos_x > 0 or (dest_x - pos_x == 0 and pos_x < Row/2) else -1
    # b = 1 if dest_y - pos_y > 0 or (dest_y - pos_y == 0 and pos_y < Row/2) else -1
    for r in range(20 * Row):
        for i in range(min(r + 1, Row)):
            j = r - i
            # print("i, j", i, j)
            for a in [-1, 1]:
                for b in [-1, 1]:
                    npos_x = pos_x + a * i
                    npos_y = pos_y + b * j
                    if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < Row and len(empty_space[(npos_x, npos_y)]) < location_size:
                        if len(initial_space[(npos_x, npos_y)]) < location_size:
                            index_set = set(range(location_size + 1))
                            for q in initial_space[(npos_x, npos_y)]:
                                index_set.discard(location_index[q])
                            for q in empty_space[(npos_x, npos_y)]:
                                index_set.discard(target_location_index[q])
                            for m in extra_move:
                                if m[2][0] == npos_x and m[2][1] == npos_y:
                                    index_set.discard(m[2][2])
                            if not index_set:
                                continue
                            loc = (npos_x, npos_y, index_set.pop())
                            new_move = (move[0], move[1], loc)
                            new_move2 = (move[0], loc, move[2])

                            for i1 in range(len(parallel_move_groups)-1):
                                if not check_movelist_conflict(parallel_move_groups[i1], new_move):
                                    comp_group[loc] = comp_group.get(loc, []) + [i1]
                                    for i2 in range(max(release_index, i1+1), len(parallel_move_groups)):
                                        if not check_movelist_conflict(parallel_move_groups[i2], new_move2):
                                            return ([i2, i1], [new_move2, new_move])
    # for access_pos in comp_group.keys():
    #     for i in comp_group[access_pos]:
    #         if (len(parallel_move_groups)-i > 2):
    #             switcher = split_move(empty_space, initial_space, (move[0], access_pos, move[2]), parallel_move_groups[(i+1):], extra_move, Row, location_size, location_index, target_location_index, release_index, iter_num-1)
    #             if switcher is not None:
    #                 swth_index, swth_pos = switcher
    #                 for idx in range(len(swth_index)):
    #                     swth_index[idx] = swth_index[idx] + i + 1
    #                 swth_index.append(i)
    #                 swth_pos.append((move[0], move[1], access_pos))
    #                 # swth_index.insert(0,i)
    #                 # swth_pos.insert(0,(move[0], move[1], loc))
    #                 return (swth_index, swth_pos)
    return None

import math
from scipy.optimize import fsolve

def compare_if_split(new_pos, move, parallel_move_groups, move_distance, add_move_num, num_of_node, cost_para, para1, para2):
    def get_distance(move):
        return move_distance[move]
    
    cost_not_split, cost_split = 0, 0
    swth_index, new_move_list = new_pos
    sim_list = []
    for i, m in enumerate(new_move_list):
        aod = parallel_move_groups[swth_index[i]]
        sim_list.append(compare_similar(m,aod))
        aod.sort(reverse = True, key = get_distance)
        cost_not_split += 200*(get_distance(aod[0])/110)**(1/2)
        aod.append(m)
        move_distance[m] = abs((m[2][0]-m[1][0])*X_SEP+(m[2][2]-m[1][2])*SITE_SEP)+abs(m[2][1]-m[1][1])*Y_SEP
        aod.sort(reverse = True, key = get_distance)
        cost_split += 200*(get_distance(aod[0])/110)**(1/2)
    a = len(parallel_move_groups) + 1
    b = add_move_num 

    def f(x):
        return b * math.log(x) + x**(-a) - x
    x0 = 0.2

    # sol = fsolve(f, x0)[0]
    # sol, info, ier, msg = fsolve(f, x0, full_output=True)

    # if ier != 1 or sol >1 or sol <0:
    #     sol = 0.85
    # else:
    #     print(sol)
    # print(sol)
    # num_aod_est = (sol-sol**(-num_of_node))/math.log(sol)
    # cost_not_split += (200*(get_distance(move)/110)**(1/2) + 2 * MUS_PER_FRM ) / (num_of_node/num_aod_est)

    cost_not_split += (200*(get_distance(move)/110)**(1/2) + 2 * MUS_PER_FRM ) / (1+(num_of_node-add_move_num)/(cost_para*a))
    # print(num_of_node/num_aod_est)
    # print((1+(num_of_node-add_move_num)/(2*a)))
    if (cost_split<=cost_not_split):
        # print(sim_list)
        # return np.mean(np.array(sim_list)) > -1 * get_distance(move)**0.5
        # return True
        if np.mean(np.array(sim_list)) > para1:
            return True
        else:
            # print(sim_list)
            return False

def find_chains_deg1(G: nx.DiGraph, min_length):
    visited = set()
    chains = []
    # print("min length", min_length)
    # print("edges", G.edges())
    # print("nodes", G.nodes())
    for node in G.nodes:
        if node in visited:
            continue

        # 找到当前链/环的起点
        # start = node
        # preds = list(G.predecessors(start))
        # if len(preds) == 1 and preds[0] not in visited:
        #     start = preds[0]
        # else:
        #     break
        if G.in_degree(node) != 0:
            continue
        
        # 从起点顺着 successor 走
        chain = [node]
        cur = node
        while True:
            visited.add(cur)
            succs = list(G.successors(cur))
            if len(succs) != 1:
                break
            nxt = succs[0]
            if nxt in chain:  # 成环
                cycle_start = chain.index(nxt)
                print("error")
                chain = chain[cycle_start:]  # 只保留环部分
                break
            chain.append(nxt)
            cur = nxt

        if len(chain) - 1 >= min_length:
            chains.append(chain)

    return chains

def update_dest(move, empty_space, initial_space, extra_move, Row, location_size, location_index, target_location_index, change_dest, move_distance, dependency_graph):
    q = move[0]
    src = move[1]
    dest = move[2]
    q_list = copy.deepcopy(empty_space[(dest[0], dest[1])])
    if len(empty_space[(dest[0], dest[1])])==2:
        q_list.remove(q)
        q2 = q_list[0]
        loc, loc2 = find_transfer_loc_2qubit(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index)
        new_move = (q, src, loc)
        new_move2 = (q2, (dest[0], dest[1], (dest[2]+1)%2), loc2)
        target_location_index[q] = loc[2]
        target_location_index[q2] = loc2[2]
        # print(empty_space[(dest[0], dest[1])])
        # for key, qubit in empty_space.items():
        #     if len(qubit)!=0:
        #         print(key, qubit)
        # print("update")
        empty_space[(dest[0], dest[1])].remove(q)
        empty_space[(dest[0], dest[1])].remove(q2)
        empty_space[(loc[0], loc[1])].append(q)
        empty_space[(loc2[0], loc2[1])].append(q2)
        # for key, qubit in empty_space.items():
        #     if len(qubit)!=0:
        #         print(key, qubit)
        change_dest[q] = (loc[0], loc[1])
        change_dest[q2] = (loc2[0], loc2[1])
        move2 = []
        for n in dependency_graph.nodes():
            if n[0] == q2:
                move2.append(n)
        # print("change_dest", change_dest)
        move_distance[new_move] = abs((loc[0] - src[0]) * X_SEP + (loc[2] - src[2]) * SITE_SEP) + abs((loc[1] - src[1]) * Y_SEP)
        move_distance[new_move2] = abs((loc2[0] - dest[0]) * X_SEP + (loc2[2] - dest[2] - 1) % 2 * SITE_SEP) + abs((loc2[1] - dest[1]) * Y_SEP)
        dependency_graph.add_node(new_move)
        dependency_graph.add_node(new_move2)
        for suc in dependency_graph.successors(move):
            dependency_graph.add_edge(new_move, suc)
        if len(move2)>0:
            for suc in dependency_graph.successors(move2[0]):
                dependency_graph.add_edge(new_move2, suc)
            dependency_graph.remove_node(move2[0])
        dependency_graph.remove_node(move)
    else:
        loc= find_transfer_loc_1qubit(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index)
        new_move = (q, src, loc)
        # print(move, new_move)
        target_location_index[q] = loc[2]
        # for key, qubit in empty_space.items():
        #     if len(qubit)!=0:
        #         print(key, qubit)
        empty_space[(dest[0], dest[1])].remove(q)
        empty_space[(loc[0], loc[1])].append(q)
        # print("update")
        # for key, qubit in empty_space.items():
        #     if len(qubit)!=0:
        #         print(key, qubit)
        change_dest[q] = (loc[0], loc[1])
        move_distance[new_move] = abs((loc[0] - src[0]) * X_SEP + (loc[2] - src[2]) * SITE_SEP) + abs((loc[1] - src[1]) * Y_SEP)
        dependency_graph.add_node(new_move)
        for suc in dependency_graph.successors(move):
            dependency_graph.add_edge(new_move, suc)
        dependency_graph.remove_node(move)
    return empty_space, target_location_index, change_dest, move_distance, dependency_graph


def coll_moves_scheduler(empty_space, initial_space, n, Row, move_distance, move_group, num_aod, move_in_qubits, move_out_qubits, 
                         qubits_not_in_storage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration, list_movement_duration, 
                         num_movement_stage, location_index, target_location_index, location_size, method, count_sum, loop_num, split_succ, 
                         split_fail, cost_para, para1, para2):
    def get_distance(move):
        # return conflict_graph.nodes[move]['move_distance']
        return move_distance[move]

    moves = move_group
    # print(moves)
    # for m1 in moves:
    #     if m1[2][2] == 1:
    #         for m2 in moves:
    #             if m2[0]!=m1[0] and m2[2] == m1[2]:
    #                 print("error")
    #                 print(moves)
    change_dest = {}
    ready_moves = []
    parallel_move_groups = []
    compatible_index = {}
    # moves.sort(key = get_distance)
    # print(moves)

    # make dependency graph and implement baseline solution (if exists loop, change destination of one move)
    dependency_graph = nx.DiGraph()
    for move in moves:
        dependency_graph.add_node(move)
    for i, move in enumerate(moves):
        for j, other_move in enumerate(moves):
            if i != j:
                if move[2] == other_move[1]:
                    dependency_graph.add_edge(other_move, move)
    count_sum = count_paths_and_loops_deg2(dependency_graph, count_sum)
    # print(count)

    loops = list(nx.simple_cycles(dependency_graph))
    loop_num += len(loops)

    # print("graph",dependency_graph.edges())
    # for k in initial_space.keys():
    #     if len(initial_space[k])!=0:
    #         print(k,initial_space[k])
    # print(initial_space)
    # print(empty_space)
    extra_move = []
    move_in_loop = []

    ###################################################################################
    # trivial task split
    if method in ['base', 'move_split', "break_chains"]:
    # if 'change_dest' not in method:
        for l in loops:
            # print("loop", l)
            max_dist = 0
            for m in l:
                if get_distance(m) > max_dist:
                    max_dist = get_distance(m)
                    move = m
            move_in_loop.append(move)
        if len(move_in_loop) != 0:
            pop_idx_list = []
            for i in range(len(move_in_loop)):
                dest1 = (move_in_loop[i][2][0], move_in_loop[i][2][1])
                for j in range(i+1, len(move_in_loop)):
                    dest2 = (move_in_loop[j][2][0], move_in_loop[j][2][1])
                    if dest1 == dest2:
                        pop_idx_list.append(i)
            for i in pop_idx_list[::-1]:
                move_in_loop.pop(i)
            return empty_space, parallel_move_groups, num_movement_stage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration, list_movement_duration, target_location_index, change_dest, move_in_loop, count_sum, loop_num, split_succ, split_fail
    ####################################################################################
    # 将loop中distance最大的move的目标位置改为其他位置
    # for l in loops:
    #     print("loop", l)
    #     max_dist = 0
    #     for m in l:
    #         if get_distance(m) > max_dist:
    #             max_dist = get_distance(m)
    #             move = m
    #     #try: 最大distance move改成随机
    #     # move = l[0]
    #     q = move[0]
    #     src = move[1]
    #     dest = move[2]
    #     loc= find_transfer_loc2(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index)
    #     new_move = (q, src, loc)
    #     new_move2 = (q, loc, dest)
    #     extra_move.append(new_move)
    #     extra_move.append(new_move2)
    #     move_distance[new_move] = abs((loc[0] - src[0]) * X_SEP + (loc[2] - src[2]) * SITE_SEP) + abs((loc[1] - src[1]) * Y_SEP)
    #     move_distance[new_move2] = abs((loc[0] - dest[0]) * X_SEP + (loc[2] - dest[2]) * SITE_SEP) + abs((loc[1] - dest[1]) * Y_SEP)
    #     dependency_graph.add_node(new_move)
    #     dependency_graph.add_node(new_move2)
    #     for suc in dependency_graph.successors(move):
    #         dependency_graph.add_edge(new_move, suc)
    #     for pre in dependency_graph.predecessors(move):
    #         dependency_graph.add_edge(pre, new_move2)
    #     dependency_graph.remove_node(move)
    #     for m in dependency_graph.nodes():
    #         if m[2] == src:
    #             dependency_graph.add_edge(new_move, m)
    #         if m[1] == dest:
    #             dependency_graph.add_edge(m, new_move2)
    #####################################################################################

    ####################################################################################
    # 将loop中distance最大的move两个qubit一起移到其他位置
    else:
        for l in loops:
            # print("loop", l)
            max_dist = 0
            for m in l:
                if get_distance(m) > max_dist:
                    max_dist = get_distance(m)
                    move = m
            # print(move)
            empty_space, target_location_index, change_dest, move_distance, dependency_graph = update_dest(move, empty_space, initial_space, extra_move, Row, location_size, location_index, target_location_index, change_dest, move_distance, dependency_graph)
    ####################################################################################

    count = count_paths_and_loops_deg2(dependency_graph)
    sorted(count.items())
    # print(count)
    threshold_length = find_threshold_key(count, 0.7)
    if "break_chains" in method:
        break_chain_move = []
        # print("edge", dependency_graph.edges())
        chains = find_chains_deg1(dependency_graph, max(3,threshold_length))
        # print("chains", chains)
        # while len(chains) > 0:
        # print(chains)
        for c in chains:
            break_idx = int((len(c)-1)/2)
            move = c[break_idx]
            break_chain_move.append(move)
        pop_idx_list = []
        

        if method == "break_chains":
            if len(break_chain_move) != 0:
                # print("move group", move_in_loop)
                return empty_space, parallel_move_groups, num_movement_stage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration, list_movement_duration, target_location_index, change_dest, break_chain_move, count_sum, loop_num, split_succ, split_fail
        # if method == "break_chains+change_dest" or method == "break_chains+change_dest+move_split":
        else:
            for i in range(len(break_chain_move)-1):
                dest1 = (break_chain_move[i][2][0], break_chain_move[i][2][1])
                for j in range(i+1, len(break_chain_move)):
                    dest2 = (break_chain_move[j][2][0], break_chain_move[j][2][1])
                    if dest1 == dest2:
                        pop_idx_list.append(i)
            for i in pop_idx_list[::-1]:
                break_chain_move.pop(i)
        # print(break_chain_move)
            for move in break_chain_move:
                # print("break chain")
                empty_space, target_location_index, change_dest, move_distance, dependency_graph = update_dest(move, empty_space, initial_space, extra_move, Row, location_size, location_index, target_location_index, change_dest, move_distance, dependency_graph)
        break_chain_move = []
        # print("edge", dependency_graph.edges())
        # new_chains = find_chains_deg1(dependency_graph, max(3,threshold_length))
        # if len(new_chains)!=0:
        #     print("new chain",new_chains)
        # chains = new_chains
            
    # 把所有被depend的move和independent的move统一放在一个pool中，然后每次从pool中选出不冲突的move组成一个parallel move group，
    # 若选出的move中有被depend的move，则把它的depend move也加入到pool中，pool中维持所有的move按move distance从小到大排序

    for move in dependency_graph.nodes():
        if dependency_graph.in_degree(move) == 0:
            if move not in ready_moves:
                ready_moves.append(move)

    ready_moves.sort(key = get_distance)

    for m in extra_move:
        compatible_index[m[0]] = 1
    ########################################################################################
    # basic
    # if 'move_split' not in method:
    if method in ['base', 'change_dest', "break_chains", "break_chains+change_dest", 'powermove']:
        while len(ready_moves) > 0:
            move = ready_moves[0]
            ready_moves.remove(move)
            # 把它的depend move也加入到pool中
            for succ in dependency_graph.successors(move):
                ready_moves.append(succ)
            ready_moves.sort(key = get_distance)
            flag = False
            release_index = len(parallel_move_groups)
            if dependency_graph.in_degree(move) == 0:
                release_index = 0
            for i in range(len(parallel_move_groups)):
                pg = parallel_move_groups[i]
                if release_index > i:
                    for pre_move in dependency_graph.predecessors(move):
                        if pre_move in pg:
                            if move in extra_move:
                                release_index = max(i, compatible_index[move[0]])
                            else:
                                release_index = i
                            break
                if release_index > i:
                    continue
                if not check_movelist_conflict(pg, move):
                    parallel_move_groups[i].append(move)
                    if move in extra_move:
                        compatible_index[move[0]] = i+1
                    flag = True
                    break

            if not flag:
                parallel_move_groups.append([move])
                if move in extra_move:
                    compatible_index[move[0]] = len(parallel_move_groups) + 1
    ########################################################################################

    ########################################################################################
    # greedily 拆分 move
    else:
        # left_move_num = dependency_graph.number_of_nodes()
        added_move_num = 0
        while len(ready_moves) > 0:
            move = ready_moves[0]
            ready_moves.remove(move)
            added_move_num+=1
            # left_move_num-=1
            # 把它的depend move也加入到pool中
            for succ in dependency_graph.successors(move):
                if succ not in ready_moves:
                    ready_moves.append(succ)
            ready_moves.sort(key = get_distance)
            flag = False
            release_index = len(parallel_move_groups)
            if dependency_graph.in_degree(move) == 0:
                release_index = 0
            for i in range(len(parallel_move_groups)):
                pg = parallel_move_groups[i]
                if release_index > i:
                    for pre_move in dependency_graph.predecessors(move):
                        if pre_move in pg:
                            if move in extra_move:
                                release_index = max(i, compatible_index[move[0]])
                            else:
                                release_index = i
                            break
                if release_index > i:
                    continue
                if not check_movelist_conflict(pg, move):
                    parallel_move_groups[i].append(move)
                    if move in extra_move:
                        compatible_index[move[0]] = i+1
                    flag = True
                    break
            if not flag:
                # print("try to split", move, release_index, parallel_move_groups)
                new_pos = split_move(empty_space, initial_space, move, parallel_move_groups, extra_move, Row, location_size, location_index, target_location_index, release_index, iter_num=2)
                if new_pos is not None:
                    # print("split move", move, new_pos)
                    # if True:
                    if compare_if_split(new_pos, move, parallel_move_groups, move_distance, added_move_num, dependency_graph.number_of_nodes(), cost_para, para1, para2):
                        split_succ += 1
                        swth_index, new_move_list = new_pos
                        dep_move = move
                        for suc in dependency_graph.successors(move):
                            dependency_graph.add_edge(new_move_list[-1], suc)
                            dep_move = suc
                        if dep_move != move:
                            dependency_graph.remove_edge(move, suc)
                        if move in extra_move:
                            compatible_index[move[0]] = swth_index[0]+1
                        for i, m in enumerate(new_move_list):
                            move_distance[m] = abs((m[2][0]-m[1][0])*X_SEP+(m[2][2]-m[1][2])*SITE_SEP)+abs(m[2][1]-m[1][1])*Y_SEP
                            parallel_move_groups[swth_index[i]].append(m)
                            extra_move.append(m)
                        # for m in dependency_graph.nodes():
                        #     if m[2] == new_move_list[-1][1]:
                        #         dependency_graph.add_edge(new_move_list[-1], m)
                        flag = True
                    else:
                        split_fail += 1
                else:
                    split_fail += 1

            if not flag:
                parallel_move_groups.append([move])
                if move in extra_move:
                    compatible_index[move[0]] = len(parallel_move_groups) + 1
    ########################################################################################

    sum = 0
    for m in parallel_move_groups:
        sum += len(m)
    if sum < len(move_group):
        print("error in scheduling moves")
        print("parallel_move_groups", parallel_move_groups)
        print("move_group", move_group)

    ms_index = 0
    while ms_index < len(parallel_move_groups):
        max_distance = 0
        for i in range(num_aod):
            if ms_index == len(parallel_move_groups):
                break
            ms = parallel_move_groups[ms_index]

            list_active_qubits = []
            for m in ms:
                list_active_qubits.append(m[0])
                if m[0] in move_in_qubits:
                    move_in_qubits.remove(m[0])
                if m[0] in move_out_qubits:
                    move_out_qubits.remove(m[0])
            cir_fidelity_atom_transfer *= pow(Fidelity_Atom_Transfer, len(list_active_qubits))
            for i in range(n):
                if i not in list_active_qubits and ((i in qubits_not_in_storage and i not in move_out_qubits) or i in move_in_qubits):
                    cir_qubit_idle_time[i] = cir_qubit_idle_time[i] + MUS_PER_FRM * 2

            ms.sort(reverse = True, key = get_distance)
            max_distance = max(max_distance, get_distance(ms[0]))
            ms_index += 1
        num_movement_stage += 1
        move_duration = 200*(max_distance /110)**(1/2)
        for i in range(n):
            if (i in qubits_not_in_storage and i not in move_out_qubits) or i in move_in_qubits:
                cir_qubit_idle_time[i] += move_duration
        list_transfer_duration.append(2 * MUS_PER_FRM)
        list_movement_duration.append(move_duration)
    return empty_space, parallel_move_groups, num_movement_stage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration, list_movement_duration, target_location_index, change_dest, move_in_loop, count_sum, loop_num, split_succ, split_fail