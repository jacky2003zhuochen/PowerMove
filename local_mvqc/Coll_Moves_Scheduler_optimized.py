import networkx as nx
import copy
from collections import defaultdict
import numpy as np
import heapq
from functools import lru_cache

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

def length_ratio_xy(v,w):
    return (min(abs(v[0]), abs(w[0]))/max(abs(v[0]), abs(w[0])) + min(abs(v[1]), abs(w[1]))/max(abs(v[1]), abs(w[1])))/2

def length_ratio_xy2(v,w):
    return (abs(v[0])/ abs(w[0]) + abs(v[1])/ abs(w[1]))/2

def length_ratio_dis_mix(v,w):
    return (abs(w[0]-v[0])+abs(w[1]-v[1]))*(abs(v[0])+abs(v[1]))/(abs(v[0])+abs(v[1]) + abs(w[0]) + abs(w[1]))**2

# 判断"加入v是否破坏整体相似度"
def compare_similar(m, L, method=length_ratio_dis_mix):
    from itertools import combinations

    def avg_pairwise(L, method):
        pairs = combinations(L, 2)
        return np.mean([method(a, b) for a,b in pairs])

    a = np.array(((m[2][0]-m[1][0])*X_SEP+(m[2][2]-m[1][2])*SITE_SEP,(m[2][1]-m[1][1])*Y_SEP))

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
X_SEP = 10
Y_SEP = 10
SITE_SEP = 2
SOTORAGE_SEP= 3
Storage_Y_SEP = 10

def compute_move_distance(source, dest, d=1):
    """Compute zone-aware move distance between source and dest positions.
    Each position is (x, y, slot_index). y >= 0 means compute zone, y < 0 means storage zone."""
    if dest[1] >= 0 and source[1] >= 0:
        return abs((dest[0] - source[0]) * X_SEP + (dest[2] - source[2]) * SITE_SEP) * d + abs(dest[1] - source[1]) * Y_SEP * d
    elif dest[1] >= 0:
        return abs((dest[0]*4 - source[0]) * X_SEP/4 + (dest[2] - source[2]) * SITE_SEP) * d + abs(dest[1]) * Y_SEP * d + Storage_Y_SEP * d + abs(source[1]+1) * SOTORAGE_SEP * d
    elif source[1] >= 0:
        return abs((dest[0] - source[0]*4) * X_SEP/4 + (dest[2] - source[2]) * SITE_SEP) * d + abs(source[1]) * Y_SEP * d + Storage_Y_SEP * d + abs(dest[1]+1) * SOTORAGE_SEP * d
    else:
        return abs((dest[0] - source[0]) * X_SEP/4 + (dest[2] - source[2]) * SITE_SEP) * d + abs(dest[1] - source[1]) * SOTORAGE_SEP * d

def check_conflict(pre_move, move, dim):
    src0 = move[1][dim]
    src1 = pre_move[1][dim]
    dst0 = move[2][dim]
    dst1 = pre_move[2][dim]

    if dim == 1:
        src0 = (src0 + 1) * Y_SEP / Storage_Y_SEP if src0 >= 0 else src0
        src1 = (src1 + 1) * Y_SEP / Storage_Y_SEP if src1 >= 0 else src1
        dst0 = (dst0 + 1) * Y_SEP / Storage_Y_SEP if dst0 >= 0 else dst0
        dst1 = (dst1 + 1) * Y_SEP / Storage_Y_SEP if dst1 >= 0 else dst1
        distance_src = src0 - src1
        distance_dst = dst0 - dst1
    else:
        src0 = src0*4 if src0 >= 0 else src0
        src1 = src1*4 if src1 >= 0 else src1
        dst0 = dst0*4 if dst0 >= 0 else dst0
        dst1 = dst1*4 if dst1 >= 0 else dst1
        distance_src = src0 - src1 + 0.8 * (move[1][2] - pre_move[1][2])
        distance_dst = dst0 - dst1 + 0.8 * (move[2][2] - pre_move[2][2])

    if distance_src == 0:
        dir_src = 0
    else:
        dir_src = distance_src // abs(distance_src)

    if distance_dst == 0:
        dir_dst = 0
    else:
        dir_dst = distance_dst // abs(distance_dst)

    return dir_dst != dir_src

def find_transfer_loc(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index, h, a1):
    pos_x = move[1][0]
    pos_y = move[1][1]
    pos_find_flag = False
    for r in range(20 * Row):
        for i in range(min(r + 1, Row)):
            j = r - i
            for a in [-1, 1]:
                for b in [-1, 1]:
                    npos_x = pos_x + a * i
                    npos_y = pos_y + b * j
                    if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < math.ceil(Row/(a1**2))-h and len(empty_space[(npos_x, npos_y)]) < location_size:
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

def find_transfer_loc2(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index, a1, h):
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
            if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < math.ceil(Row/(a1**2))-h and len(empty_space[(npos_x, npos_y)]) < location_size:
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

def find_transfer_loc_2qubit(empty_space, src2, initial_space, move, extra_move, Row, location_size, location_index, target_location_index, confliction_graph, h, a1):
    pos_x = move[1][0]
    pos_y = move[1][1]
    if pos_y<0:
        pos_x = move[1][0]//4
        pos_y = 0
    pos_x2 = src2[0]
    pos_y2 = src2[1]
    if pos_y2<0:
        pos_x2 = src2[0]//4
        pos_y2 = 0
    dest_x = move[2][0]
    dest_y = move[2][1]

    # make the longest length move as small as possble
    for r in range(20 * Row):
        accessible_loc = []
        for i in range(min(r + 1, Row)):
            j = r - i
            for a in [1, -1]:
                for b in [1, -1]:
                    npos_x = (pos_x + pos_x2) // 2 + a * i
                    npos_y = (pos_y + pos_y2) // 2 + b * j

                    if npos_x >= 0 and npos_x < Row and npos_y >= 0 and npos_y < math.ceil(Row/(a1**2))-h and len(empty_space[(npos_x, npos_y)]) == 0: #####################
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
                            accessible_loc.append((loc, loc2))
        if len(accessible_loc) != 0:
            min_conflict = 20*confliction_graph.number_of_nodes()
            min_loc = ((0,0,0),(0,0,0))
            for new_loc in accessible_loc:
                new_move1 = (move[0], move[1], new_loc[0])
                new_move2 = (move[0], (move[2][0], move[2][1], (move[2][2]+1)%2), new_loc[1])
                conflict_num = 0
                for node in confliction_graph.nodes():
                    if check_conflict(node, new_move1, 0) or check_conflict(node, new_move1, 1):
                        conflict_num += 1
                    if check_conflict(node, new_move2, 0) or check_conflict(node, new_move2, 1):
                        conflict_num += 1
                min_conflict = min(min_conflict, conflict_num)
                if min_conflict == conflict_num:
                    min_loc = new_loc
            return min_loc[0], min_loc[1]
    print("skip")
    return move[2], move[2]

def find_transfer_loc_1qubit(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index, confliction_graph, a1):
    pos_x = move[1][0]
    pos_y = move[1][1]
    if pos_y > 0:
        pos_x = move[1][0]*4
    pos_find_flag = False
    for r in range(40 * Row):
        accessible_loc = []
        for i in range(min(r + 1, 4*Row)):
            j = r - i
            for a in [1, -1]:
                npos_x = pos_x + a * i
                npos_y = -j-1
                if npos_x >= 0 and npos_x < 4*Row and npos_y >= -math.ceil(Row/(2*a1**2)) and len(empty_space[(npos_x, npos_y)]) == 0 and len(initial_space[(npos_x, npos_y)]) == 0:
                    loc = (npos_x, npos_y, 0)
                    accessible_loc.append(loc)
        if len(accessible_loc) != 0:
            min_conflict = 10*confliction_graph.number_of_nodes()
            min_loc = (0,0,0)
            for new_loc in accessible_loc:
                new_move = (move[0], move[1], new_loc)
                conflict_num = 0
                for node in confliction_graph.nodes():
                    if check_conflict(node, new_move, 0) or check_conflict(node, new_move, 1):
                        conflict_num += 1
                min_conflict = min(min_conflict, conflict_num)
                if min_conflict == conflict_num:
                    min_loc = new_loc
            return min_loc
    return loc

def find_relay_point(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index, a1, h):
    """Find a relay point for splitting a move into 2 moves: src -> relay -> dest.
    Searches first in compute zone (single slot), then in storage zone.
    Returns (x, y, slot_index) or None if no relay point found."""
    src_x, src_y = move[1][0], move[1][1]

    # Build set of positions occupied by extra_move destinations
    extra_dest_slots = {}
    for m in extra_move:
        key = (m[2][0], m[2][1])
        if key not in extra_dest_slots:
            extra_dest_slots[key] = set()
        extra_dest_slots[key].add(m[2][2])

    # Phase 1: Search compute zone for single-slot availability
    col1 = math.ceil(Row / (a1 ** 2)) - h
    best_dist = float('inf')
    best_loc = None

    # Search outward from source position
    if src_y < 0:
        center_x = src_x // 4
        center_y = 0
    else:
        center_x = src_x
        center_y = src_y

    for r in range(20 * Row):
        found_in_ring = False
        for i in range(min(r + 1, Row)):
            j = r - i
            for da in [1, -1]:
                for db in [1, -1]:
                    npos_x = center_x + da * i
                    npos_y = center_y + db * j
                    if npos_x < 0 or npos_x >= Row or npos_y < 0 or npos_y >= col1:
                        continue
                    pos_key = (npos_x, npos_y)
                    empty_cell = empty_space.get(pos_key, [])
                    initial_cell = initial_space.get(pos_key, [])
                    # Need at least 1 available slot (not 2 like 2-qubit case)
                    if len(empty_cell) < location_size and len(initial_cell) < location_size:
                        index_set = set(range(location_size))
                        for q in initial_cell:
                            index_set.discard(location_index[q])
                        for q in empty_cell:
                            index_set.discard(target_location_index[q])
                        if pos_key in extra_dest_slots:
                            index_set -= extra_dest_slots[pos_key]
                        if index_set:
                            loc = (npos_x, npos_y, next(iter(index_set)))
                            return loc
        if r > 5 * Row:
            break  # Don't search too far in compute zone

    # Phase 2: Search storage zone
    if src_y > 0:
        storage_center_x = src_x * 4
    else:
        storage_center_x = src_x

    col2 = math.ceil(Row / (2 * a1 ** 2))
    for r in range(40 * Row):
        for i in range(min(r + 1, 4 * Row)):
            j = r - i
            for da in [1, -1]:
                npos_x = storage_center_x + da * i
                npos_y = -j - 1
                if npos_x < 0 or npos_x >= 4 * Row or npos_y < -col2:
                    continue
                pos_key = (npos_x, npos_y)
                empty_cell = empty_space.get(pos_key, [])
                initial_cell = initial_space.get(pos_key, [])
                if len(empty_cell) == 0 and len(initial_cell) == 0:
                    # Check extra_move destinations
                    if pos_key in extra_dest_slots and 0 in extra_dest_slots[pos_key]:
                        continue
                    return (npos_x, npos_y, 0)

    return None

def check_movelist_conflict(move_list, new_move):
    if not move_list:
        return False

    Y_RATIO = Y_SEP / Storage_Y_SEP

    nm_src_x = new_move[1][0]
    nm_dst_x = new_move[2][0]
    nm_sx = nm_src_x * 4 if nm_src_x >= 0 else nm_src_x
    nm_dx = nm_dst_x * 4 if nm_dst_x >= 0 else nm_dst_x
    nm_sz = new_move[1][2]
    nm_dz = new_move[2][2]
    nm_sy = new_move[1][1]
    nm_dy = new_move[2][1]
    nm_sy_t = (nm_sy + 1) * Y_RATIO if nm_sy >= 0 else nm_sy
    nm_dy_t = (nm_dy + 1) * Y_RATIO if nm_dy >= 0 else nm_dy

    for pre_move in move_list:
        pm_sx = pre_move[1][0]
        pm_dx = pre_move[2][0]
        pm_sx_t = pm_sx * 4 if pm_sx >= 0 else pm_sx
        pm_dx_t = pm_dx * 4 if pm_dx >= 0 else pm_dx

        dsx = nm_sx - pm_sx_t + 0.8 * (nm_sz - pre_move[1][2])
        ddx = nm_dx - pm_dx_t + 0.8 * (nm_dz - pre_move[2][2])

        if (dsx > 0) != (ddx > 0) or (dsx < 0) != (ddx < 0):
            return True

        pm_sy = pre_move[1][1]
        pm_dy = pre_move[2][1]
        pm_sy_t = (pm_sy + 1) * Y_RATIO if pm_sy >= 0 else pm_sy
        pm_dy_t = (pm_dy + 1) * Y_RATIO if pm_dy >= 0 else pm_dy

        dsy = nm_sy_t - pm_sy_t
        ddy = nm_dy_t - pm_dy_t

        if (dsy > 0) != (ddy > 0) or (dsy < 0) != (ddy < 0):
            return True

    return False

@lru_cache(maxsize=64)
def _get_grid_points_aligned_cached(row: int, h: int, a1: int) -> tuple:
    upper_pts = [
        (i * 12.0, j * 10.0)
        for i in range(row) for j in range(math.ceil(row/(a1**2))-h)
    ]
    lower_ny = math.ceil(row/(2*a1**2))
    lower_pts = [
        (i * 3.0, -10.0 - j * 3.0)
        for i in range(4 * row) for j in range(lower_ny)
    ]
    return tuple(upper_pts + lower_pts)

def get_grid_points_aligned(row: int, h: int, a1: int) -> list:
    return list(_get_grid_points_aligned_cached(row, h, a1))


def split_move(empty_space, initial_space, move, parallel_move_groups, extra_move, Row, location_size, location_index, target_location_index, release_index, iter_num, h, a):
    if iter_num == 0:
        return None

    pos_x, pos_y = move[1][0], move[1][1]
    dest_x, dest_y = move[2][0], move[2][1]
    comp_group = {}

    # 预计算 extra_move（这个保留，因为 extra_move 通常较小，且后续按需查询方便）
    extra_move_positions = {}
    for m in extra_move:
        pos = (m[2][0], m[2][1])
        if pos not in extra_move_positions:
            extra_move_positions[pos] = set()
        extra_move_positions[pos].add(m[2][2])

    # 惰性求值缓存
    position_occupied_cache = {}

    max_r = 20 * Row # 限制搜索范围 (原代码保留变量)
    directions = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    all_indices = set(range(location_size + 1))

    # 优化写法，减少变量重复赋值
    mid_pos_x = (12 * pos_x if pos_y >= 0 else 3 * pos_x) + (12 * dest_x if dest_y >= 0 else 3 * dest_x)
    mid_pos_x /= 2
    mid_pos_y = (10 * pos_y if pos_y >= 0 else -7 - 3 * pos_y) + (10 * dest_y if dest_y >= 0 else -7 - 3 * dest_y)
    mid_pos_y /= 2

    # 第一轮搜索：寻找直接可用的中间位置
    points = get_grid_points_aligned(Row, h, a)
    # 去掉平方项里不必要的 abs() 调用
    points.sort(key=lambda p: abs(p[0] - mid_pos_x) + (p[1] - mid_pos_y)**2)
    if len(points) > 200:
        points = points[:200]

    for p in points:
        npos_x = int(p[0] / 12) if p[1] >= 0 else int(p[0] / 3)
        npos_y = int(p[1] / 10) if p[1] >= 0 else -1 - int((-10 - p[1]) / 3)
        pos_key = (npos_x, npos_y)

        empty_cell = empty_space.get(pos_key, [])
        initial_cell = initial_space.get(pos_key, [])
        size = 1 if p[1] < 0 else 2

        if len(empty_cell) < size and len(initial_cell) < size:
            # 惰性求值：如果该位置还没算过占用情况，才进行计算
            if pos_key not in position_occupied_cache:
                occupied = set()
                # 直接复用上面已经拿到的 initial_cell 和 empty_cell，省去两次 .get()
                for q in initial_cell:
                    occupied.add(location_index.get(q, -1))
                for q in empty_cell:
                    occupied.add(target_location_index.get(q, -1))
                if pos_key in extra_move_positions:
                    occupied.update(extra_move_positions[pos_key])
                position_occupied_cache[pos_key] = occupied
            else:
                occupied = position_occupied_cache[pos_key]

            available_indices = all_indices - occupied

            if available_indices:
                loc = (npos_x, npos_y, next(iter(available_indices)))
                new_move = (move[0], move[1], loc)
                new_move2 = (move[0], loc, move[2])

                # 合并冲突检查：将第一个循环成功的结果缓存，避免第二个循环重复检测新移动方案
                valid_i1s = []
                for i1 in range(len(parallel_move_groups) - 1):
                    if not check_movelist_conflict(parallel_move_groups[i1], new_move):
                        valid_i1s.append(i1)  # 记录通过检测的索引
                        if not check_movelist_conflict(parallel_move_groups[i1+1], new_move2):
                            return ([i1 + 1, i1], [new_move2, new_move])

                # 直接遍历刚刚记录的可用索引，跳过之前必定失败的 check_movelist_conflict
                for i1 in valid_i1s:
                    if loc not in comp_group:
                        comp_group[loc] = []
                    comp_group[loc].append(i1)

                    # 检查第二个移动的冲突
                    if i1 + 2 < len(parallel_move_groups):
                        for i2 in range(max(release_index, i1 + 2), len(parallel_move_groups)):
                            if not check_movelist_conflict(parallel_move_groups[i2], new_move2):
                                return ([i2, i1], [new_move2, new_move])

    # 第二轮搜索：递归尝试
    _max_rec = 10
    _rec_tried = 0
    for access_pos, group_indices in comp_group.items():
        for i in group_indices:
            _rec_tried += 1
            if _rec_tried > _max_rec:
                return comp_group
            if len(parallel_move_groups) - i > 2:
                # 使用切片视图而不是实际切片，避免复制数据
                sub_groups = parallel_move_groups[(i + 1):]
                switcher = split_move(
                    empty_space, initial_space,
                    (move[0], access_pos, move[2]),
                    sub_groups, extra_move, Row, location_size,
                    location_index, target_location_index,
                    release_index, iter_num - 1, h, a
                )

                if isinstance(switcher, tuple):
                    swth_index, swth_pos = switcher
                    # 批量更新索引
                    for idx in range(len(swth_index)):
                        swth_index[idx] += i + 1
                    swth_index.append(i)
                    swth_pos.append((move[0], move[1], access_pos))
                    return swth_index, swth_pos

    return comp_group


import math

def compare_if_split(confliction_graph, new_pos, move, parallel_move_groups, move_distance, add_move_num, num_of_node, cost_para, para1):
    def get_distance(move):
        return move_distance[move]

    # return True
    cost_not_split, cost_split = 0, 0
    swth_index, new_move_list = new_pos
    sim_list = []
    for i, m in enumerate(new_move_list):
        aod = (parallel_move_groups[swth_index[i]]).copy()
        aod.sort(reverse = True, key = get_distance)
        cost_not_split += 200*(get_distance(aod[0])/110)**(1/2)
        aod.append(m)
        source = m[1]
        dest = m[2]
        move_distance[m] = compute_move_distance(source, dest)

        aod.sort(reverse = True, key = get_distance)
        cost_split += 200*(get_distance(aod[0])/110)**(1/2)
    a = len(parallel_move_groups) + 1
    b = add_move_num

    def f(x):
        return b * math.log(x) + x**(-a) - x
    x0 = 0.2

    cost_not_split += (200*(get_distance(move)/110)**(1/2) + 2 * MUS_PER_FRM ) * cost_para
    if (cost_split<=cost_not_split):
        return True
        if np.mean(np.array(sim_list)) < para1:
            return True
        else:
            return False

def compare_if_half_split(confliction_graph,new_pos, move, parallel_move_groups, move_distance, cost_para2, para2):
    def get_distance(move):
        return move_distance[move]

    return None
    cost_not_split = 0
    cost_half_split_dic = {}
    sim_list = []
    for loc, idx_list in new_pos.items():
        for idx in idx_list:
            cost_half_split = 0
            aod = (parallel_move_groups[idx]).copy()
            new_move1 = (move[0], move[1], loc)
            new_move2 = (move[0], loc, move[2])
            aod.append(new_move1)

            move_distance[new_move1] = compute_move_distance(move[1], loc)
            move_distance[new_move2] = compute_move_distance(loc, move[2])

            aod.sort(reverse = True, key = get_distance)
            cost_half_split += 200*(get_distance(aod[0])/110)**(1/2)
            cost_half_split += (200*(get_distance(new_move2)/110)**(1/2) + 2 * MUS_PER_FRM ) * cost_para2
            cost_half_split_dic[(loc,idx)] = cost_half_split
    min_item = min(cost_half_split_dic.items(), key=lambda x: x[1])
    loc, idx = min_item[0]
    aod = (parallel_move_groups[idx]).copy()
    new_move1 = (move[0], move[1], loc)
    new_move2 = (move[0], loc, move[2])
    sim_list.append(compare_similar(new_move1,aod))
    cost_half_split = min_item[1]

    cost_not_split += (200*(get_distance(move)/110)**(1/2) + 2 * MUS_PER_FRM ) * cost_para2
    cost_not_split += 200*(get_distance(aod[0])/110)**(1/2)
    if (cost_half_split<=cost_not_split):
        return loc, idx
        if np.mean(np.array(sim_list)) < para2:
            return loc,idx
        else:
            return None
    return None

def find_chains_deg1(G: nx.DiGraph, min_length):
    visited = set()
    chains = []
    for node in G.nodes:
        if G.in_degree(node) != 0:
            continue

        # 从起点顺着 successor 走
        chain = [node]
        cur = node
        while True:
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

def update_dest(move, empty_space, initial_space, extra_move, Row, location_size, location_index, target_location_index, change_dest, move_distance, dependency_graph, h, a1):
    q, src, dest = move[0], move[1], move[2]
    dest_key = (dest[0], dest[1])

    # 使用局部变量减少字典访问
    dest_empty = empty_space[dest_key]
    q_list = dest_empty.copy()  # 使用浅拷贝替代深拷贝

    if len(dest_empty) == 2:
        q_list.remove(q)
        q2 = q_list[0]
        src2_x, src2_y, src2_idx = dest[0], dest[1], (dest[2] + 1) % 2
        for n in dependency_graph.nodes():
            if n[0] == q2:
                src2_x, src2_y, src2_idx = n[1]
                break

        # 批量查找转移位置
        loc, loc2 = find_transfer_loc_2qubit(empty_space, (src2_x, src2_y), initial_space, move, extra_move, Row, location_size, location_index, target_location_index, dependency_graph, h, a1)
        if loc == move[2] and loc2 == move[2]:
            # Fallback: split original move into 2 moves via relay point
            relay_loc = find_relay_point(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index, a1, h)
            if relay_loc is None:
                # Truly no space available anywhere — keep original skip behavior
                print("skip (no relay point found)")
                return empty_space, target_location_index, change_dest, move_distance, dependency_graph

            relay_key = (relay_loc[0], relay_loc[1])

            # Create 2 new moves: src -> relay -> original dest
            new_move_relay = (q, src, relay_loc)
            new_move_back = (q, relay_loc, dest)

            # Update empty_space: remove q from dest, add to relay
            dest_empty.remove(q)
            empty_space[relay_key].append(q)

            # Update change_dest and target_location_index
            change_dest[q] = relay_key
            target_location_index[q] = relay_loc[2]

            # Compute move distances
            move_distance[new_move_relay] = compute_move_distance(src, relay_loc)
            move_distance[new_move_back] = compute_move_distance(relay_loc, dest)

            # Update dependency graph
            dependency_graph.add_nodes_from([new_move_relay, new_move_back])
            dependency_graph.add_edge(new_move_relay, new_move_back)

            successors = list(dependency_graph.successors(move))
            dependency_graph.add_edges_from([(new_move_back, suc) for suc in successors])

            dependency_graph.remove_node(move)

            # Track as extra moves
            extra_move.append(new_move_relay)
            extra_move.append(new_move_back)

            return empty_space, target_location_index, change_dest, move_distance, dependency_graph

        loc_key = (loc[0], loc[1])
        loc2_key = (loc2[0], loc2[1])

        new_move = (q, src, loc)
        new_move2 = (q2, (src2_x, src2_y, src2_idx), loc2)
        # 批量更新目标位置索引
        target_location_index[q] = loc[2]
        target_location_index[q2] = loc2[2]

        # 批量更新空位状态
        dest_empty.remove(q)
        dest_empty.remove(q2)
        empty_space[loc_key].append(q)
        empty_space[loc2_key].append(q2)

        # 批量更新目标变更
        change_dest[q] = loc_key
        change_dest[q2] = loc2_key

        # 预计算移动距离
        move_distance[new_move] = compute_move_distance(src, loc)
        move_distance[new_move2] = compute_move_distance((src2_x, src2_y, src2_idx), loc2)

        # 优化依赖图更新
        # 找到q2相关的所有移动
        move2_nodes = [n for n in dependency_graph.nodes() if n[0] == q2]

        # 批量添加新节点和边
        dependency_graph.add_nodes_from([new_move, new_move2])

        # 获取后继节点并批量添加边
        successors = list(dependency_graph.successors(move))
        dependency_graph.add_edges_from([(new_move, suc) for suc in successors])

        if move2_nodes:
            move2_node = move2_nodes[0]  # 假设每个量子比特只有一个移动节点
            successors2 = list(dependency_graph.successors(move2_node))
            dependency_graph.add_edges_from([(new_move2, suc) for suc in successors2])
            dependency_graph.remove_node(move2_node)

        dependency_graph.remove_node(move)

    else:
        # 单量子比特情况
        loc = find_transfer_loc_1qubit(empty_space, initial_space, move, extra_move, Row, location_size, location_index, target_location_index, dependency_graph, a1)
        loc_key = (loc[0], loc[1])
        new_move = (q, src, loc)
        # 批量更新状态
        target_location_index[q] = loc[2]
        dest_empty.remove(q)
        empty_space[loc_key].append(q)
        change_dest[q] = loc_key

        # 预计算移动距离
        source = src
        dest = loc
        d = 1
        move_distance[new_move] = abs((dest[0] - source[0]*4) * X_SEP/4 + (dest[2] - source[2]) * SITE_SEP) * d + abs(source[1]) * Y_SEP * d + Storage_Y_SEP * d + abs(dest[1]+1) * SOTORAGE_SEP * d
        dist = abs((loc[0] - src[0]) * X_SEP + (loc[2] - src[2]) * SITE_SEP) + abs((loc[1] - src[1]) * Y_SEP)
        move_distance[new_move] = dist

        # 优化依赖图更新
        dependency_graph.add_node(new_move)
        successors = list(dependency_graph.successors(move))
        dependency_graph.add_edges_from([(new_move, suc) for suc in successors])
        dependency_graph.remove_node(move)

    return empty_space, target_location_index, change_dest, move_distance, dependency_graph

def move_qubit(move_list):
    mq_list = []
    for move in move_list:
        q = move[0]
        mq_list.append(q)
    return mq_list

def coll_moves_scheduler(empty_space, initial_space, n, Row, move_distance, move_group, num_aod, move_in_qubits, move_out_qubits,
                         qubits_not_in_storage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration, list_movement_duration,
                         num_movement_stage, location_index, target_location_index, location_size, method, count_sum, loop_num, split_fail, split_succ, cost_para, para1, cost_para2, para2, thre, h, a1, iter_num=2):
    def get_distance(move):
        return move_distance[move]

    moves = move_group
    change_dest = {}
    ready_moves = []
    parallel_move_groups = []
    compatible_index = {}

    # make dependency graph and implement baseline solution (if exists loop, change destination of one move)
    dependency_graph = nx.DiGraph()
    dependency_graph.add_nodes_from(moves)
    ms_index = 0
    span = (100000,0)

    # 使用更高效的方式构建边
    move_dest_map = {move: move[2] for move in moves}
    edges_to_add = []
    for i, move in enumerate(moves):
        dest = move_dest_map[move]
        for j, other_move in enumerate(moves):
            if i != j and dest == other_move[1]:
                edges_to_add.append((other_move, move))
    dependency_graph.add_edges_from(edges_to_add)
    count_sum = count_paths_and_loops_deg2(dependency_graph, count_sum)

    loops = list(nx.simple_cycles(dependency_graph))
    loop_num += len(loops)

    extra_move = []
    move_in_loop = []

    ###################################################################################
    # trivial task split
    if 'change_dest' not in method:
        for l in loops:
            max_dist = 0
            for m in l:
                if m[2][1] >= 0:
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
            return empty_space, parallel_move_groups, num_movement_stage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration, list_movement_duration, target_location_index, change_dest, move_in_loop, count_sum, loop_num, ms_index, 0, 0, span,0
    ####################################################################################

    ####################################################################################
    # 将loop中distance最大的move两个qubit一起移到其他位置
    else:
        break_loop_move = []
        for l in loops:
            max_dist = 0
            for m in l:
                if m[2][1] < 0:
                    move = m
                    break
                if get_distance(m) > max_dist:
                    max_dist = get_distance(m)
                    move = m
            break_loop_move.append(move)
        pop_idx_list = []
        for i in range(len(break_loop_move)-1):
            dest1 = (break_loop_move[i][2][0], break_loop_move[i][2][1])
            for j in range(i+1, len(break_loop_move)):
                dest2 = (break_loop_move[j][2][0], break_loop_move[j][2][1])
                if dest1 == dest2:
                    pop_idx_list.append(i)
        for i in pop_idx_list[::-1]:
            dependency_graph.remove_node(break_loop_move[i])
            break_loop_move.pop(i)
        for move in break_loop_move:
            empty_space, target_location_index, change_dest, move_distance, dependency_graph = update_dest(move, empty_space, initial_space, extra_move, Row, location_size, location_index, target_location_index, change_dest, move_distance, dependency_graph, h, a1)
    ####################################################################################


    pick_drop_times = 0
    pick_drop_times += dependency_graph.number_of_nodes()

    count = count_paths_and_loops_deg2(dependency_graph)
    sorted(count.items())
    threshold_length = find_threshold_key(count, 0.7)
    thre = max(4,threshold_length)
    if "break_chains" in method:
        break_chain_move = []
        chains = find_chains_deg1(dependency_graph, thre)
        while True:
            if len(chains) == 0:
                break
            for c in chains:
                break_idx = int((len(c)-1)/2)
                move = c[break_idx]
                break_chain_move.append(move)
            pop_idx_list = []
            for i in range(len(break_chain_move)-1):
                dest1 = (break_chain_move[i][2][0], break_chain_move[i][2][1])
                for j in range(i+1, len(break_chain_move)):
                    dest2 = (break_chain_move[j][2][0], break_chain_move[j][2][1])
                    if dest1 == dest2:
                        pop_idx_list.append(i)
            for i in pop_idx_list[::-1]:
                break_chain_move.pop(i)

            if method == "break_chains":
                if len(break_chain_move) != 0:
                    return empty_space, parallel_move_groups, num_movement_stage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration, list_movement_duration, target_location_index, change_dest, break_chain_move, count_sum, loop_num, ms_index, 0, 0, span, 0
            else:
                for move in break_chain_move:
                    empty_space, target_location_index, change_dest, move_distance, dependency_graph = update_dest(move, empty_space, initial_space, extra_move, Row, location_size, location_index, target_location_index, change_dest, move_distance, dependency_graph, h, a1)
                break
            break_chain_move = []
            chains = find_chains_deg1(dependency_graph, thre) # 长度大于等于max(4,threshold_length)都会被break

    count_final = count_paths_and_loops_deg2(dependency_graph)
    sorted(count_final.items())
    max_length = list(count_final.keys())[-1]

    # 把所有被depend的move和independent的move统一放在一个pool中，然后每次从pool中选出不冲突的move组成一个parallel move group，
    # 若选出的move中有被depend的move，则把它的depend move也加入到pool中，pool中维持所有的move按move distance从小到大排序

    for move in dependency_graph.nodes():
        if dependency_graph.in_degree(move) == 0:
            if move not in ready_moves:
                ready_moves.append(move)

    # Convert ready_moves to a heap for O(log n) pop/push instead of O(n log n) sort
    _tie = 0
    def _heap_push(heap, move):
        nonlocal _tie
        heapq.heappush(heap, (get_distance(move), _tie, move))
        _tie += 1

    ready_heap = []
    for m in ready_moves:
        _heap_push(ready_heap, m)
    ready_moves_set = set(ready_moves)

    added_move_num = 0
    extra_move_set = set(extra_move)

    for m in extra_move:
        compatible_index[m[0]] = 1
    ########################################################################################
    # basic
    if 'move_split' not in method:
        while ready_heap:
            _, _, move = heapq.heappop(ready_heap)
            for succ in dependency_graph.successors(move):
                _heap_push(ready_heap, succ)
            flag = False
            release_index = len(parallel_move_groups)
            if dependency_graph.in_degree(move) == 0:
                release_index = 0
            for i in range(len(parallel_move_groups)):
                pg = parallel_move_groups[i]
                if release_index > i:
                    for pre_move in dependency_graph.predecessors(move):
                        if pre_move in pg:
                            if move in extra_move_set:
                                release_index = max(i, compatible_index[move[0]])
                            else:
                                release_index = i
                            break
                if release_index > i:
                    continue
                if not check_movelist_conflict(pg, move):
                    parallel_move_groups[i].append(move)
                    if move in extra_move_set:
                        compatible_index[move[0]] = i+1
                    flag = True
                    break

            if not flag:
                parallel_move_groups.append([move])
                if move in extra_move_set:
                    compatible_index[move[0]] = len(parallel_move_groups) + 1
    ########################################################################################

    ########################################################################################
    # greedily 拆分 move
    else:
        confliction_graph = nx.Graph()
        node_list = list(dependency_graph.nodes())
        for move in node_list:
            confliction_graph.add_node(move)
        for i in range(len(node_list)):
            for j in range(i+1, len(node_list)):
                move = node_list[i]
                other_move = node_list[j]
                if check_conflict(move, other_move, 0) or check_conflict(move, other_move, 1):
                    confliction_graph.add_edge(move, other_move)
        while ready_heap:
            _, _, move = heapq.heappop(ready_heap)
            ready_moves_set.discard(move)
            added_move_num+=1
            for succ in dependency_graph.successors(move):
                if succ not in ready_moves_set:
                    ready_moves_set.add(succ)
                    _heap_push(ready_heap, succ)
            flag = False
            release_index = len(parallel_move_groups)
            if dependency_graph.in_degree(move) == 0:
                release_index = 0
            for i in range(len(parallel_move_groups)):
                pg = parallel_move_groups[i]
                if release_index > i:
                    for pre_move in dependency_graph.predecessors(move):
                        if pre_move in pg:
                            if move in extra_move_set:
                                release_index = max(i, compatible_index[move[0]])
                            else:
                                release_index = i
                            break
                if release_index > i:
                    continue
                if not check_movelist_conflict(pg, move):
                    parallel_move_groups[i].append(move)
                    if move in extra_move_set:
                        compatible_index[move[0]] = i+1
                    flag = True
                    break
            if not flag:
                new_pos = split_move(empty_space, initial_space, move, parallel_move_groups, extra_move, Row, location_size, location_index, target_location_index, release_index, iter_num=iter_num, h=h, a=a1)
                if isinstance(new_pos, tuple) and new_pos:
                    if True:
                        split_succ += 1
                        swth_index, new_move_list = new_pos
                        dep_move = move
                        for suc in dependency_graph.successors(move):
                            dependency_graph.add_edge(new_move_list[-1], suc)
                            dep_move = suc
                        if dep_move != move:
                            dependency_graph.remove_edge(move, suc)
                        if move in extra_move_set:
                            compatible_index[move[0]] = swth_index[0]+1
                        for i, m in enumerate(new_move_list):
                            source = m[1]
                            dest = m[2]
                            move_distance[m] = compute_move_distance(source, dest)
                            parallel_move_groups[swth_index[i]].append(m)
                            extra_move.append(m)
                        flag = True
                    else:
                        split_fail += 1

                elif isinstance(new_pos, dict):
                    if len(new_pos) != 0:
                        loc_idx = compare_if_half_split(confliction_graph, new_pos, move, parallel_move_groups, move_distance, cost_para2, para2)
                        if loc_idx is not None:
                            idx = loc_idx[1]
                            loc = loc_idx[0]
                            if move in extra_move_set:
                                compatible_index[move[0]] = idx+1
                            split_succ += 1
                            new_move1 = (move[0], move[1], loc)
                            new_move2 = (move[0], loc, move[2])
                            move_distance[new_move1] = abs((new_move1[2][0]-new_move1[1][0])*X_SEP+(new_move1[2][2]-new_move1[1][2])*SITE_SEP)+abs(new_move1[2][1]-new_move1[1][1])*Y_SEP
                            move_distance[new_move2] = abs((new_move2[2][0]-new_move2[1][0])*X_SEP+(new_move2[2][2]-new_move2[1][2])*SITE_SEP)+abs(new_move2[2][1]-new_move2[1][1])*Y_SEP
                            parallel_move_groups[idx].append(new_move1)
                            extra_move.append(new_move1)
                            parallel_move_groups.append([new_move2])
                            dep_move = move
                            for suc in dependency_graph.successors(move):
                                dependency_graph.add_edge(new_move1, suc)
                                dep_move = suc
                            if dep_move != move:
                                dependency_graph.remove_edge(move, suc)
                            flag = True
                        else:
                            split_fail += 1
                    else:
                        split_fail += 1

            if not flag:
                parallel_move_groups.append([move])
                if move in extra_move_set:
                    compatible_index[move[0]] = len(parallel_move_groups) + 1

            confliction_graph.remove_node(move)

    ########################################################################################
    sum = 0
    for m in parallel_move_groups:
        sum += len(m)

    if sum < len(move_group):
        print("error in scheduling moves")
        print("parallel_move_groups", parallel_move_groups)
        print("move_group", move_group)

    pre_pick_qubits = []
    pick_qubits = []
    AOD_moveTime_list = []
    AOD_move_num_list = []
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
            cir_fidelity_atom_transfer *= pow(Fidelity_Atom_Transfer, len(list_active_qubits)*2)
            for i in range(n):
                if i not in list_active_qubits:
                    cir_qubit_idle_time[i] = cir_qubit_idle_time[i] + MUS_PER_FRM * 2

            ms.sort(reverse = True, key = get_distance)
            max_distance = max(max_distance, get_distance(ms[0]))
            ms_index += 1

            pre_pick_qubits = pick_qubits
            pick_qubits = [m[0] for m in ms]
            if len(pre_pick_qubits) != 0:
                for q in pick_qubits:
                    if q in pre_pick_qubits:
                        sum -= 1
                        cir_fidelity_atom_transfer /= Fidelity_Atom_Transfer

        num_movement_stage += 1
        move_duration = 200*(max_distance /110)**(1/2)
        span = (min(span[0],move_duration), max(span[1], move_duration))
        AOD_moveTime_list.append(move_duration)
        AOD_move_num_list.append(len(ms))
        for i in range(n):
            if i not in list_active_qubits:
                cir_qubit_idle_time[i] += move_duration
        transf_time = 2 * MUS_PER_FRM
        list_transfer_duration.append(transf_time)
        list_movement_duration.append(move_duration)

    return empty_space, parallel_move_groups, num_movement_stage, cir_qubit_idle_time, cir_fidelity_atom_transfer, list_transfer_duration, list_movement_duration, target_location_index, change_dest, move_in_loop, count_sum, loop_num, ms_index, np.sum(np.array(AOD_moveTime_list)), AOD_move_num_list, span, max_length
