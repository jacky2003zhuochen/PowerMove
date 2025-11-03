import numpy as np
import json, ast
import matplotlib.pyplot as plt

# 假设有 3 个类别，每组 3 个值
group_num = 5   # 组数
bar_per_group = 4  # 每组柱子数量

# method_list = ["base", "powermove", "break_chains+change_dest+move_split_double_half_sim10_"]#, "break_chains+change_dest+move_split_double_half_sim6_"]
method_list = ["base", "slm3_double_half", "slm4_double_half","break_chains+change_dest+move_split_double_half_sim10_", "break_chains+change_dest+move_split_slm3_double_half", "break_chains+change_dest+move_split_slm4_double_half"]
# method_name_list = ['baseline', 'unlimited resource', 'Our method', 'tuned 2']
method_name_list = ['baseline', 'SLM=3', 'SLM=4', 'tuned SLM2', 'tuned SLM3', 'tuned SLM4']
idx = -1
bench_dic = {'ghz':[-3,-2], 'cat':[-3,-2], 'ising':[-4,-3], 'wstate':[-3,-2], 'qft':[-3,-2]}
# bench_dic = {'ghz':idx, 'cat':idx, 'ising':idx, 'wstate':idx, 'qft':idx}
N_Bench_list = {'ghz':[23,40,78,127,255], 'qft':[4,18,29,63,160,320], 'cat':[4,22,35,65,130,260], 'ising':[10,26,34,42,66,98,420], 'wstate':[3,27,36,76,118,380], 'rca':[6, 10, 20, 30, 40, 50, 60, 80, 100, 150, 250]}
qoao_regu_dic = {4:[-3,-2], 7:[-2,-3]}
qoao_rand_dic = {0.2:-3,0.5:-2}
Qaoa_rand_Qubit_List = [10, 20, 30 ,50, 100]
qaoa_regu_qubit_dic = {4:[30, 40, 50, 60, 80, 100], 7:[10, 20, 30, 50, 100]}

def safe_load(line):
    line = line.strip()
    try:
        return json.loads(line)
    except json.JSONDecodeError:
        return ast.literal_eval(line) 


all_duration_list = []
all_fidelity_list = []
for benchm, bench_index_list in bench_dic.items():
    for bench_index in bench_index_list:
        duration_list = []
        fidelity_list = []
        for method in method_list:
            # if 'break_chains+change_dest+move_split' in method:
            if method == "break_chains+change_dest+move_split_double_half_sim10_":
                for thre in [0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1]:
                    method_thre = method + str(thre)
                    min_duration = 10**20
                    max_fid = 0
                    with open(f"data/compare_{benchm}_{method_thre}.txt", "r", encoding="utf-8") as f:
                        duration_n = safe_load(f.readline())[bench_index]
                        min_duration = min(min_duration, duration_n)
                        for _ in range(5):
                            _ = safe_load(f.readline())
                        fidelity_n = safe_load(f.readline())[bench_index]
                        max_fid = max(max_fid, fidelity_n)
                duration_list.append(min_duration)
                fidelity_list.append(max_fid)
            else:
                with open(f"data/compare_{benchm}_{method}.txt", "r", encoding="utf-8") as f:
                    duration_n = json.loads(f.readline().strip())[bench_index]
                    for _ in range(5):
                        _ = safe_load(f.readline())
                    fidelity_n = safe_load(f.readline())[bench_index]
                duration_list.append(duration_n)
                fidelity_list.append(fidelity_n)
        # if benchm == 'qft':
        #     for i in range(len(duration_list)):
        #         duration_list[i] = duration_list[i]/2
        # if benchm in ['ghz', 'cat', 'ising']:
        #     for i in range(len(duration_list)):
        #         duration_list[i] = duration_list[i]*4
        all_duration_list.append(duration_list)
        all_fidelity_list.append(fidelity_list)

# for degree, bench_index_list in qoao_regu_dic.items():
#     for bench_index in bench_index_list:
#         duration_list = []
#         fidelity_list = []
#         for method in method_list:    
#             if 'break_chains+change_dest+move_split' in method:
#                 for thre in [0.2,0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1]:
#                     method_thre = method + str(thre)
#                     min_duration = 10**20
#                     max_fid = 0
#                     with open(f"data/qaoa_regular{degree}_no_storage_compare_{method_thre}.txt", "r", encoding="utf-8") as f:
#                         # qoao_regu_dic[degree] = json.loads(f.readline().strip())[bench_index]
#                         _ = json.loads(f.readline().strip())[bench_index]
#                         duration_n = json.loads(f.readline().strip())[bench_index]
#                         min_duration = min(min_duration, duration_n)
#                         for _ in range(5):
#                             _ = safe_load(f.readline())
#                         fidelity_n = safe_load(f.readline())[bench_index]
#                         max_fid = max(max_fid, fidelity_n)
#                 duration_list.append(min_duration)
#                 fidelity_list.append(max_fid)
#             else:
#                 with open(f"data/qaoa_regular{degree}_no_storage_compare_{method}.txt", "r", encoding="utf-8") as f:
#                     # qoao_regu_dic[degree] = json.loads(f.readline().strip())[bench_index]
#                     _ = json.loads(f.readline().strip())[bench_index]
#                     duration_n = json.loads(f.readline().strip())[bench_index]
#                     for _ in range(5):
#                         _ = safe_load(f.readline())
#                     fidelity_n = safe_load(f.readline())[bench_index]
#                 fidelity_list.append(fidelity_n)
#                 duration_list.append(duration_n)
#         # if degree == 7:
#         #     for i in range(len(duration_list)):
#         #         duration_list[i] = duration_list[i]/4
#         # if degree == 4:
#         #     for i in range(len(duration_list)):
#         #         duration_list[i] = duration_list[i]/4
#         all_duration_list.append(duration_list)
#         all_fidelity_list.append(fidelity_list)

# for p, rand_index in qoao_rand_dic.items():
#     qoao_rand_dic[p] = Qaoa_rand_Qubit_List[rand_index]
#     fidelity_list = []
#     duration_list = []
#     for method in method_list:    
#         if 'break_chains+change_dest+move_split' in method:
#             for thre in [0.2,0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1]:
#                 method_thre = method + str(thre)
#                 min_duration = 10**20
#                 max_fid = 0
#                 with open(f"data/qaoa_rand_no_storage{p}_compare_{method_thre}.txt", "r", encoding="utf-8") as f:
#                     duration_n = json.loads(f.readline().strip())[rand_index]
#                     min_duration = min(min_duration, duration_n)
#                     for _ in range(5):
#                         _ = safe_load(f.readline())
#                     max_fid = max(max_fid, fidelity_n)
#                     fidelity_n = safe_load(f.readline())[bench_index]
#             duration_list.append(min_duration)
#             fidelity_list.append(max_fid)
#         else:
#             with open(f"data/qaoa_rand_no_storage{p}_compare_{method}.txt", "r", encoding="utf-8") as f:
#                 duration_n = json.loads(f.readline().strip())[rand_index]
#                 for _ in range(5):
#                     _ = safe_load(f.readline())
#                 fidelity_n = safe_load(f.readline())[bench_index]
#             fidelity_list.append(fidelity_n)
#             duration_list.append(duration_n)
#     # if p == 0.2:
#     #     for i in range(len(duration_list)):
#     #         duration_list[i] = duration_list[i]
#     # if p == 0.5:
#     #     for i in range(len(duration_list)):
#     #         duration_list[i] = duration_list[i]/2
#     all_duration_list.append(duration_list)
#     all_fidelity_list.append(fidelity_list)

# print(all_duration_list)
# print(all_fidelity_list)
duration_array = np.array(all_duration_list)
duration_array = duration_array / duration_array[:, [0]]
fidelity_array = np.array(all_fidelity_list)
fidelity_array = fidelity_array / fidelity_array[:, [0]]
print(all_duration_list)
# x 坐标（每组之间留一点间距）
x = np.arange(len(all_duration_list))

# 每个柱的宽度
bar_width = 0.15

# 每组的 3 个柱分别平移
offsets = np.arange(len(method_list)) * bar_width

colors = [
    '#1f77b4',  # 蓝色
    '#ff7f0e',  # 橙色
    '#2ca02c',  # 绿色
    '#d62728',  # 红色
    '#9467bd',  # 紫色
    '#8c564b'   # 棕色
]
# 画柱子
fig, ax = plt.subplots(figsize=(18, 5))
for i in range(len(method_list)):
    ax.bar(x + offsets[i], duration_array[:, i], width=bar_width, color=colors[i], label=method_name_list[i])
    # ax.bar(x + offsets[i], fidelity_array[:, i], width=bar_width, color=colors[i], label=method_name_list[i])

# 美化 x 轴
ax.set_xticks(x + bar_width)
ax.set_xticklabels([str(bench) + '_n' + str(N_Bench_list[bench][bench_idx])
          for bench, bench_idx_list in bench_dic.items()
          for bench_idx in bench_idx_list])
# ax.set_xticklabels(["qaoa_regular "+str(bench)+'_n'+str(bench_idx) for bench, bench_idx in qoao_regu_dic.items()])
# ax.set_xticklabels([str(bench) + '_n' + str(N_Bench_list[bench][bench_idx])
#           for bench, bench_idx_list in bench_dic.items()
#           for bench_idx in bench_idx_list] + ["qRegu"+str(bench)+'_n'+str(qaoa_regu_qubit_dic[bench][bench_idx])
#           for bench, bench_idx_list in qoao_regu_dic.items()
#           for bench_idx in bench_idx_list]+ ["qRand"+str(bench)+'_n'+str(bench_idx) for bench, bench_idx in qoao_rand_dic.items()])
# ax.set_xticklabels("qaoa regular "+str(degree) for degree in range(3,11))
ax.set_xlabel('Benchmark')
ax.set_ylabel('Execution Time')
# ax.set_ylabel('Fidelity')
ax.legend()

plt.tight_layout()
# plt.savefig(f"fig/table_fig_time_regular{idx}.png")
# plt.savefig(f"fig/table_fig_time_nor2.png")
plt.savefig(f"fig/table_fig_time_SLM.png")
# plt.savefig(f"fig/table_fig_fidelity_nor.png")
