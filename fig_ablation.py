from matplotlib.ticker import LogLocator, FuncFormatter
import matplotlib.pyplot as plt
import numpy as np
import json, ast

method_list = ["base", "powermove", "break_chains+change_dest+move_split_double_half_sim6_"]#, "break_chains+change_dest+move_split_double_half_sim6_"]
method_name_list = ['baseline', 'unlimited resource', 'Our method', 'tuned 2']
idx = -1
bench_dic = {'ghz':5, 'cat':6, 'ising':6, 'wstate':5, 'qft':2}
bench_list = ['ghz', 'cat', 'ising', 'wstate', 'qft']
# bench_dic = {'ghz':idx, 'cat':idx, 'ising':idx, 'wstate':idx, 'qft':idx}
N_Bench_list = {'ghz':[23,40,78,127,255], 'qft':[4,18], 'cat':[4,22,35,65,130,260], 'ising':[26,34,42,66,98], 'wstate':[3,27,36,76,118], 'rca':[6, 10, 20, 30, 40, 50, 60, 80, 100, 150, 250]}
# N_Bench_list = {'ghz':[23,40,78,127,255], 'qft':[4,18,29,63,160,320], 'cat':[4,22,35,65,130,260], 'ising':[10,26,34,42,66,98,420], 'wstate':[3,27,36,76,118,380], 'rca':[6, 10, 20, 30, 40, 50, 60, 80, 100, 150, 250]}
qoao_regu_dic = {4:-2, 7:-3}
qoao_rand_dic = {0.2:-3,0.5:-2}
Qaoa_rand_Qubit_List = [10, 20, 30 ,50, 100]
bench_sim_para = {'ghz':0.7, 'cat':0.9, 'ising':0.7, 'wstate':0.8, 'qft':0.8}

def safe_load(line):
    line = line.strip()
    try:
        return json.loads(line)
    except json.JSONDecodeError:
        return ast.literal_eval(line) 
    
method = "break_chains+change_dest+move_split_double_half_sim10_"
# method = "base"
# method = "powermove"
storage_cir_fidelity_list = []
storage_cir_fidelity_1q_gate_list = []
storage_cir_fidelity_2q_gate_list = []
storage_cir_fidelity_2q_gate_for_idle_list = []
storage_cir_fidelity_atom_transfer_list = []
storage_cir_fidelity_coherence_list = []

start_idx = 0
for benchm, bench_index in bench_dic.items():
    fidelity_list = []
    if benchm == 'ising':
        start_idx = 1
    else:
        start_idx = 0
    thre = bench_sim_para[benchm]
    if method == "break_chains+change_dest+move_split_double_half_sim10_":
        method = method + str(thre)
    with open(f"data/compare_{benchm}_{method}.txt", "r", encoding="utf-8") as f:
        for _ in range(6):
            _ = safe_load(f.readline())
        # storage_cir_fidelity_1q_gate_list.append(safe_load(f.readline())[start_idx:len(N_Bench_list[benchm])])
        storage_cir_fidelity_list.append(safe_load(f.readline())[start_idx:bench_index])
        storage_cir_fidelity_1q_gate_list.append(safe_load(f.readline())[start_idx:bench_index])
        storage_cir_fidelity_2q_gate_list.append(safe_load(f.readline())[start_idx:bench_index])
        storage_cir_fidelity_2q_gate_for_idle_list.append(safe_load(f.readline())[start_idx:bench_index])
        storage_cir_fidelity_atom_transfer_list.append(safe_load(f.readline())[start_idx:bench_index])
        storage_cir_fidelity_coherence_list.append(safe_load(f.readline())[start_idx:bench_index])


fig, axs = plt.subplots(1, len(bench_list), sharey=True, figsize=(15, 2.2))
for i in range(len(bench_list)): 
    storage_fidelity_2q_gate_for_idle_list = np.multiply(storage_cir_fidelity_2q_gate_list[i], storage_cir_fidelity_2q_gate_for_idle_list[i])
    storage_fidelity_transfer_list = np.multiply(storage_fidelity_2q_gate_for_idle_list, storage_cir_fidelity_atom_transfer_list[i])
    storage_fidelity_decoherence_list = np.multiply(storage_fidelity_transfer_list, storage_cir_fidelity_coherence_list[i])
    print(benchm, storage_fidelity_decoherence_list, )
    if i == len(bench_list)-1:
        axs[i].fill_between(N_Bench_list[bench_list[i]], storage_fidelity_2q_gate_for_idle_list, 1, color='green', hatch = '.', alpha = 0.3, label = 'two-qubit gate')
        axs[i].fill_between(N_Bench_list[bench_list[i]], storage_fidelity_transfer_list, storage_fidelity_2q_gate_for_idle_list, hatch = '*', color='purple', alpha = 0.3, label = 'transfer')
        axs[i].fill_between(N_Bench_list[bench_list[i]], storage_fidelity_decoherence_list, storage_fidelity_transfer_list, hatch = '\\', color='orange', alpha = 0.3, label = 'decoherence')
    else:
        axs[i].fill_between(N_Bench_list[bench_list[i]], storage_fidelity_2q_gate_for_idle_list, 1, color='green', hatch = '.', alpha = 0.3)
        axs[i].fill_between(N_Bench_list[bench_list[i]], storage_fidelity_transfer_list, storage_fidelity_2q_gate_for_idle_list, hatch = '*', color='purple', alpha = 0.3)
        axs[i].fill_between(N_Bench_list[bench_list[i]], storage_fidelity_decoherence_list, storage_fidelity_transfer_list, hatch = '\\', color='orange', alpha = 0.3)

    axs[i].set_yscale('log')
    axs[i].yaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0], numticks=3))

    axs[i].set_title(f'{bench_list[i]}', fontsize = 18)

    axs[i].set_xlabel('#Qubits', fontsize = 15)

plt.figtext(0.5, -0.22, 'Our method', ha='center', va='center', fontsize=18)
# plt.figtext(0.5, -0.22, 'Baseline', ha='center', va='center', fontsize=18)
# plt.figtext(0.5, -0.22, 'Unlimited resource', ha='center', va='center', fontsize=18)
fig.legend(fontsize = 13, loc="lower right", bbox_to_anchor=(0.83, 0.12), borderaxespad=0.)
# plt.savefig("fig/qfidelity_aba_powermove.png",  bbox_inches='tight')
plt.savefig("fig/qfidelity_aba_tune.png",  bbox_inches='tight')
# plt.savefig("fig/qfidelity_aba_base.png",  bbox_inches='tight')
# plt.show()