import matplotlib.pyplot as plt
import numpy as np
import json


# method_list = ['base', "move_split", "change_dest", "change_dest+move_split"]
# method_list = ['base', 'move_split', 'change_dest', 'change_dest+move_split', "break_chains", "break_chains+change_dest", "break_chains+change_dest+move_split"]
# method_list = ['base', 'move_split',  'change_dest+move_split', "break_chains+change_dest+move_split", 'move_split_double',  'change_dest+move_split_double', "break_chains+change_dest+move_split_double"]
method_list = ['base', 'move_split',  'change_dest+move_split', "break_chains+change_dest+move_split", 'move_split_triple',  'change_dest+move_split_triple', "break_chains+change_dest+move_split_triple"]
# method_list = ['base', "break_chains", "break_chains+change_dest", "break_chains+change_dest+move_split", "double_break_chains", "double_break_chains+change_dest", "double_break_chains+change_dest+move_split"]
# method_list = ['base', "dynamic_min3_ratio07_break_chains", "dynamic_min3_ratio07_break_chains+change_dest", "dynamic_min3_ratio07_break_chains+change_dest+move_split", "dynamic_break_chains", "dynamic_break_chains+change_dest", "dynamic_break_chains+change_dest+move_split"]
# method_list = ['base', "break_chains+change_dest+move_split", "break_chains+change_dest+move_split_sim22", "break_chains+change_dest+move_split_sim2"]
# method_list = ['base', 'move_split', 'change_dest', 'change_dest+move_split', "break_chains", "break_chains+change_dest", "break_chains+change_dest+move_split"]

# method = 'break_chains+change_dest+move_split'
N_Qubit_Dic = {}
P_List = [0.1, 0.2, 0.3, 0.4, 0.5]
# N_Qubit_List = [5, 10, 20, 30 ,50, 100]
# N_Qubit_List = [6, 10, 20, 30, 40, 50, 60, 80, 100]
N_Qubit_List = [10, 30, 50, 70, 100, 150, 200]
# N_Qubit_List = [6, 15, 30, 40, 50, 70, 80]
# N_Qubit_List = [6, 10, 20, 30, 40, 50, 60, 80, 100, 150, 250] #rand
# N_Qubit_List = [10, 26, 34, 42, 66, 98, 420] #Ising
# N_Qubit_List = [4, 22, 35, 65, 130, 260] #cat
# N_Qubit_List = [17]
duration_list = {}
base1 = 0
base2 = 0
opt1 = 0
opt2 = 0

bench_list = ['ghz', 'cat', 'ising', 'wstate', 'qft']#'bwt', #'vqe_uccsd']
# bench_list = ['hs']
# # bench_list = ['adder', 'bv', 'cc', 'dnn', 'ghz', 'knn', 'multiplier', 'qft', 'cat', 'ising', 'qugan', 'square_root', 'swap_test', 'wstate', 'rca', ]
N_Bench_list = {'hs':[10, 30, 50, 70, 100, 150, 200],'adder':[10,28,64,118,433], 'bv':[14,19,30,70,140,280], 'bwt':[21,37,57,97,177], 'cc':[12,32,64,151,301], 'dnn':[8,16,33,51], 'ghz':[23,40,78,127,255], 'knn':[25,31,41,67,129,341], 'multiplier':[15,45,75,350,400], 'qft':[4,18,29,63,160,320], 'cat':[4,22,35,65,130,260], 'ising':[10,26,34,42,66,98,420], 'qugan':[39, 71, 111, 395], 'square_root':[18,45,60], 'swap_test':[25,41,83,115,361], 'vqe_uccsd':[4,6,8,28], 'wstate':[3,27,36,76,118,380], 'rca':[6, 10, 20, 30, 40, 50, 60, 80, 100, 150, 250]}

for benchm in bench_list:
    duration_list = {}
    N_Qubit_List = N_Bench_list[benchm]
    for method in method_list:
        # duration_list[method] = []
    # for degree in range(1,11):
    # for p in P_List:
        # opt1 = 0
        # opt2 = 0
        # path = f"data/compare_{benchm}_break_chains+change_dest+move_split.txt"
        # with open(path, "r", encoding="utf-8") as f:
        #     duration_list['base'] = json.loads(f.readline().strip())
        # plt.plot(N_Qubit_List[-3:], duration_list['base'][-3:], marker='o', label=f"base")
        # base1 = duration_list['base'][-1]
        # opt1 = duration_list['base'][-1]
        # for method in method_list:
        # for min_len in [2,3,4,5]:
        # for thre in [0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6]:
        # for thre in [0.1,0.2,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9]:
        # for thre in [0.2,0.3,0.4]:
        # for thre in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]:
        # for thre in [0.5, 0.9, 1, 1.1, 1.2, 1.4, 1.7, 2, 3, 5, 2.5,7,10,15,20,25]:
        # for thre in [0.1,0.2,0.5,0.7,1,2,3,5,7,10,20]:
        # for thre in [-9,-8,-7,-6,-5,-4,-3, -1]:
            # path = f"data/qaoa_regular{degree}_no_storage_compare_{method}_sim6_{thre}.txt"
            # path = f"data/qaoa_rand_no_storage{p}_compare_{method}_sim5_{thre}.txt"
            # path = f"data/compare_{benchm}_{method}_sim_{thre}.txt"
        path = f"data/compare_{benchm}_{method}.txt"
        # path = f"data/compare_{benchm}_{method}_{min_len}tune{thre}.txt"
        # path = f"data/compare_{benchm}_break_chains+change_dest+move_split_sim_{thre}.txt"
        # path = f"data/qaoa_regular{degree}_no_storage_compare_{method}_{min_len}tune{thre}.txt"
        # path = f"data/qaoa_rand_no_storage{p}_compare_{method}_{thre}_sim2.txt"
        try:
            with open(path, "r", encoding="utf-8") as f:
                # N_Qubit_Dic[method] = json.loads(f.readline().strip())
                # N_Qubit_Dic[str(thre)] = json.loads(f.readline().strip())
                duration_list[method] = json.loads(f.readline().strip())
                # duration_list[str(thre)] = json.loads(f.readline().strip())
                # if len(duration_list[method]) == 0:
                #     break
                if method == 'base':
                    base1 = duration_list[method][-1]
                    base2 = duration_list[method][-2]
                    opt1 = duration_list[method][-1]
                    opt2 = duration_list[method][-2]
                # print(first_line)
        except FileNotFoundError:
            continue
        # if len(duration_list[method]) == 0:
        #     break

        # if thre == 0.1:
        #     min_val = duration_list[str(thre)][-1]
        #     min_idx = thre
        #     min_val_2 = duration_list[str(thre)][-2]
        #     min_idx_2 = thre
        # else:
        #     if duration_list[str(thre)][-1] <= min_val:
        #         min_val = duration_list[str(thre)][-1]
        #         min_idx = thre
        #     if duration_list[str(thre)][-2] <= min_val_2:
        #         min_val_2 = duration_list[str(thre)][-2]
        #         min_idx_2 = thre
        # opt1 = min(opt1,duration_list[str(thre)][-1])
        opt1 = min(opt1,duration_list[method][-1])
        opt2 = min(opt2,duration_list[method][-2])
        # plt.plot(N_Qubit_Dic[method], duration_list[method], marker='o', label=f"{method}")
        # plt.plot(N_Qubit_Dic[str(thre)], duration_list[str(thre)], marker='o', label=f"{str(thre)}")
        # plt.plot(N_Qubit_List[-2:], duration_list[method][-2:], marker='o', label=f"{method}")
        plt.plot(N_Qubit_List, duration_list[method], marker='o', label=f"{method}")
        # plt.plot(N_Qubit_List[-3:], duration_list[str(thre)][-3:], marker='o', label=f"{thre}")
        # plt.plot(N_Qubit_List, duration_list[str(thre)], marker='o', label=f"{thre}")
    # plt.title(f"QAOA Regular{degree}")

    # if base1 == 0:
    #     continue
    # plt.text(0.5, -0.14, f"at q={N_Qubit_List[-1]}, best para is {min_idx}, {min_val}",
    #     ha='center', va='top', transform=plt.gca().transAxes)
    # plt.text(0.5, -0.20, f"at q={N_Qubit_List[-2]}, best para is {min_idx_2}, {min_val_2}",
    #         ha='center', va='top', transform=plt.gca().transAxes)
    plt.text(0.5, -0.14, f"at q={N_Qubit_List[-1]}, improve {int(1000*(base1/opt1-1))/10}%",
        ha='center', va='top', transform=plt.gca().transAxes)
    plt.text(0.5, -0.20, f"at q={N_Qubit_List[-2]}, improve {int(1000*(base2/opt2-1))/10}%",
            ha='center', va='top', transform=plt.gca().transAxes)
    plt.subplots_adjust(bottom=0.2)

    # plt.tight_layout()
    plt.title(f"{benchm}")
    # plt.title(f"hs")
    # plt.title(f"Regular{degree}")
    # plt.title(f"Random{p}")
    plt.legend()
    plt.xlabel("qubit")
    plt.ylabel("duration list")
    # plt.text(1, -0.2, f"at q={N_Qubit_List[-1]}, improve {int(1000*(base1/opt1-1))/10}%", ha='center', va='top')
    # plt.text(1, -0.4, f"at q={N_Qubit_List[-2]}, improve {int(1000*(base2/opt2-1))/10}%", ha='center', va='top')
    plt.savefig(f"fig/compare_{benchm}_triple.png")
    # plt.savefig(f"fig/compare_regular{degree}_sim6.png")
    # plt.savefig(f"fig/compare_regular{degree}_double_break_chain_{min_len}tune.png")
    # plt.savefig(f"fig/compare_rand{p}_sim5.png")
    # plt.savefig(f"fig/compare_random{p}_double_break_chain_{min_len}tune.png")
    # print(base1,base2,opt1,opt2)
    # print(f"at q={N_Qubit_List[-1]}, improve {int(1000*(base1/opt1-1))/10}%")
    # print(f"at q={N_Qubit_List[-2]}, improve {int(1000*(base2/opt2-1))/10}%")
    plt.clf()


    # plt.savefig("compare_move_duration.png")
    # plt.show()