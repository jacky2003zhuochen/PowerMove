import matplotlib.pyplot as plt
import numpy as np
import json

# x = [6, 10, 20, 30, 40, 50, 60, 80, 100]
# y00 = [176, 288, 480, 528, 656, 880, 880, 1088, 1600]
# y10 = [176, 288, 432, 544, 656, 896, 880, 1088, 1568]
# y20 = [176, 288, 416, 528, 656, 864, 880, 1088, 1552]
# y30 = [176, 272, 368, 400, 592, 704, 688, 896, 1184]
# y40 = [176, 272, 352, 416, 592, 688, 688, 896, 1104]
# y50 = [176, 272, 416, 448, 608, 720, 704, 896, 1248]

# y0 = [1212.8209899243668, 2333.672909855848, 4001.894247355939, 4790.167400944261, 6606.610026603597, 9503.14517493982, 9148.530009321965, 13179.125176931759, 20711.539051548625]
# # y1 = [1033.6022217794846, 1786.6280995352304, 3168.4089754341326, 4324.891520937293, 5773.98219533513, 8637.678248104052, 8417.286258476486, 12063.67297626926, 18183.34379645612]
# y2 = [1033.6022217794846, 1786.6280995352304, 3015.489403562491, 4240.205385455603, 5773.98219533513, 8296.956659370637, 8417.286258476486, 12063.67297626926, 17911.536746576625]
# y2_no_opt = [1212.8209899243668, 2120.0537654490568, 3518.863869622499, 4861.665052379752, 6606.610026603597, 9175.124701672808, 9148.530009321965, 13179.125176931759, 19979.48432399213]
# y3 = [1033.6022217794846, 1740.5879189178493, 2751.234474231636, 3347.075244771684, 5293.847313042831, 6849.931647807649, 6749.430194348092, 9940.9518732568, 14150.128188643966]
# y4 = [1033.6022217794846, 1740.5879189178493, 2609.0016938833655, 3523.0866554974887, 5293.847313042831, 6932.719862964788, 6749.430194348092, 9940.9518732568, 13064.5053539443]
# y5 = [1033.6022217794846, 1740.5879189178493, 2999.6997342451837, 3660.574493621247, 5371.984437042873, 6769.141889141052, 6931.371515154442, 9948.152412237177, 14507.562493706791]
# y6 = [1209.6022217794846, 2012.5879189178493, 3415.6997342451837, 4108.574493621247, 5885.847313042831, 7489.141889141052, 7437.430194348092, 10836.9518732568, 15709.154015367978]
# # y1 = [1212.8209899243668, 2077.0554956729275, 3654.023140706029, 4962.351187861441, 6606.610026603597, 9418.26018026288, 9148.530009321965, 13179.125176931759, 20353.276142234958]
# plt.plot(x, y0, marker='o', label='trivial task split')
# # plt.plot(x, y1, marker='o', label='trivial break', linestyle="solid")
# plt.plot(x, [a+b for a,b in zip(y2,y20)], marker='o', label='change dest')
# # plt.plot(x, y2_no_opt, marker='o', label='change dest')
# # plt.plot(x, y3, marker='o', label='trivial break + split move', linestyle="solid")
# plt.plot(x, [a+b for a,b in zip(y4,y40)], marker='o', label='change dest + split move')
# # plt.plot(x, [a+b for a,b in zip(y5,y50)], marker='o', label="trivial + move split")
# plt.plot(x, y6, marker='o', label="trivial + move split")




method_list = ['base', "move_split", "change_dest", "change_dest+move_split"]
N_Qubit_Dic = {}
P_List = [0.1, 0.2, 0.3, 0.4, 0.5]
N_Qubit_List = [5, 10, 20, 30 ,50, 100]
# N_Qubit_List = [6, 15, 30, 40, 50, 70, 80]
# N_Qubit_List = [6, 10, 20, 30, 40, 50, 60, 80, 100]
duration_list = {}
# for method in method_list:
#     duration_list[method] = []
# for degree in range(1,11):
for p in P_List:
# for p in [0.5]:
    for method in method_list:
        # path = f"data/qaoa_regular{degree}_no_storage_compare_{method}.txt"
        # path = f"data/qaoa_random_no_storage{p}_compare_{method}.txt"
        path = f"data/qaoa_random_new_qLayout_no_storage{p}_compare_{method}.txt"
        # path = f"data/new_compare_vqe_{method}.txt"
        try:
            with open(path, "r", encoding="utf-8") as f:
                # N_Qubit_Dic[method] = json.loads(f.readline().strip())
                duration_list[method] = json.loads(f.readline().strip())
                # print(first_line)
        except FileNotFoundError:
            continue
        # if method == 'base' or method == 'move_split':
        #     marker = '^'
        # else:
        #     marker = 'v'
        marker = 'o'
        # plt.plot(N_Qubit_Dic[method], duration_list[method], marker=marker, label=f"{method}")
        plt.plot(N_Qubit_List, duration_list[method], marker=marker, label=f"{method}")
    plt.title(f"QAOA Rand")
    plt.legend()
    plt.xlabel("qubit")
    plt.ylabel("duration list")
    plt.savefig(f"fig/new_compare_qaoa_rand{p}_new_qLayout.png")
    plt.clf()



# plt.savefig("compare_move_duration.png")
# plt.show()