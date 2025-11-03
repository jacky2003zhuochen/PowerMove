import matplotlib.pyplot as plt

# data = {0: 108, 1: 93, 2: 99, 3: 54, 4: 45, 5: 40, 6: 14, 7: 15, 8: 12, 9: 13, 10: 11, 11: 5, 12: 7, 13: 1, 14: 1, 15: 2, 16: 4, 17: 4, 18: 1, 20: 3, 23: 2, 30: 1, 35: 1}

data_list = [
    # {0: 108, 1: 44, 2: 22, 3: 13, 4: 5, 5: 8, 6: 3, 7: 3, 8: 2, 9: 1},
    {0: 125, 1: 104, 2: 120, 3: 63, 4: 53, 5: 33, 6: 18, 7: 16, 8: 13, 9: 13, 10: 7, 11: 12, 12: 10, 13: 6, 14: 2, 15: 3, 16: 5, 17: 3, 18: 1, 19: 1, 21: 4, 22: 1, 23: 2, 24: 2, 25: 1, 26: 2, 27: 1, 28: 1},
    {0: 67, 1: 35, 2: 18, 3: 13, 4: 4, 5: 8, 6: 1, 8: 1, 11: 1, 18: 1}
]

titles = [ "QAOA Random0.5-n100", "QAOA Regular7-n80"]

# 创建一行多列子图
fig, axes = plt.subplots(1, len(data_list), figsize=(6, 2.5))

# 避免只有一个图时 axes 不是数组
if len(data_list) == 1:
    axes = [axes]

# 逐个绘制
for i, (data, ax) in enumerate(zip(data_list, axes)):
    keys = list(data.keys())
    values = list(data.values())
    bar_width = 0.6
    bars = ax.bar(keys, values, width=bar_width, color='#1f77b4', alpha=0.8)
    ax.set_title(titles[i], fontsize=10)
    ax.grid(axis='y', linestyle='--', alpha=0.4)
    ax.set_xlabel("Length", fontsize=9)
    if i == 0:
        ax.set_ylabel("Number of chains", fontsize=9)


# 添加总标题
fig.suptitle("Length Distribution of Chains", fontsize=15)

# 调整布局
fig.tight_layout(rect=[0, 0, 1, 0.99])


plt.grid(axis='y', linestyle='--', alpha=0.5)
# plt.show()
plt.savefig("fig/chain_length.png",  bbox_inches='tight')