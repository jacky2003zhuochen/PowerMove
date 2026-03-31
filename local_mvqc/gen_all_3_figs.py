"""
Generate all 3 figures with the SAME benchmark settings:
  1. rebuttal_main_fig4_new_ablation_qft63_v3.png  (T_exe / fidelity / comp time)
  2. fig11_zone_arch_new.png  (max chain / avg move per AOD / AOD rounds)
  3. fid_aba_new.png  (2q error / atom transfer / decoherence)

Uses the same benchmarks/params as gen_final_fig_v2.py.
"""
import sys, os, math, json, time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
os.chdir(ROOT)
sys.path.insert(0, os.path.join(ROOT, 'local_mvqc'))
sys.path.insert(1, ROOT)
from Construct_Circuit import *
from PowerMove_optimized import *
sys.path.insert(0, 'benchmarks/benchmarks/new_qv_rcs')
from qv_builder import build_standard_qv_circuit_aligned_cz
from rcs_builder import build_dense_cz_rcs, strip_measurements_and_idle_qubits
from qiskit.qasm3 import load as qasm3_load

# ---------------------------------------------------------------------------
# Loaders (same as gen_final_fig_v2.py)
# ---------------------------------------------------------------------------
def load_qasm(benchm, n):
    base = f'benchmarks/QASMBench-master/{benchm}'
    for p in [f'{base}/{benchm}_n{n}/{benchm}_n{n}_transpiled.qasm', f'{base}/{benchm}_n{n}/{benchm}_n{n}.qasm']:
        if os.path.exists(p):
            try: return get_cz_blocks(transpile(QuantumCircuit.from_qasm_file(p), basis_gates=['u1','u2','u3','cz','id'], optimization_level=2))
            except: pass
    for p in [f'{base}/{benchm}_n{n}/{benchm}_n{n}.qasm', f'{base}/{benchm}_n{n}.qasm']:
        if os.path.exists(p):
            return get_cz_blocks(transpile(qasm3_load(p), basis_gates=['u1','u2','u3','cz','id'], optimization_level=2))

def load_qaoa_regu(n, degree):
    with open(f'benchmarks/benchmarks/qaoa/q{n}_regular_{degree}.txt') as f: return [eval(f.read())]
def load_qaoa_rand(n, p):
    with open(f'benchmarks/benchmarks/qaoa/q{n}_random_{p}.txt') as f: return [eval(f.read())]
def load_qv(n=100, depth=8, seed=2024):
    circ = build_standard_qv_circuit_aligned_cz(n, depth=depth, seed=seed, add_measurements=True, insert_barriers=True)
    return get_cz_blocks(strip_measurements_and_idle_qubits(transpile(circ, basis_gates=['rz','rx','ry','cz'], optimization_level=2, seed_transpiler=seed)))
def load_rcs(n=100, depth=15, seed=23):
    circ = build_dense_cz_rcs(n, depth=depth, seed=seed, add_measurements=True, matching='random', insert_barriers=True)
    return get_cz_blocks(strip_measurements_and_idle_qubits(transpile(circ, basis_gates=['rz','rx','ry','cz'], optimization_level=2, seed_transpiler=seed)))

def run(cz_blocks, n, method, a_val, it, cp, sp, cp2, sp2):
    return mvqc(cz_blocks, math.ceil(math.sqrt(n/2)), n, True, 1, 1, method, cp, sp, cp2, sp2, 0.7, 0, a_val, iter_num=it)

# ---------------------------------------------------------------------------
# Benchmarks (same settings as current fig4 v2)
# (name, loader, n, a, best_iter, cp, sp, cp2, sp2)
# ---------------------------------------------------------------------------
BENCHMARKS = [
    ('cat_260',       lambda: load_qasm('cat',260),     260, 6, 1, 0.5,1.0,0.3,1.0),
    ('ghz_255',       lambda: load_qasm('ghz',255),     255, 6, 1, 0.3,1.0,1.0,1.0),
    ('wstate_380',    lambda: load_qasm('wstate',380),  380, 7, 2, 0.5,0.4,0.3,0.15),
    ('ising_420',     lambda: load_qasm('ising',420),   420, 7, 2, 1.0,1.0,0.3,1.0),
    ('qft_63',        lambda: load_qasm('qft',63),       63, 6, 2, 0.4,0.5,1.5,0.4),
    ('qaoa_regu_100', lambda: load_qaoa_regu(100,8),    100, 8, 2, 0.6,0.45,0.9,0.4),
    ('qv_100',        lambda: load_qv(100,8,2024),        100, 1, 2, 0.4,0.5,1.5,0.4),
    ('rcs_100',       lambda: load_rcs(100,15,23),      100, 1, 2, 0.4,0.5,1.5,0.4),
]
LABELS = ['cat 260','ghz 255','wstate 380','ising 420','qft 63',
          'qaoa 100','qv 100','rcs 100']

# ZAC data
with open('/scr/dataset/jixuan/Chenzhuo/ZAC-main/result/fig11_collected/zac_best_texe_fig12.json') as f:
    zac12 = json.load(f)
with open('/scr/dataset/jixuan/Chenzhuo/ZAC-main/result/fig11_collected/zac_best_texe_fig13.json') as f:
    zac13 = json.load(f)

# ---------------------------------------------------------------------------
# Run all benchmarks and collect full data (with caching)
# ---------------------------------------------------------------------------
CACHE_PATH = 'local_mvqc/all_3_figs_data.json'
USE_CACHE = os.path.exists(CACHE_PATH) and '--rerun' not in sys.argv

if USE_CACHE:
    print(f'\n  Loading cached data from {CACHE_PATH}', flush=True)
    print(f'  (use --rerun flag to force re-running benchmarks)', flush=True)
    with open(CACHE_PATH) as _f:
        data = json.load(_f)
else:
    data = {}

for name, loader, n, a, best_it, cp, sp, cp2, sp2 in BENCHMARKS:
    if name in data and USE_CACHE:
        print(f'\n  {name}: cached (PM={data[name]["base_texe"]:.0f} O={data[name]["orbit_texe"]:.0f})', flush=True)
        continue
    print(f'\n  {name} (n={n}, a={a}, iter={best_it})', flush=True)
    cz_blocks = loader()
    n_cz = sum(len(b) for b in cz_blocks)
    print(f'    {n_cz} CZ gates', flush=True)

    # T_exe/fid/metrics at best_iter
    rb = run(cz_blocks,n,'base',a,best_it,cp,sp,cp2,sp2)
    ros = run(cz_blocks,n,'break_chains+change_dest',a,best_it,cp,sp,cp2,sp2)
    ro = run(cz_blocks,n,'break_chains+change_dest+move_split',a,best_it,cp,sp,cp2,sp2)

    # Comp time at iter=1
    t0=time.time(); run(cz_blocks,n,'base',a,1,cp,sp,cp2,sp2); cb=time.time()-t0
    t0=time.time(); run(cz_blocks,n,'break_chains+change_dest',a,1,cp,sp,cp2,sp2); cos_=time.time()-t0
    t0=time.time(); run(cz_blocks,n,'break_chains+change_dest+move_split',a,1,cp,sp,cp2,sp2); co=time.time()-t0

    d = {
        'a': a, 'best_iter': best_it, 'n': n, 'n_cz': n_cz,
        # T_exe
        'base_texe': rb[0]+rb[1], 'orbs_texe': ros[0]+ros[1], 'orbit_texe': ro[0]+ro[1],
        # Fidelity total
        'base_fid': rb[2], 'orbs_fid': ros[2], 'orbit_fid': ro[2],
        # Fidelity breakdown: exic = fid_2q * fid_2q_idle, trans = fid_transfer, coh = fid_coherence
        'base_exic': rb[4]*rb[5], 'base_trans': rb[6], 'base_coh': rb[7],
        'orbs_exic': ros[4]*ros[5], 'orbs_trans': ros[6], 'orbs_coh': ros[7],
        'orbit_exic': ro[4]*ro[5], 'orbit_trans': ro[6], 'orbit_coh': ro[7],
        # Fig11 metrics
        # PowerM max_chain: max key from count dict (before break chain)
        # ORBIT max_chain: result[15] (after break chain)
        'base_max_chain': max(rb[9].keys()) if rb[9] else 0, 'base_avg_move_aod': rb[13], 'base_aod_rounds': rb[11],
        'orbs_max_chain': ros[15], 'orbs_avg_move_aod': ros[13], 'orbs_aod_rounds': ros[11],
        'orbit_max_chain': ro[15], 'orbit_avg_move_aod': ro[13], 'orbit_aod_rounds': ro[11],
        # Comp time
        'comp_base': cb, 'comp_orbs': cos_, 'comp_orbit': co,
    }
    data[name] = d
    print(f'    PM={d["base_texe"]:.0f} O*={d["orbs_texe"]:.0f} O={d["orbit_texe"]:.0f} sp={d["base_texe"]/d["orbit_texe"]:.3f}x', flush=True)
    print(f'    chain: PM={d["base_max_chain"]} O*={d["orbs_max_chain"]} O={d["orbit_max_chain"]}  '
          f'avg_move: PM={d["base_avg_move_aod"]:.2f} O={d["orbit_avg_move_aod"]:.2f}  '
          f'aod: PM={d["base_aod_rounds"]} O={d["orbit_aod_rounds"]}', flush=True)

# Save data
json.dump(data, open('local_mvqc/all_3_figs_data.json','w'), indent=2)
print('\nData saved to local_mvqc/all_3_figs_data.json', flush=True)

names = [b[0] for b in BENCHMARKS]
OUT = 'rebuttal/rebuttal_fig'

# =========================================================================
# Helper
# =========================================================================
def gm(v):
    v=np.asarray(v,dtype=float); v=v[v>0]; return np.exp(np.mean(np.log(v)))

c_zac='#55A868'; c_pm='#DD8452'; c_os='#C44E52'; c_ob='#4C72B0'

# =========================================================================
# FIGURE 1: rebuttal_main_fig4 (T_exe / f_out / T_comp) — 4 bars
# =========================================================================
print('\n=== Generating Fig 1 (T_exe / f_out / T_comp) ===', flush=True)

t_zac = np.array([zac12.get(n,{}).get('T_exe',np.nan) for n in names])
t_pm = np.array([data[n]['base_texe'] for n in names])
t_os = np.array([data[n]['orbs_texe'] for n in names])
t_ob = np.array([data[n]['orbit_texe'] for n in names])
f_zac = np.array([zac12.get(n,{}).get('fidelity',np.nan) for n in names])
f_pm = np.array([data[n]['base_fid'] for n in names])
f_os = np.array([data[n]['orbs_fid'] for n in names])
f_ob = np.array([data[n]['orbit_fid'] for n in names])
cc_zac = np.array([zac12.get(n,{}).get('compilation_time',np.nan) for n in names])
cc_pm = np.array([data[n]['comp_base'] for n in names])
cc_os = np.array([data[n]['comp_orbs'] for n in names])
cc_ob = np.array([data[n]['comp_orbit'] for n in names])

nt_zac=t_zac/t_ob; nt_pm=t_pm/t_ob; nt_os=t_os/t_ob; nt_ob_=np.ones_like(t_ob)
nf_zac=np.ones_like(f_zac); nf_pm=f_pm/f_zac; nf_os=f_os/f_zac; nf_ob=f_ob/f_zac

print('ORBIT vs ZAC T_exe:  %.3fx' % gm(t_zac/t_ob))
print('ORBIT vs PM T_exe:   %.3fx' % gm(t_pm/t_ob))
print('ORBIT* vs PM T_exe:  %.3fx' % gm(t_pm/t_os))

x=np.arange(len(LABELS)); w=0.17
off=np.array([-(1.5*w),-(0.5*w),(0.5*w),(1.5*w)])
fig,(ax1,ax2,ax3)=plt.subplots(3,1,figsize=(16,13.5),gridspec_kw={'height_ratios':[1,1,1]})
FS=38
FS_TICK=32
FS_LEG=28
FS_XLAB=30
ax1.bar(x+off[0],nt_zac,w,label='ZAC',color=c_zac,edgecolor='white',lw=.5)
ax1.bar(x+off[1],nt_pm,w,label='PowerM',color=c_pm,edgecolor='white',lw=.5)
ax1.bar(x+off[2],nt_os,w,label='ORBIT*',color=c_os,edgecolor='white',lw=.5)
ax1.bar(x+off[3],nt_ob_,w,label='ORBIT',color=c_ob,edgecolor='white',lw=.5)
ax1.set_ylabel('(a) Norm. $T_{exe}$',fontsize=FS,labelpad=10)
ax1.set_yscale('log'); ax1.axhline(1,color='k',lw=.8,ls='--',alpha=.5)
ax1.tick_params(axis='y',labelsize=FS_TICK); ax1.tick_params(axis='x',labelbottom=False)
ax1.legend(loc='upper right',fontsize=FS_LEG,ncol=4,handletextpad=0.2,columnspacing=0.6)
ax2.bar(x+off[0],nf_zac,w,color=c_zac,edgecolor='white',lw=.5)
ax2.bar(x+off[1],nf_pm,w,color=c_pm,edgecolor='white',lw=.5)
ax2.bar(x+off[2],nf_os,w,color=c_os,edgecolor='white',lw=.5)
ax2.bar(x+off[3],nf_ob,w,color=c_ob,edgecolor='white',lw=.5)
ax2.set_ylabel('(b) Norm. $f_{out}$',fontsize=FS,labelpad=10)
ax2.set_yscale('log'); ax2.axhline(1,color='k',lw=.8,ls='--',alpha=.5)
ax2.set_yticks([1e0, 1e10, 1e20]); ax2.set_yticklabels(['$10^{0}$','$10^{10}$','$10^{20}$'])
ax2.tick_params(axis='y',labelsize=FS_TICK); ax2.tick_params(axis='x',labelbottom=False)
for i in range(len(x)):
    if not np.isnan(cc_zac[i]): ax3.bar(x[i]+off[0],cc_zac[i],w,color=c_zac,edgecolor='white',lw=.5)
    ax3.bar(x[i]+off[1],cc_pm[i],w,color=c_pm,edgecolor='white',lw=.5)
    ax3.bar(x[i]+off[2],cc_os[i],w,color=c_os,edgecolor='white',lw=.5)
    ax3.bar(x[i]+off[3],cc_ob[i],w,color=c_ob,edgecolor='white',lw=.5)
ax3.set_ylabel('(c) $T_{comp}$ (s)',fontsize=FS,labelpad=10)
ax3.set_yscale('log'); ax3.tick_params(axis='y',labelsize=FS_TICK)
ax3.set_xticks(x); ax3.set_xticklabels(LABELS,rotation=20,ha='right',fontsize=FS_XLAB,rotation_mode='anchor')
plt.subplots_adjust(hspace=0.10); plt.tight_layout()
plt.savefig(f'{OUT}/rebuttal_main_fig4_new_ablation_qft63_v3.png',dpi=200,bbox_inches='tight')
plt.savefig(f'{OUT}/rebuttal_main_fig4_new_ablation_qft63_v3.pdf',bbox_inches='tight')
plt.close()
print('Saved rebuttal_main_fig4_new_ablation_qft63_v3', flush=True)

# =========================================================================
# FIGURE 2: fig11_zone_arch_new — LINE PLOT (3 rows × 10 columns)
# Each column = one benchmark, x-axis = #Qubits, 2 lines: PowerM (--) / ORBIT (-)
# =========================================================================
print('\n=== Generating Fig 2 (Motion Ablation — line plots) ===', flush=True)

# Use tuned per-size data if available, else fall back to multisize
_fig11_path = 'local_mvqc/fig11_tuned_per_size.json' if os.path.exists('local_mvqc/fig11_tuned_per_size.json') else 'local_mvqc/fig11_multisize_data.json'
print(f'  Using fig11 data from: {_fig11_path}', flush=True)
with open(_fig11_path) as _f:
    fig11_data = json.load(_f)

import string
BENCH_FAMILIES = ['cat','ghz','wstate','ising','qft','qaoa_regu','qv','rcs']
BENCH_TITLES = ['cat','ghz','wstate','ising','qft','qaoa','qv','rcs']
COLORS_FIG11 = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
                '#8c564b','#e377c2','#17becf','#bcbd22','#7f7f7f']

n_col = len(BENCH_FAMILIES)
COL_W = 3.2
fig2, axs = plt.subplots(3, n_col, figsize=(COL_W * n_col, 9.5), squeeze=False, sharey='row')

TITLE_FS=32; LEGEND_FS=18; TICK_FS=22; LABEL_FS=28; LW=2.5; MS=6

for col, (fam, title) in enumerate(zip(BENCH_FAMILIES, BENCH_TITLES)):
    color = COLORS_FIG11[col]
    fd = fig11_data.get(fam, {'sizes':[],'base':[],'orbit':[]})
    xs = fd['sizes']

    for label, key, ls in [('PowerM','base','--'), ('ORBIT','orbit','-')]:
        chain_vals = [d['max_chain'] for d in fd[key]]
        avg_vals = [d['avg_move_aod'] for d in fd[key]]
        aod_vals = [d['aod_rounds'] for d in fd[key]]

        if chain_vals:
            axs[0, col].plot(xs, chain_vals, marker='o', linestyle=ls, color=color, linewidth=LW, markersize=MS, label=label)
        if avg_vals:
            axs[1, col].plot(xs, avg_vals, marker='o', linestyle=ls, color=color, linewidth=LW, markersize=MS, label=label)
        if aod_vals:
            axs[2, col].plot(xs, aod_vals, marker='o', linestyle=ls, color=color, linewidth=LW, markersize=MS)

    axs[0, col].set_title(title, fontsize=TITLE_FS, pad=10)
    axs[0, col].legend(loc='upper left', fontsize=LEGEND_FS, frameon=True)

    for row in range(3):
        axs[row, col].grid(axis='y', linestyle='--', linewidth=0.8, alpha=0.5)
        axs[row, col].tick_params(axis='y', labelsize=TICK_FS)
        axs[row, col].tick_params(axis='x', labelsize=TICK_FS-2)
        sublabel = f'({string.ascii_lowercase[row]}{col+1})'
        axs[row, col].text(0.5, -0.18, sublabel, transform=axs[row,col].transAxes,
                           fontsize=TICK_FS, ha='center', va='top', clip_on=False)
    axs[2, col].set_xlabel('#Qubits', fontsize=LABEL_FS-2, labelpad=30)

# Align y-axis labels vertically using fixed x coordinate
_ylabel_x = -0.30
axs[0,0].set_ylabel('Max Chain\nLength', fontsize=LABEL_FS, rotation=90)
axs[1,0].set_ylabel('Aver. #Move\nper AOD', fontsize=LABEL_FS, rotation=90)
axs[2,0].set_ylabel('AOD Rounds\n(Log)', fontsize=LABEL_FS, rotation=90)
axs[0,0].yaxis.set_label_coords(_ylabel_x, 0.5)
axs[1,0].yaxis.set_label_coords(_ylabel_x, 0.5)
axs[2,0].yaxis.set_label_coords(_ylabel_x, 0.5)
for col in range(n_col):
    axs[2, col].set_yscale('log')

plt.subplots_adjust(hspace=0.45, wspace=0.35)
plt.tight_layout()
plt.savefig(f'{OUT}/fig11_zone_arch_new.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{OUT}/fig11_zone_arch_new.pdf', bbox_inches='tight')
plt.close()
print('Saved fig11_zone_arch_new', flush=True)

# --- Also generate nolog version with identical layout ---
fig2b, axs2b = plt.subplots(3, n_col, figsize=(COL_W * n_col, 9.5), squeeze=False, sharey='row')
for col, (fam, title) in enumerate(zip(BENCH_FAMILIES, BENCH_TITLES)):
    color = COLORS_FIG11[col]
    fd = fig11_data.get(fam, {'sizes':[],'base':[],'orbit':[]})
    xs = fd['sizes']
    for label, key, ls in [('PowerM','base','--'), ('ORBIT','orbit','-')]:
        chain_vals = [d['max_chain'] for d in fd[key]]
        avg_vals = [d['avg_move_aod'] for d in fd[key]]
        aod_vals = [d['aod_rounds'] for d in fd[key]]
        if chain_vals:
            axs2b[0, col].plot(xs, chain_vals, marker='o', linestyle=ls, color=color, linewidth=LW, markersize=MS, label=label)
        if avg_vals:
            axs2b[1, col].plot(xs, avg_vals, marker='o', linestyle=ls, color=color, linewidth=LW, markersize=MS, label=label)
        if aod_vals:
            axs2b[2, col].plot(xs, aod_vals, marker='o', linestyle=ls, color=color, linewidth=LW, markersize=MS)
    axs2b[0, col].set_title(title, fontsize=TITLE_FS, pad=10)
    axs2b[0, col].legend(loc='upper left', fontsize=LEGEND_FS, frameon=True)
    for row in range(3):
        axs2b[row, col].grid(axis='y', linestyle='--', linewidth=0.8, alpha=0.5)
        axs2b[row, col].tick_params(axis='y', labelsize=TICK_FS)
        axs2b[row, col].tick_params(axis='x', labelsize=TICK_FS-2)
        sublabel = f'({string.ascii_lowercase[row]}{col+1})'
        axs2b[row, col].text(0.5, -0.18, sublabel, transform=axs2b[row,col].transAxes,
                             fontsize=TICK_FS, ha='center', va='top', clip_on=False)
    axs2b[2, col].set_xlabel('#Qubits', fontsize=LABEL_FS-2, labelpad=30)
    # NO log on row 3
axs2b[0,0].set_ylabel('Max Chain\nLength', fontsize=LABEL_FS, rotation=90)
axs2b[1,0].set_ylabel('Aver. #Move\nper AOD', fontsize=LABEL_FS, rotation=90)
axs2b[2,0].set_ylabel('AOD Rounds', fontsize=LABEL_FS, rotation=90)
axs2b[0,0].yaxis.set_label_coords(_ylabel_x, 0.5)
axs2b[1,0].yaxis.set_label_coords(_ylabel_x, 0.5)
axs2b[2,0].yaxis.set_label_coords(_ylabel_x, 0.5)
plt.subplots_adjust(hspace=0.45, wspace=0.35)
plt.tight_layout()
plt.savefig(f'{OUT}/fig11_zone_arch_new_nolog.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{OUT}/fig11_zone_arch_new_nolog.pdf', bbox_inches='tight')
plt.close()
print('Saved fig11_zone_arch_new_nolog', flush=True)

# =========================================================================
# FIGURE 3: fid_aba_new — inverted bars (matching fid_aba_v2.py style)
# =========================================================================
print('\n=== Generating Fig 3 (Fidelity Ablation) ===', flush=True)
from matplotlib.ticker import LogLocator, NullLocator, NullFormatter, LogFormatterMathtext

# ZAC fidelity breakdown from fig13
# For ZAC exic, use orbit's exic (circuit-dependent, same for all compilers)
zac_exic_arr = [data[n]['orbit_exic'] for n in names]  # same circuit = same 2q error
zac_trans_arr = [zac13.get(n,{}).get('fidelity_atom_transfer',np.nan) for n in names]
zac_coh_arr = [zac13.get(n,{}).get('fidelity_coherence',np.nan) for n in names]

pm_exic_arr = [data[n]['base_exic'] for n in names]
pm_trans_arr = [data[n]['base_trans'] for n in names]
pm_coh_arr = [data[n]['base_coh'] for n in names]
ob_exic_arr = [data[n]['orbit_exic'] for n in names]
ob_trans_arr = [data[n]['orbit_trans'] for n in names]
ob_coh_arr = [data[n]['orbit_coh'] for n in names]

datasets_aba = [
    [pm_exic_arr, ob_exic_arr],
    [pm_trans_arr, ob_trans_arr],
    [pm_coh_arr, ob_coh_arr],
]
titles_aba = ['(a) 2q error comparison', '(b) atom transfer error comparison', '(c) decoherence comparison']
colors_aba = [c_pm, c_ob]
cases_aba = ['PowerM', 'ORBIT']

fig3, axes3 = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
FS3 = 22; w3 = 0.18
offsets3 = np.array([-0.10, 0.10])
legend_handles = []

for i, ax in enumerate(axes3):
    ax.set_title(titles_aba[i], pad=6, fontsize=FS3)
    ax.set_yscale('log')
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,), numticks=10))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis='y', labelsize=FS3-2)
    subplot_pos_vals = []

    for j, case in enumerate(cases_aba):
        vals = np.array(datasets_aba[i][j], dtype=float)
        valid = np.isfinite(vals) & (vals > 0)
        subplot_pos_vals.extend(vals[valid])

        heights = np.where(valid, 1 - vals, 0)
        bottoms = np.where(valid, vals, 1)

        rects = ax.bar(x + offsets3[j], height=heights, bottom=bottoms, width=w3,
                       color=colors_aba[j], linewidth=0.6)
        if i == 0:
            legend_handles.append(rects)

    ax.grid(axis='y', linestyle='--', alpha=0.3)
    if subplot_pos_vals:
        ymin = min(subplot_pos_vals) / 100.0
        ymax = 1.05
        ax.set_ylim(bottom=max(ymin, 1e-40), top=ymax)
        import matplotlib.ticker as ticker
        exp_lo = int(np.floor(np.log10(max(ymin, 1e-40))))
        exp_hi = 0
        step = max(1, int(np.ceil((exp_hi - exp_lo) / 2)))
        tick_exps = list(range(exp_lo, exp_hi + 1, step))
        if 0 not in tick_exps:
            tick_exps.append(0)
        tick_exps = sorted(tick_exps)[:3]
        yticks = [10.0 ** e for e in tick_exps if max(ymin, 1e-40) <= 10.0 ** e <= ymax]
        if len(yticks) < 2:
            yticks = [10.0 ** exp_lo, 10.0 ** 0]
        ax.set_yticks(yticks)
        ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))

axes3[-1].set_xticks(x)
axes3[-1].set_xticklabels(LABELS, fontsize=FS3-2, rotation=20, ha='right', rotation_mode='anchor')
fig3.text(0.01, 0.5, 'Fidelity (Log Scale)', va='center', rotation='vertical', fontsize=FS3)
axes3[0].legend(legend_handles, cases_aba,
                loc='upper right', fontsize=FS3-2, frameon=True, edgecolor='black',
                ncol=1, columnspacing=0.5, handletextpad=0.3)

plt.subplots_adjust(left=0.12, right=0.97, top=0.95, bottom=0.16, hspace=0.28)
plt.savefig(f'{OUT}/fid_aba_new.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{OUT}/fid_aba_new.pdf', bbox_inches='tight')
plt.close()
print('Saved fid_aba_new', flush=True)

# Print full data table
print('\n' + '='*120, flush=True)
print('FULL DATA TABLE', flush=True)
print('='*120, flush=True)
print('%-18s %10s %10s %10s | %6s %6s %6s | %6s %6s %6s | %5s %5s %5s' %
      ('Bench','PM T','O* T','ORB T','PM ch','O* ch','O ch','PM am','O* am','O am','PM ar','O* ar','O ar'), flush=True)
print('-'*120, flush=True)
for n in names:
    d = data[n]
    print('%-18s %10.0f %10.0f %10.0f | %6d %6d %6d | %6.2f %6.2f %6.2f | %5d %5d %5d' %
          (n, d['base_texe'],d['orbs_texe'],d['orbit_texe'],
           d['base_max_chain'],d['orbs_max_chain'],d['orbit_max_chain'],
           d['base_avg_move_aod'],d['orbs_avg_move_aod'],d['orbit_avg_move_aod'],
           d['base_aod_rounds'],d['orbs_aod_rounds'],d['orbit_aod_rounds']), flush=True)

print('\n%-18s %12s %12s %12s | %12s %12s %12s | %12s %12s %12s' %
      ('Bench','PM exic','PM trans','PM coh','O exic','O trans','O coh','ZAC exic','ZAC trans','ZAC coh'), flush=True)
print('-'*140, flush=True)
for i,n in enumerate(names):
    d = data[n]
    print('%-18s %12.4e %12.4e %12.4e | %12.4e %12.4e %12.4e | %12.4e %12.4e %12.4e' %
          (n, d['base_exic'],d['base_trans'],d['base_coh'],
           d['orbit_exic'],d['orbit_trans'],d['orbit_coh'],
           zac_exic_arr[i],zac_trans_arr[i],zac_coh_arr[i]), flush=True)

print('\nAll done!', flush=True)
