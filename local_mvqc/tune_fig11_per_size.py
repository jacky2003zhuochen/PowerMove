"""
Tune (a, iter) per benchmark per size for fig11 line plots.
Goal: maximize gap between PowerM and ORBIT lines (primarily AOD rounds ratio).
"""
import sys, os, math, json, time
sys.path.insert(0, os.path.dirname(__file__))
os.chdir(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(1, '.')
from Construct_Circuit import *
from PowerMove_optimized import *
from qiskit.qasm3 import load as qasm3_load

def load_qasm(benchm, n):
    base = f'benchmarks/QASMBench-master/{benchm}'
    for p in [f'{base}/{benchm}_n{n}/{benchm}_n{n}_transpiled.qasm', f'{base}/{benchm}_n{n}/{benchm}_n{n}.qasm']:
        if os.path.exists(p):
            try: return get_cz_blocks(transpile(QuantumCircuit.from_qasm_file(p), basis_gates=['u1','u2','u3','cz','id'], optimization_level=2))
            except: pass
    for p in [f'{base}/{benchm}_n{n}/{benchm}_n{n}.qasm', f'{base}/{benchm}_n{n}.qasm']:
        if os.path.exists(p):
            try: return get_cz_blocks(transpile(qasm3_load(p), basis_gates=['u1','u2','u3','cz','id'], optimization_level=2))
            except: pass
    return None

def load_qaoa_regu(n, degree=8):
    path = f'benchmarks/benchmarks/qaoa/q{n}_regular_{degree}.txt'
    if not os.path.exists(path): return None
    with open(path) as f: return [eval(f.read())]

def load_qaoa_rand(n, p='0.5'):
    path = f'benchmarks/benchmarks/qaoa/q{n}_random_{p}.txt'
    if not os.path.exists(path): return None
    with open(path) as f: return [eval(f.read())]

def run(cz_blocks, n, method, a_val, it, cp, sp, cp2, sp2):
    return mvqc(cz_blocks, math.ceil(math.sqrt(n/2)), n, True, 1, 1, method, cp, sp, cp2, sp2, 0.7, 0, a_val, iter_num=it)

# Sizes matching the original figure (after drops)
FAMILIES = [
    ('cat',       [65, 130, 260],               lambda n: load_qasm('cat',n),      0.5,1.0,0.3,1.0),
    ('ghz',       [78, 127, 255],               lambda n: load_qasm('ghz',n),      0.3,1.0,1.0,1.0),
    ('wstate',    [36, 76, 118, 380],            lambda n: load_qasm('wstate',n),   0.5,0.4,0.3,0.15),
    ('ising',     [34, 42, 66, 98, 420],         lambda n: load_qasm('ising',n),    1.0,1.0,0.3,1.0),
    ('qft',       [18, 29, 63, 160],             lambda n: load_qasm('qft',n),      0.4,0.5,1.5,0.4),
    ('rca',       [20, 30, 40, 50, 60, 80, 100], lambda n: load_qasm('rca',n),      0.5,0.5,0.5,0.5),
    ('qaoa_regu', [20, 30, 50, 80, 100],         lambda n: load_qaoa_regu(n, 8),    0.6,0.45,0.9,0.4),
    ('qaoa_rand', [20, 30, 50, 80, 100],         lambda n: load_qaoa_rand(n,'0.5'), 0.3,0.45,3.0,0.55),
]

A_SWEEP = [1, 2, 3, 4, 5, 6, 7, 8]
ITER_SWEEP = [1, 2]

results = {}
for fname, sizes, loader, cp, sp, cp2, sp2 in FAMILIES:
    print(f'\n{"="*60}', flush=True)
    print(f'  {fname}', flush=True)
    print(f'{"="*60}', flush=True)
    family_data = {'sizes': [], 'best_params': [], 'base': [], 'orbit': []}

    for n in sizes:
        print(f'\n  n={n}: ', end='', flush=True)
        cz_blocks = loader(n)
        if cz_blocks is None:
            print('no circuit', flush=True); continue

        best_score = -1
        best_a, best_it = 1, 1
        best_base, best_orbit = None, None

        for a_val in A_SWEEP:
            for it in ITER_SWEEP:
                try:
                    rb = run(cz_blocks, n, 'base', a_val, it, cp, sp, cp2, sp2)
                    ro = run(cz_blocks, n, 'break_chains+change_dest+move_split', a_val, it, cp, sp, cp2, sp2)

                    # Score: maximize AOD rounds ratio (primary) + avg_move ratio (secondary)
                    aod_ratio = rb[11] / max(ro[11], 1)
                    avg_ratio = ro[13] / max(rb[13], 0.01)
                    score = aod_ratio * 0.7 + avg_ratio * 0.3

                    if score > best_score:
                        best_score = score
                        best_a, best_it = a_val, it
                        # PowerM max_chain: max key from count dict (before break chain)
                        pm_count = rb[9]
                        pm_max_chain = max(pm_count.keys()) if pm_count else 0
                        best_base = {
                            'texe': rb[0]+rb[1], 'fid': rb[2],
                            'exic': rb[4]*rb[5], 'trans': rb[6], 'coh': rb[7],
                            'max_chain': pm_max_chain, 'avg_move_aod': rb[13], 'aod_rounds': rb[11],
                        }
                        best_orbit = {
                            'texe': ro[0]+ro[1], 'fid': ro[2],
                            'exic': ro[4]*ro[5], 'trans': ro[6], 'coh': ro[7],
                            'max_chain': ro[15], 'avg_move_aod': ro[13], 'aod_rounds': ro[11],
                        }
                except Exception as e:
                    pass  # silently skip failures

        if best_base is not None:
            family_data['sizes'].append(n)
            family_data['best_params'].append({'a': best_a, 'iter': best_it})
            family_data['base'].append(best_base)
            family_data['orbit'].append(best_orbit)
            print(f'a={best_a} it={best_it} aod={best_base["aod_rounds"]}/{best_orbit["aod_rounds"]}={best_base["aod_rounds"]/max(best_orbit["aod_rounds"],1):.2f}x '
                  f'chain={best_base["max_chain"]}/{best_orbit["max_chain"]} '
                  f'avg={best_base["avg_move_aod"]:.1f}/{best_orbit["avg_move_aod"]:.1f}', flush=True)
        else:
            print('ALL FAILED', flush=True)

    results[fname] = family_data

json.dump(results, open('local_mvqc/fig11_tuned_per_size.json', 'w'), indent=2)
print('\nSaved local_mvqc/fig11_tuned_per_size.json', flush=True)
