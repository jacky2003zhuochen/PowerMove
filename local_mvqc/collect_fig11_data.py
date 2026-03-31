"""
Collect data for fig11 line plots: run multiple qubit sizes per benchmark.
Uses same (a, cost_params) as gen_all_3_figs.py for each benchmark family.
"""
import sys, os, math, json, time
sys.path.insert(0, os.path.dirname(__file__))
os.chdir(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(1, '.')
from Construct_Circuit import *
from PowerMove_optimized import *
from qiskit.qasm3 import load as qasm3_load
sys.path.insert(0, 'benchmarks/benchmarks/new_qv_rcs')
from qv_builder import build_standard_qv_circuit_aligned_cz
from rcs_builder import build_dense_cz_rcs, strip_measurements_and_idle_qubits

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

def load_qv(n, depth=8, seed=2024):
    circ = build_standard_qv_circuit_aligned_cz(n, depth=depth, seed=seed, add_measurements=True, insert_barriers=True)
    return get_cz_blocks(strip_measurements_and_idle_qubits(transpile(circ, basis_gates=['rz','rx','ry','cz'], optimization_level=2, seed_transpiler=seed)))

def load_rcs(n, depth=15, seed=23):
    circ = build_dense_cz_rcs(n, depth=depth, seed=seed, add_measurements=True, matching='random', insert_barriers=True)
    return get_cz_blocks(strip_measurements_and_idle_qubits(transpile(circ, basis_gates=['rz','rx','ry','cz'], optimization_level=2, seed_transpiler=seed)))

def run(cz_blocks, n, method, a_val, it, cp, sp, cp2, sp2):
    return mvqc(cz_blocks, math.ceil(math.sqrt(n/2)), n, True, 1, 1, method, cp, sp, cp2, sp2, 0.7, 0, a_val, iter_num=it)

# Benchmark families: (name, sizes, loader_func, a, best_iter, cp, sp, cp2, sp2)
FAMILIES = [
    ('cat',       [22, 35, 65, 130, 260],        lambda n: load_qasm('cat',n),      6, 1, 0.5,1.0,0.3,1.0),
    ('ghz',       [23, 40, 78, 127, 255],         lambda n: load_qasm('ghz',n),      6, 1, 0.3,1.0,1.0,1.0),
    ('wstate',    [27, 36, 76, 118, 380],          lambda n: load_qasm('wstate',n),   7, 2, 0.5,0.4,0.3,0.15),
    ('ising',     [26, 34, 42, 66, 98, 420],       lambda n: load_qasm('ising',n),    7, 2, 1.0,1.0,0.3,1.0),
    ('qft',       [18, 29, 63],                    lambda n: load_qasm('qft',n),      6, 2, 0.4,0.5,1.5,0.4),
    ('rca',       [10, 20, 30, 40, 50, 60, 80],    lambda n: load_qasm('rca',n),      1, 2, 0.5,0.5,0.5,0.5),
    ('qaoa_regu', [20, 40, 50, 60, 80, 100],       lambda n: load_qaoa_regu(n, 8),    8, 2, 0.6,0.45,0.9,0.4),
    ('qaoa_rand', [10, 20, 30, 50, 80, 100],       lambda n: load_qaoa_rand(n,'0.5'), 1, 2, 0.3,0.45,3.0,0.55),
    ('qv',        [20, 50, 80, 100],                lambda n: load_qv(n, depth=8, seed=2024),  1, 2, 0.4,0.5,1.5,0.4),
    ('rcs',       [20, 50, 80, 100],                lambda n: load_rcs(n, depth=15, seed=23),  1, 2, 0.4,0.5,1.5,0.4),
]

results = {}
for fname, sizes, loader, a, it, cp, sp, cp2, sp2 in FAMILIES:
    print(f'\n=== {fname} (a={a}, iter={it}) ===', flush=True)
    family_data = {'sizes': [], 'base': [], 'orbit': []}
    for n in sizes:
        print(f'  n={n}...', end='', flush=True)
        try:
            cz_blocks = loader(n)
            if cz_blocks is None:
                print(' no circuit', flush=True); continue
            rb = run(cz_blocks, n, 'base', a, it, cp, sp, cp2, sp2)
            ro = run(cz_blocks, n, 'break_chains+change_dest+move_split', a, it, cp, sp, cp2, sp2)
            family_data['sizes'].append(n)
            family_data['base'].append({
                'texe': rb[0]+rb[1], 'fid': rb[2],
                'exic': rb[4]*rb[5], 'trans': rb[6], 'coh': rb[7],
                'max_chain': max(rb[9].keys()) if rb[9] else 0, 'avg_move_aod': rb[13], 'aod_rounds': rb[11],
            })
            family_data['orbit'].append({
                'texe': ro[0]+ro[1], 'fid': ro[2],
                'exic': ro[4]*ro[5], 'trans': ro[6], 'coh': ro[7],
                'max_chain': ro[15], 'avg_move_aod': ro[13], 'aod_rounds': ro[11],
            })
            print(f' chain={rb[15]}/{ro[15]} aod={rb[11]}/{ro[11]} avg={rb[13]:.1f}/{ro[13]:.1f}', flush=True)
        except Exception as e:
            print(f' FAILED: {str(e)[:60]}', flush=True)
    results[fname] = family_data

json.dump(results, open('local_mvqc/fig11_multisize_data.json', 'w'), indent=2)
print('\nSaved local_mvqc/fig11_multisize_data.json', flush=True)
