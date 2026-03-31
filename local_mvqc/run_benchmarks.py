"""
Run 9 zone-architecture benchmarks from rebuttal_main_fig4_data.json.
Methods: "base" (PowerMove) and "break_chains+change_dest+move_split" (ORBIT).
Uses the optimized versions of Coll_Moves_Scheduler and PowerMove.

Usage:
    cd PowerMove_AE
    python local_mvqc/run_benchmarks.py
"""

import sys, os, math, json, time

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
os.chdir(project_root)
sys.path.insert(0, os.path.join(project_root, 'local_mvqc'))
sys.path.insert(1, project_root)

from Construct_Circuit import *
from PowerMove_optimized import *

# QV/RCS builders
sys.path.insert(0, os.path.join(project_root, 'benchmarks', 'new_qv_rcs'))
from qv_builder import build_dense_cz_qv_like_circuit
from rcs_builder import build_dense_cz_rcs

# ---------------------------------------------------------------------------
# Circuit loaders (from scripts/zone_arch_fig11_data.py)
# ---------------------------------------------------------------------------

def load_qasm_circuit(benchm, n):
    try:
        circ = QuantumCircuit.from_qasm_file(
            f"benchmarks/QASMBench-master/{benchm}/{benchm}_n{n}/"
            f"{benchm}_n{n}_transpiled.qasm")
    except Exception:
        circ = QuantumCircuit.from_qasm_file(
            f"benchmarks/QASMBench-master/{benchm}/{benchm}_n{n}/"
            f"{benchm}_n{n}.qasm")
    tc = transpile(circ, basis_gates=["u1", "u2", "u3", "cz", "id"],
                   optimization_level=2)
    return get_cz_blocks(tc)

def load_qv_circuit(n):
    circ = build_dense_cz_qv_like_circuit(n, depth=n, seed=42,
                                          add_measurements=False,
                                          insert_barriers=False)
    tc = transpile(circ, basis_gates=["u1", "u2", "u3", "cz", "id"],
                   optimization_level=2)
    return get_cz_blocks(tc)

def load_rcs_circuit(n):
    circ = build_dense_cz_rcs(n, depth=n, seed=42,
                              add_measurements=False,
                              insert_barriers=False)
    tc = transpile(circ, basis_gates=["u1", "u2", "u3", "cz", "id"],
                   optimization_level=2)
    return get_cz_blocks(tc)

def load_qaoa_rand_circuit(n, P=0.5):
    path = f"benchmarks/qaoa/rand/q{n}_p{P}/i0.txt"
    with open(path, "r") as f:
        gates = eval(f.read())
    return [gates]

def load_qaoa_regu_circuit(n, degree=10, index=0):
    path = f"benchmarks/qaoa/regular/q{n}_regular{degree}/i{index}.txt"
    with open(path, "r") as f:
        gates = eval(f.read())
    return [gates]

# ---------------------------------------------------------------------------
# Benchmark configuration
# ---------------------------------------------------------------------------

def run_mvqc_storage(cz_blocks, n, method, cost_para, sim_para,
                     cost_para2, sim_para2, a_val, iter_num=2, d=1, h=0):
    Row = math.ceil(math.sqrt(n / 2))
    return mvqc(cz_blocks, Row, n, True, d, 1, method,
                cost_para, sim_para, cost_para2, sim_para2, 0.7, h, a_val,
                iter_num=iter_num)

# Benchmark name -> (loader_function, qubit_count)
BENCH_LOADERS = {
    'cat_260':       lambda: load_qasm_circuit('cat', 260),
    'ghz_255':       lambda: load_qasm_circuit('ghz', 255),
    'wstate_380':    lambda: load_qasm_circuit('wstate', 380),
    'ising_420':     lambda: load_qasm_circuit('ising', 420),
    'qft_160':       lambda: load_qasm_circuit('qft', 160),
    'qaoa_rand_100': lambda: load_qaoa_rand_circuit(100, P=0.5),
    'qaoa_regu_100': lambda: load_qaoa_regu_circuit(100, degree=10, index=0),
    'qv_100':        lambda: load_qv_circuit(100),
    'rcs_100':       lambda: load_rcs_circuit(100),
}

METHOD_BASE  = "base"
METHOD_ORBIT = "break_chains+change_dest+move_split"

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Load parameters
    with open(os.path.join(project_root, 'local_mvqc', 'rebuttal_main_fig4_data.json'), 'r') as f:
        fig4_data = json.load(f)

    params = fig4_data['parameters']
    ref_texe = fig4_data['T_exe']
    ref_fid = fig4_data['fidelity']
    bench_list = fig4_data['benchmarks']

    results = {
        'benchmarks': bench_list,
        'T_exe': {'PowerM': [], 'ORBIT': []},
        'fidelity': {'PowerM': [], 'ORBIT': []},
        'runtime_sec': {'PowerM': [], 'ORBIT': []},
    }

    for idx, bench in enumerate(bench_list):
        p = params[bench]
        n = p['n']
        a_val = p['a']
        cp, sp = p['cost_para'], p['sim_para']
        cp2, sp2 = p['cost_para2'], p['sim_para2']
        it = p['iter_num']

        print(f"\n{'='*60}")
        print(f"Benchmark: {bench} (n={n}, Row={p['Row']}, a={a_val})")
        print(f"{'='*60}")

        # Load circuit
        print(f"  Loading circuit...")
        t0 = time.time()
        cz_blocks = BENCH_LOADERS[bench]()
        load_time = time.time() - t0
        print(f"  Circuit loaded in {load_time:.1f}s")

        # Run PowerMove (base)
        print(f"\n  Running PowerMove (base)...")
        t0 = time.time()
        try:
            result_base = run_mvqc_storage(cz_blocks, n, METHOD_BASE,
                                           cp, sp, cp2, sp2, a_val, iter_num=it)
            base_texe = result_base[0] + result_base[1]
            base_fid = result_base[2]
            base_time = time.time() - t0
            print(f"    T_exe = {base_texe:.2f} us  (ref: {ref_texe['PowerM'][idx]})")
            print(f"    fidelity = {base_fid:.6e}  (ref: {ref_fid['PowerM'][idx]})")
            print(f"    runtime = {base_time:.1f}s")
        except Exception as e:
            print(f"    FAILED: {e}")
            import traceback; traceback.print_exc()
            base_texe, base_fid, base_time = None, None, None

        results['T_exe']['PowerM'].append(base_texe)
        results['fidelity']['PowerM'].append(base_fid)
        results['runtime_sec']['PowerM'].append(base_time)

        # Run ORBIT
        print(f"\n  Running ORBIT (break_chains+change_dest+move_split)...")
        t0 = time.time()
        try:
            result_orbit = run_mvqc_storage(cz_blocks, n, METHOD_ORBIT,
                                            cp, sp, cp2, sp2, a_val, iter_num=it)
            orbit_texe = result_orbit[0] + result_orbit[1]
            orbit_fid = result_orbit[2]
            orbit_time = time.time() - t0
            print(f"    T_exe = {orbit_texe:.2f} us  (ref: {ref_texe['ORBIT'][idx]})")
            print(f"    fidelity = {orbit_fid:.6e}  (ref: {ref_fid['ORBIT'][idx]})")
            print(f"    runtime = {orbit_time:.1f}s")
        except Exception as e:
            print(f"    FAILED: {e}")
            import traceback; traceback.print_exc()
            orbit_texe, orbit_fid, orbit_time = None, None, None

        results['T_exe']['ORBIT'].append(orbit_texe)
        results['fidelity']['ORBIT'].append(orbit_fid)
        results['runtime_sec']['ORBIT'].append(orbit_time)

    # Print summary table
    print(f"\n\n{'='*100}")
    print(f"SUMMARY")
    print(f"{'='*100}")
    print(f"{'Benchmark':<18} {'PowerM T_exe':>14} {'(ref)':>14} {'ORBIT T_exe':>14} {'(ref)':>14} {'Speedup':>8}")
    print(f"{'-'*18} {'-'*14} {'-'*14} {'-'*14} {'-'*14} {'-'*8}")
    for i, bench in enumerate(bench_list):
        pt = results['T_exe']['PowerM'][i]
        ot = results['T_exe']['ORBIT'][i]
        rp = ref_texe['PowerM'][i]
        ro = ref_texe['ORBIT'][i]
        pt_s = f"{pt:.0f}" if pt is not None else "FAIL"
        ot_s = f"{ot:.0f}" if ot is not None else "FAIL"
        sp = f"{pt/ot:.2f}x" if pt is not None and ot is not None and ot > 0 else "N/A"
        print(f"{bench:<18} {pt_s:>14} {rp:>14.0f} {ot_s:>14} {ro:>14.0f} {sp:>8}")

    print(f"\n{'Benchmark':<18} {'PowerM fid':>14} {'(ref)':>14} {'ORBIT fid':>14} {'(ref)':>14}")
    print(f"{'-'*18} {'-'*14} {'-'*14} {'-'*14} {'-'*14}")
    for i, bench in enumerate(bench_list):
        pf = results['fidelity']['PowerM'][i]
        of = results['fidelity']['ORBIT'][i]
        rpf = ref_fid['PowerM'][i]
        rof = ref_fid['ORBIT'][i]
        pf_s = f"{pf:.4e}" if pf is not None else "FAIL"
        of_s = f"{of:.4e}" if of is not None else "FAIL"
        print(f"{bench:<18} {pf_s:>14} {rpf:>14.4e} {of_s:>14} {rof:>14.4e}")

    # Save results
    output_path = os.path.join(project_root, 'local_mvqc', 'benchmark_results.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_path}")

if __name__ == '__main__':
    main()
