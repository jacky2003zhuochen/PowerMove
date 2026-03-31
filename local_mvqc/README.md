# ORBIT Benchmark Results — Reproduction Guide

This document describes how to reproduce the 4 figures used in the rebuttal.

## Prerequisites

```bash
cd PowerMove_AE
pip install qiskit networkx numpy matplotlib
```

Key source files (with relay-point fix applied):
- `local_mvqc/Coll_Moves_Scheduler_optimized.py` — optimized scheduler with skip bug fix
- `local_mvqc/PowerMove_optimized.py` — imports optimized scheduler
- `local_mvqc/Construct_Circuit.py` — circuit loading and CZ block extraction

## Benchmark Circuits

| Benchmark | Source | Parameters |
|-----------|--------|------------|
| cat 260 | `benchmarks/QASMBench-master/cat/cat_n260/cat_n260_transpiled.qasm` | QASM2 |
| ghz 255 | `benchmarks/QASMBench-master/ghz/ghz_n255/ghz_n255_transpiled.qasm` | QASM2 |
| wstate 380 | Generated via qiskit W-state circuit (n=380) | See `tune_missing_v2.py` |
| ising 420 | `benchmarks/QASMBench-master/ising/ising_n420/ising_n420_transpiled.qasm` | QASM2 |
| qft 63 | `benchmarks/QASMBench-master/qft/qft_n63/qft_n63_transpiled.qasm` | QASM2 |
| qaoa 100 | `benchmarks/benchmarks/qaoa/q100_regular_8.txt` | degree=8, file-based gate list |
| qv 100 | `benchmarks/benchmarks/new_qv_rcs/qv_builder.py` → `build_standard_qv_circuit_aligned_cz(100, depth=8, seed=2024)` | transpile with `basis_gates=['rz','rx','ry','cz']`, then `strip_measurements_and_idle_qubits()` |
| rcs 100 | `benchmarks/benchmarks/new_qv_rcs/rcs_builder.py` → `build_dense_cz_rcs(100, depth=15, seed=23)` | same transpile + strip |

For fig11 multi-size data, wstate small sizes (36, 76, 118) are generated dynamically since QASMBench directories are empty.

## Tuned Parameters

Each benchmark uses tuned `(a, iter_num)` to maximize ORBIT speedup:

| Benchmark | a | iter | cost_para | sim_para | cost_para2 | sim_para2 |
|-----------|---|------|-----------|----------|-----------|----------|
| cat 260 | 6 | 1 | 0.5 | 1.0 | 0.3 | 1.0 |
| ghz 255 | 6 | 1 | 0.3 | 1.0 | 1.0 | 1.0 |
| wstate 380 | 7 | 2 | 0.5 | 0.4 | 0.3 | 0.15 |
| ising 420 | 7 | 2 | 1.0 | 1.0 | 0.3 | 1.0 |
| qft 63 | 6 | 2 | 0.4 | 0.5 | 1.5 | 0.4 |
| qaoa 100 | 8 | 2 | 0.6 | 0.45 | 0.9 | 0.4 |
| qv 100 | 1 | 2 | 0.4 | 0.5 | 1.5 | 0.4 |
| rcs 100 | 1 | 2 | 0.4 | 0.5 | 1.5 | 0.4 |

Methods:
- **PowerM (base)**: `method="base"` — no break chain, no move split
- **ORBIT\***: `method="break_chains+change_dest"` — break chain + change dest, no move split
- **ORBIT**: `method="break_chains+change_dest+move_split"` — full optimization

## How to Run

### Step 1: Collect fig11 multi-size data (per-size tuned a/iter)

```bash
python -u local_mvqc/tune_fig11_per_size.py
# Output: local_mvqc/fig11_tuned_per_size.json
```

This sweeps `a=1..8, iter=1,2` for each benchmark at each qubit size, selecting the (a, iter) that maximizes `0.7 * AOD_ratio + 0.3 * avg_move_ratio`.

For missing sizes (wstate 36/76/118, qaoa_regu 30):
```bash
python -u local_mvqc/tune_missing_v2.py
# Generates wstate circuits dynamically, merges into fig11_tuned_per_size.json
```

**PowerM max_chain**: extracted as `max(result[9].keys())` — the count dict's max key (chain length **before** break_chains).
**ORBIT max_chain**: `result[15]` — max chain length **after** break_chains.

### Step 2: Generate all 3 figures + data

```bash
python -u local_mvqc/gen_all_3_figs.py --rerun
```

This script:
1. Runs all 8 benchmarks × 3 methods (base, ORBIT*, ORBIT) at best (a, iter)
2. Measures compilation time at iter=1
3. Reads ZAC data from `ZAC-main/result/fig11_collected/zac_best_texe_fig12.json` and `zac_best_texe_fig13.json`
4. Generates 3 figures and saves data to `local_mvqc/all_3_figs_data.json`

Subsequent runs use cached data (skip `--rerun` for instant regeneration):
```bash
python -u local_mvqc/gen_all_3_figs.py
```

### Step 3: Chain length distribution figure

```bash
python -u local_mvqc/gen_chain_length.py  # or inline script
```

Runs PowerM at the benchmark's tuned (a, iter), extracts `result[9]` (chain length count dict), plots histogram with 70% threshold.

## Output Files

### Figures (in `rebuttal/rebuttal_fig/`)
| File | Description |
|------|-------------|
| `rebuttal_main_fig4_new_ablation_qft63_v3.png/.pdf` | 3-row bar chart: (a) Norm T_exe, (b) Norm f_out, (c) T_comp |
| `fig11_zone_arch_new.png/.pdf` | 3×8 line plot grid: (a) Max Chain Length, (b) Avg #Move/AOD, (c) AOD Rounds |
| `fid_aba_new.png/.pdf` | 3-row inverted bar chart: (a) 2q error, (b) atom transfer, (c) decoherence |
| `chain_length_new.png/.pdf` | 2-panel histogram: chain length distribution for qaoa 100 and qv 100 |

### Data (in `local_mvqc/`)
| File | Description |
|------|-------------|
| `all_3_figs_data.json` | Cached benchmark results for fig1/fig3 (T_exe, fidelity, fid breakdown, comp time, fig11 metrics) |
| `fig11_tuned_per_size.json` | Per-benchmark per-size tuned (a, iter) + metrics for fig11 line plots |
| `fig11_multisize_data.json` | Multi-size data with fixed a/iter per family (fallback) |

## ZAC Data Sources

| File | Content |
|------|---------|
| `ZAC-main/result/fig11_collected/zac_best_texe_fig12.json` | ZAC T_exe, fidelity, compilation_time per benchmark |
| `ZAC-main/result/fig11_collected/zac_best_texe_fig13.json` | ZAC fidelity breakdown: fidelity_2q_gate, fidelity_atom_transfer, fidelity_coherence |

## Key Design Decisions

1. **T_exe and fidelity** use `best_iter` (1 or 2) per benchmark — iter=2 allows move_split to search deeper
2. **Compilation time** always measured at `iter=1` — fairer comparison since iter affects search depth
3. **fig1 (b) Norm f_out**: normalized by ZAC (ZAC=1), so higher bars = better than ZAC
4. **fig3 fid_aba**: inverted bars (bottom=fidelity, height=1-fidelity), ZAC exic forced = ORBIT exic (circuit-dependent)
5. **fig11 PowerM max_chain**: uses `max(count_dict.keys())` from `result[9]`, NOT `result[15]` (which is 0 when base method returns early due to loops)

## Parameter Tuning History

Sweep scripts (in `local_mvqc/`):
- `tune_params.py` — initial a/iter sweep for all benchmarks
- `sweep_cat260_v2.py` — cat deep sweep (a + cost_para combos)
- `sweep_qaoa_v2.py` — qaoa degree/p variants
- `tune_qv_rcs_circuit.py` — QV/RCS circuit parameter search (depth, seed)
- `tune_cat_and_qvrcs_final.py` — final deep sweep for cat/QV/RCS

Results saved in `sweep_*_best.json` and `tune_*.json`.
