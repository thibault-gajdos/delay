# Scripts

Code for the simulations, the model–data comparison, the figures and the statistics of
the manuscript (`../manuscript/latex/delay.tex`). Every script here is used in the paper.

## Setup

Python, through uv:

```sh
cd scripts
uv sync
export DELAY_DATA="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"   # unpacked OSF archive
```

R, through renv, in `scripts/R` (packages: afex, emmeans, dplyr, tidyr):

```sh
cd scripts/R
Rscript -e 'renv::restore()'      # or renv::init() and renv::install(...) the first time
```

`rerun.sh` sets `ARCHIVE`, the archive path read by `run_stats.R`, to `DELAY_DATA`.

## Order of execution

1. `./run_sweep.sh [P]` — the main sweep: 1200 parameter vectors, 40000 trials per
   condition, both models, no leak. Writes `results/shards/`. Can be interrupted and
   relaunched, with any shard count: completed vectors are skipped.
2. `./run_leak.sh [P]` — the leak grid: γ_x, γ_q ∈ {0, 0.25, 0.5, 1} s⁻¹, the 15
   combinations other than (0, 0), on the vectors admissible in either model without leak.
   Writes `results/leak/`. Same resume rule.
3. `uv run ./rerun.sh` — everything downstream, deterministic given steps 1 and 2: the
   comparison, the tables, the figures and the statistics on the data.

Then `../manuscript/latex/make.sh` builds the PDF.

## Files

Model and data (`delaymodel/`):

| file | role |
|---|---|
| `sim.py` | simulation of both models (Methods, "Theoretical framework") |
| `empirical.py` | per-trial data from the OSF archive: coherence, duration, accuracy, RT, signal and response times, onset of the response-related activation, PMA times |
| `targets.py` | per-participant observed measures: accuracy, RT, PMT, MT, PMA-rate profiles in the stimulus-onset and delay windows, their slopes and curvatures |

Simulations:

| file | role |
|---|---|
| `fused_sweep.py` | one parameter vector: simulates both models, applies the admissibility criteria, writes the simulated measures and their Monte-Carlo standard errors (one JSON line per vector) |
| `run_sweep.sh`, `run_leak.sh` | parallel launch of `fused_sweep.py` |
| `shardio.py`, `repair_shards.py` | reading shard files; removing an incomplete last line after a crash |
| `sweep_summary.py` | concatenates the shards into `results/sweep_screened.jsonl` |
| `leak_index.py` | the list of vectors re-run with leak (`results/leak/adm_idx.json`) |

Comparison (Methods, "Model–data comparison"; Table 2; Figure 6):

| file | role |
|---|---|
| `effects.py` | definition of every effect compared, observed and simulated, with standard errors |
| `classify.py` | classification of every admissible vector: reproduces / does not reproduce / uncertain; sign when not fixed by derivation. Writes `results/classification.json` |
| `leak_summary.py` | the same verdicts at every leak combination. Writes `results/leak_summary.json` |
| `table_comparison.py` | Table 2 (`../manuscript/latex/table_comparison.tex`) |
| `fig_effects.py` | Figure 6 (`fig08.png`) |

Figures:

| file | figure |
|---|---|
| `representative.py` | chooses the vector Figure 2B is drawn at (`results/representative.json`) |
| `figures.py` | Figure 1 (`fig01.png`, `fig02.png`) and Figure 2B (`fig04.png`) |
| `fig_behavior.py` | Figure 3 (`fig05.png`) |
| `fig_delay_profiles.py` | Figure 5 (`fig07.png`) |
| `fig_approx_error.py` | Figure 7 (`fig09.png`) |

Figures 2A (`fig03.png`), 4 (`fig06.png`) and S1 (`fig10.png`) are not generated here.

Statistics on the data (Methods, "Statistical analysis"; Supplementary Tables):

| file | role |
|---|---|
| `R/run_stats.R` | every ANOVA, from the archive's export files. Writes `R/stats_tables.txt` and `R/stats_results.csv` |
| `stats_tables_tex.py` | the supplementary tables (`table_behav.tex`, `table_pma.tex`, `table_mbl.tex`, `table_lrp.tex`) |

`rerun.sh` runs steps 3 in this order: concatenation, comparison, leak summary, Table 2,
representative vector, figures, statistics.
