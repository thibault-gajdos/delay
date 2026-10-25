#!/bin/bash
# rerun.sh -- every analysis downstream of the simulations, in order.
#
# The simulations themselves are not here: run_sweep.sh and run_leak.sh produce them
# (README). Everything below is deterministic given their output.
#
# Usage, from scripts/, on the machine that holds the data:
#   export DELAY_DATA="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"
#   uv run ./rerun.sh
set -eu
cd "$(dirname "$0")"

: "${DELAY_DATA:?set DELAY_DATA to the unpacked OSF archive}"
: "${ARCHIVE:=$DELAY_DATA}"
export DELAY_DATA ARCHIVE

step () { printf '\n=== %s\n' "$1"; }

step "concatenated sweep (results/sweep_screened.jsonl)"
python3 sweep_summary.py

step "model-data comparison (classification of every admissible parameter vector)"
python3 classify.py

step "leak grid"
python3 leak_summary.py

step "comparison table (manuscript/latex/table_comparison.tex)"
python3 table_comparison.py

step "representative parameter vector (the point Figure 2B is drawn at)"
python3 representative.py > /dev/null && echo "wrote results/representative.json"

step "figures (manuscript/latex/figures/fig01,02,04,05,07,08,09.png)"
python3 figures.py
python3 fig_behavior.py
python3 fig_delay_profiles.py
python3 fig_approx_error.py
python3 fig_effects.py

step "statistics on the data"
( cd R && Rscript run_stats.R > /dev/null && echo "wrote R/stats_tables.txt and R/stats_results.csv" )
python3 stats_tables_tex.py

step "done"
