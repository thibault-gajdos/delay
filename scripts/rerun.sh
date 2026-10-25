#!/bin/bash
# rerun.sh -- every analysis downstream of the simulations, in order.
#
# The simulations themselves are not here: run_sweep.sh and run_leak.sh produce them
# (README, "Order of execution"). Everything below is deterministic given their output.
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

step "concatenated sweep and its summary"
python3 sweep_summary.py > /dev/null && echo "wrote results/sweep_screened.jsonl, results/sweep_summary.json"

step "model-data comparison (classification of every admissible draw)"
python3 classify.py

step "leak grid"
python3 leak_summary.py

step "effect panel (ranges and signs over the admissible draws)"
python3 effect_panel.py results/sweep_screened.jsonl | tail -3

step "representative parameter vector (the point Figure 2B is drawn at)"
python3 representative.py > /dev/null && echo "wrote results/representative.json"

step "the CSV bridge for R (figures_space.py reads it too)"
python3 export_for_R.py

step "figures"
python3 figures.py
python3 figures_space.py results/sweep_screened.jsonl
python3 mkfig.py

step "statistics on the data"
( cd R && Rscript run_stats.R > /dev/null && echo "wrote R/stats_tables.txt and R/stats_results.csv" )

step "done"
