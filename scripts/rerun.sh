#!/bin/bash
# rerun.sh -- everything downstream of the sweep, in order, from a clean state.
#
# Not here, because each is long and is run sharded across cores by hand (see README,
# "Order of execution"): the sweep itself (run_sweep.sh), the window summaries the effect
# panel needs (augment_windows.py), the sharded re-simulation that effect_sizes.py needs
# for pool_noise.py, ramp_check.py, and the leak reruns (run_leak.sh). Everything else is deterministic given their outputs and takes about
# two minutes in total.
#
# Usage:  ./rerun.sh            (from scripts/)
# Requires DELAY_DATA and, for the R step, ARCHIVE.
set -e
cd "$(dirname "$0")"

: "${DELAY_DATA:?set DELAY_DATA to the unpacked OSF archive}"
: "${ARCHIVE:=$DELAY_DATA}"
export DELAY_DATA ARCHIVE

step () { printf '\n=== %s\n' "$1"; }

step "effect sizes, stage 1 (distributions and criterion counts)"
uv run python3 effect_sizes.py results/sweep_screened.jsonl 4 4000 0 > /dev/null

step "effect panel (the eight critical effects, against the data)"
uv run python3 effect_panel.py results/sweep_screened.jsonl | tail -12

step "representative parameter vector (the point the illustrative figures use)"
uv run python3 representative.py | grep -E '"draw"|"beta"' || true

step "the CSV bridge for R (figures_space.py reads pma_counts.csv, so this comes first)"
uv run python3 export_for_R.py

step "figures"
uv run python3 figures.py
uv run python3 figures_space.py results/sweep_screened.jsonl
uv run python3 mkfig.py

step "statistics"
( cd R && Rscript run_stats.R > /dev/null && echo "wrote R/stats_tables.txt and R/stats_results.csv" )

step "done"
ls -la results/*.json results/*.jsonl plots/*.pdf | awk '{print "  ", $NF}'
