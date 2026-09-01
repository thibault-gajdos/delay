#!/bin/bash
# rerun.sh -- everything downstream of the sweep, in order, from a clean state.
#
# The sweep itself and the frontier are NOT here: both are long and both are run sharded
# across cores, so they are launched by hand (see README, "Order of execution"). Everything
# else is deterministic given their outputs and takes about two minutes in total.
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
uv run python3 effect_panel.py | tail -12

step "representative parameter vector (the point the illustrative figures use)"
uv run python3 representative.py | grep -E '"draw"|"beta"' || true

step "figures"
uv run python3 figures.py
uv run python3 figures_space.py

step "the CSV bridge for R"
uv run python3 export_for_R.py

step "statistics"
( cd R && Rscript run_stats.R --with-sim > /dev/null && echo "wrote R/stats_tables.txt and R/stats_results.csv" )

step "done"
ls -la results/*.json results/*.jsonl plots/*.pdf | awk '{print "  ", $NF}'
