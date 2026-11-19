#!/bin/sh
# Build delay.pdf from delay.tex (pdflatex, biber, pdflatex twice).
# Needs biblatex and biber. The tables table_*.tex are written by scripts/rerun.sh.
set -e
cd "$(dirname "$0")"

pdflatex -interaction=nonstopmode delay.tex >/dev/null
biber delay
pdflatex -interaction=nonstopmode delay.tex >/dev/null
pdflatex -interaction=nonstopmode delay.tex >/dev/null

echo
echo "errors    $(grep -c '^!' delay.log)"
echo "warnings  $(grep -c 'LaTeX Warning' delay.log)"
echo "pages     $(pdfinfo delay.pdf | awk '/Pages/{print $2}')"
