#!/bin/sh
# Build the PDF from delay.tex.
#
# delay.tex is now hand-edited directly (2026-09-03: the docx is no longer the source of
# truth -- see git log). This script therefore no longer regenerates anything from the
# docx; it only runs the LaTeX/biber passes. docx2tex.py, omml2tex.py, checkbib.py and
# cites.py are kept for reference (they produced the original conversion) but are not
# part of this build any more, and must not be run against delay.tex -- there is no
# "body.tex" for them to write into, and pointing them at delay.tex would clobber every
# hand edit made since the conversion.
#
#   ./make.sh
#
# Needs a TeX installation with biblatex and biber -- biber, not bibtex.
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
