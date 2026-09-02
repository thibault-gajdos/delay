#!/bin/sh
# Build the LaTeX copy from the Word manuscript.
#
# Three passes over the docx, not one, and the order matters: the converter can only
# rewrite "(Kiani et al., 2008)" into a \parencite once checkbib.py has matched every
# reference to a Zotero key, and checkbib.py can only do that once the converter has
# written out the reference list Word carries. So: convert, match, convert again.
#
#   ./make.sh [manuscript.docx] [bibliography.bib]
#
# Needs python-docx, ImageMagick's convert, and a TeX installation with biblatex and
# biber -- biber, not bibtex. On a machine without them the first pdflatex pass stops at
# \usepackage{biblatex}; there is no partial build to salvage.
set -e
cd "$(dirname "$0")"

DOCX=${1:-../Achard_et_al_R1_v3_tracked-changes_2026-09-02.docx}
BIB=${2:-$HOME/mnt/tools/bib_thib_zotero.bib}

# A snapshot of the library travels with the document: biber must read the same entries
# the citation keys were matched against, and the live export moves under us.
[ "$(readlink -f "$BIB")" = "$(readlink -f bib.bib)" ] || cp "$BIB" bib.bib

echo '--- 1/3  docx -> body.tex, references_verbatim.tex, figures/'
python3 docx2tex.py "$DOCX"

echo '--- 2/3  references -> Zotero keys'
python3 checkbib.py references_verbatim.tex bib.bib

echo '--- 3/3  docx -> body.tex, citations as \parencite'
python3 docx2tex.py "$DOCX"

echo '--- latex'
pdflatex -interaction=nonstopmode main.tex >/dev/null
biber main
pdflatex -interaction=nonstopmode main.tex >/dev/null
pdflatex -interaction=nonstopmode main.tex >/dev/null

echo
echo "errors    $(grep -c '^!' main.log)"
echo "warnings  $(grep -c 'LaTeX Warning' main.log)"
echo "pages     $(pdfinfo main.pdf | awk '/Pages/{print $2}')"
