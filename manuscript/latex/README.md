# LaTeX copy of the R1 v3 working manuscript

Generated from the Word tracked-changes copy. `body.tex` and `figures/` are **output** —
they are overwritten on every run, so edits belong in the docx, not here.

```sh
python3 docx2tex.py ../Achard_et_al_R1_v3_tracked-changes_2026-09-02.docx
pdflatex main.tex && pdflatex main.tex        # twice, for the page references
```

Requires python-docx, pdflatex with `amsmath`, `ulem`, `lineno`, `microtype`, and
ImageMagick's `convert` (Word embeds three figures as TIFF, which pdflatex cannot read).

| file | role |
|---|---|
| `main.tex` | preamble; `\input{body}` |
| `body.tex` | generated |
| `figures/` | generated: 11 images extracted from the docx |
| `omml2tex.py` | Word OMML equations to LaTeX |
| `docx2tex.py` | document body: text, change marking, tables, figures |

## Change marking

Colour carries the revision state, as in the Word copy: blue is new since R1_v3, red is
inherited, red struck is proposed for deletion. For a clean copy, add to the preamble:

```latex
\renewcommand{\new}[1]{#1}\renewcommand{\old}[1]{#1}\renewcommand{\del}[1]{}
```

## Notes

Section numbering is off. The appendix labels its own items A.1-A.10 and the prose cites
those labels; LaTeX section numbers would compete with them.

Equation numbers come from the docx via `\tag`, not from LaTeX's counter, because three
equations are displayed both in Results and in the appendix and the prose cites the
numbers literally.
