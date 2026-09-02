# LaTeX copy of the R1 v3 working manuscript

Generated from the Word tracked-changes copy. `body.tex`, `references_verbatim.tex`,
`citekeys.json` and `figures/` are **output** — they are overwritten on every run, so
edits belong in the docx, not here.

```sh
./make.sh                                   # defaults to the docx one level up
./make.sh ../some-other.docx ~/path/to.bib
```

Requires python-docx, ImageMagick's `convert` (Word embeds three figures as TIFF, which
pdflatex cannot read), and a TeX installation with **biblatex and biber** — biber, not
bibtex. MacTeX has both. The stripped TeX in the Cowork desktop VM has neither, so the
build has to run on the Mac.

| file | role |
|---|---|
| `main.tex` | preamble; `\input{body}` |
| `make.sh` | the whole build |
| `docx2tex.py` | document body: text, change marking, tables, figures |
| `omml2tex.py` | Word OMML equations to LaTeX |
| `checkbib.py` | references against Zotero; writes `citekeys.json` |
| `cites.py` | plain-text APA citations to `\parencite` / `\textcite` |
| `bib.bib` | snapshot of the Zotero export the build used |
| `references_verbatim.tex` | the APA list as Word carried it, for comparison |

## Why the docx is converted three times

The converter can only turn "(Kiani et al., 2008)" into a `\parencite` once every
reference has been matched to a Zotero key, and that matching reads the reference list,
which the converter produces. So: convert, match, convert again. `make.sh` does it.

## Citations

Word held citations as literal text, so nothing connected them to the reference list —
a citation with no reference, or a reference nobody cites, was invisible. They are now
`\parencite` / `\textcite` against Zotero keys, and the list is printed by biblatex from
`bib.bib`, so both become errors.

Consequences worth knowing:

- **Only cited works are printed.** Four references that only the struck (proposed for
  deletion) passages cited no longer appear: Afacan-Seref 2018, Corbett 2023, Murphy
  2016, van Vugt 2014. Citations inside `\del{...}` are deliberately left as plain text.
- **The reference list is uniformly black.** It is generated, so it carries no change
  marking; blue and red still mark the body.
- Multiple citations are sorted alphabetically by biblatex, which is APA's rule and
  occasionally reorders what Word had.

## Change marking

Colour carries the revision state: blue is new since R1_v3, red is inherited, red struck
is proposed for deletion. Links are set `hidelinks` so a citation takes the colour of the
sentence around it rather than hyperref's blue, which would read as new material. For a
clean copy:

```latex
\renewcommand{\new}[1]{#1}\renewcommand{\old}[1]{#1}\renewcommand{\del}[1]{}
```

## Notes

Section numbering is off: the appendix labels its own items A.1–A.10 and the prose cites
those labels, so LaTeX section numbers would compete with them.

Equation numbers come from the docx via `\tag`, not LaTeX's counter, because three
equations are displayed both in Results and in the appendix and the prose cites the
numbers literally.
