# LaTeX copy of the R1 v3 working manuscript

**As of 2026-09-03, `delay.tex` is hand-edited directly and is the source of truth.**
It started life converted from the Word tracked-changes copy (`docx2tex.py`, below), but
the docx is no longer kept in sync and must not be treated as authoritative — edits
belong in `delay.tex`, not the docx. `main.tex` and `body.tex` no longer exist: they were
merged into the single `delay.tex` on 2026-09-03 (preamble + document body in one file).

```sh
./make.sh
```

Requires a TeX installation with **biblatex and biber** — biber, not bibtex. MacTeX
has both. The stripped TeX in the Cowork desktop VM has neither, so the build has to run
on the Mac.

| file | role |
|---|---|
| `delay.tex` | the whole document: preamble + body. Hand-edited. |
| `make.sh` | compiles `delay.tex` (pdflatex + biber + pdflatex x2) |
| `bib.bib` | snapshot of the Zotero export the original conversion used |

The files below produced the original docx -> LaTeX conversion and are kept for
reference, but are **not part of the build any more** and must not be run against
`delay.tex` — there is no separate `body.tex` for them to write into, and running them
would overwrite hand edits with a regeneration from the (now stale) docx.

| file | historical role |
|---|---|
| `docx2tex.py` | document body: text, change marking, tables, figures |
| `omml2tex.py` | Word OMML equations to LaTeX |
| `checkbib.py` | references against Zotero; wrote `citekeys.json` |
| `cites.py` | plain-text APA citations to `\parencite` / `\textcite` |
| `references_verbatim.tex` | the APA list as Word carried it, for comparison |

## Why the docx used to be converted three times (historical)

The converter could only turn "(Kiani et al., 2008)" into a `\parencite` once every
reference had been matched to a Zotero key, and that matching read the reference list,
which the converter produced. So: convert, match, convert again. No longer relevant now
that `delay.tex` is edited directly, but kept here in case the docx pipeline is ever
revived.

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

Colour carries the revision state, exactly as the docx does it: **`#1F4E9E` blue** for
material new since R1_v3, **`#FF0000` red** for text inherited from it, red struck through
for proposed deletions. Those are the only two colours the text uses; a pixel audit of the
rendered PDF finds nothing else but their antialiasing.

It covers prose, **equations** (`m:r` runs, which an earlier version of the converter
never looked at, so all 218 blue equations came out black) and **headings** (38 of the 58
are coloured in the docx; an earlier version stripped the macros there to keep `\section`
arguments simple, and lost the marking on all of them). A sample of 258 runs checked
against the docx run by run now agrees on every one.

Two departures worth knowing:

- The **figures use their own palette** — `#3182BD` and `#08306B` for blue, `#C1272D` and
  `#C00000` for red — because Model 1 and Model 2 are colour-coded as data, independently
  of the revision marking. So the document really does contain two blues; only one of them
  means "new".
- 21 struck runs are `#B00000` in the docx, a darker red, and are drawn here in the single
  deletion red. They are all strikethrough, so nothing is lost but the shade.

Links are `hidelinks`, so a citation takes the colour of its sentence rather than
hyperref's blue, which would read as new material. For a clean copy:

```latex
\renewcommand{\new}[1]{#1}\renewcommand{\old}[1]{#1}\renewcommand{\del}[1]{}
```

## Preamble

Follows `note4.tex`: A4, `onehalfspacing`, its geometry (3 cm sides), `babel`, `mathtools`,
`booktabs`, `\allowdisplaybreaks`, and its macro block — `\R \E \P \Var \Cov \dd \RT
\bigO \1 \indep \sgn \poscite`. The converter emits `\Var`, `\Cov` and `\dd` for the
operators Word stores as upright runs, so the generated maths reads like the note's.

Two deliberate departures:

- **No `\mathtoolsset{showonlyrefs}`.** It prints a number only on an equation some
  `\eqref` points at. Here the numbers come from the docx through `\tag` and the prose
  cites them as literal text — "Equation (21)" — never with `\eqref`. Every number in the
  paper would vanish and every reference to one would dangle.
- **`hidelinks` instead of `colorlinks` with blue citations.** This copy uses colour to
  mark revisions; a blue `\cite` inside a red sentence reads as new material. Links stay
  live and take the colour of the text around them.

`\P` is the probability measure, as in the note, so the pilcrows in Table 2's footnotes
are `\textparagraph` rather than `\P`.

## Equation numbers

Every equation carries `\label{eq:N}` and every reference in the prose is `\eqref` (59 of
them, none left as literal text), so a reference prints whatever the equation actually
got and a reference to an equation that no longer exists is a warning rather than a wrong
number.

LaTeX numbers the equations, in the order they appear: a clean 1–27. Three of them —
Kac–Rice and its two corollaries — are displayed twice, in Results and again in the
appendix where they are derived. They are `\label`led once, on the first appearance, and
the repeat echoes that number through `\tag{\ref{...}}` rather than taking one of its own:
the same equation keeps one number.

That means one of the two places shows a number out of local order — the appendix now
runs 22, 23, 24, then 11 — which is unavoidable when an equation appears twice. Setting
`WORD_NUMBERS = True` in `docx2tex.py` puts the docx's own numbers back instead, which
moves the out-of-order number into Results and keeps this copy and the Word copy agreeing
on what "(21)" means. The references do not care either way; they are symbolic.

## Notes

Section numbering is off: the appendix labels its own items A.1–A.10 and the prose cites
those labels, so LaTeX section numbers would compete with them.
