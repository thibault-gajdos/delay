# Manuscript

`delay.tex` is the manuscript (preamble and body). It is edited by hand.

```sh
./make.sh
```

Needs a TeX installation with biblatex and biber (MacTeX has both).

| file | role |
|---|---|
| `delay.tex` | the manuscript |
| `bib.bib` | references |
| `make.sh` | pdflatex, biber, pdflatex twice |
| `figures/` | figure files; fig05, fig07 and fig08 are written by `scripts/rerun.sh` |
| `table_comparison.tex` | Table 2, model–data comparison; written by `scripts/table_comparison.py` |
| `table_behav.tex`, `table_pma.tex`, `table_mbl.tex`, `table_lrp.tex` | supplementary ANOVA tables; written by `scripts/stats_tables_tex.py` |

Do not edit the `table_*.tex` files by hand: `scripts/rerun.sh` overwrites them.
