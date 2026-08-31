;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "note4"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("eurosym" "") ("bbm" "") ("setspace" "onehalfspacing") ("authblk" "") ("pgf" "") ("xcolor" "") ("mathtools" "") ("inputenc" "utf8") ("color" "") ("calrsfs" "") ("amsfonts" "") ("amssymb" "") ("amsmath" "") ("amsthm" "") ("cellspace" "") ("caption" "") ("booktabs" "") ("multirow" "") ("listings" "") ("babel" "english") ("biblatex" "style=apa") ("geometry" "top=2cm" "bottom=2.2cm" "left=3cm" "right=3cm") ("hyperref" "breaklinks=true" "colorlinks=true" "linkcolor=blue" "citecolor=blue" "filecolor=blue" "urlcolor=blue")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "_estim_body"
    "_tables"
    "article"
    "art11"
    "eurosym"
    "bbm"
    "setspace"
    "authblk"
    "pgf"
    "xcolor"
    "mathtools"
    "inputenc"
    "color"
    "calrsfs"
    "amsfonts"
    "amssymb"
    "amsmath"
    "amsthm"
    "cellspace"
    "caption"
    "booktabs"
    "multirow"
    "listings"
    "babel"
    "biblatex"
    "geometry"
    "hyperref")
   (TeX-add-symbols
    '("refnote" 1)
    '("NN" 1)
    '("RR" 1)
    '("poscite" 1)
    "R"
    "E"
    "Var"
    "Cov"
    "dd"
    "RT"
    "bigO"
    "indep"
    "sgn"
    "Lmot"
    "Lrs"
    "Lmo"
    "toff"
    "trs"
    "sresp")
   (LaTeX-add-labels
    "sec:model"
    "eq:acc"
    "eq:m"
    "eq:maint"
    "eq:z"
    "eq:resp"
    "sec:pred"
    "eq:nu"
    "eq:slope"
    "sec:map"
    "sec:data"
    "sec:estim"
    "sec:reanalysis"
    "sec:reanalysis-tables"
    "sec:corrected"
    "sec:behav"
    "sec:sim"
    "sec:calres"
    "sec:notfollow"
    "sec:responses"
    "sec:figs"
    "fig:traj"
    "fig:single"
    "fig:rates"
    "fig:shape"
    "app:obs"
    "app:nu"
    "eq:kernel"
    "eq:stat"
    "eq:KR1"
    "eq:KR2"
    "sec:methods"
    "sec:software")
   (LaTeX-add-environments
    '("proof" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-bibliographies
    "~/thib/projects/tools/bib_thib_zotero")
   (LaTeX-add-color-definecolors
    "refcol"
    "newcol")
   (LaTeX-add-amsthm-newtheorems
    "definition"
    "axiom"
    "theorem"
    "proposition"
    "remark"))
 :latex)

