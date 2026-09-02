# -*- coding: utf-8 -*-
"""Converts the working manuscript to LaTeX, preserving its two-colour change marking.

The docx is a tracked-changes copy: text present in R1_v3 keeps its original colour (black,
or the red the authors themselves used), text added since is blue, and a few struck runs are
proposed deletions. That marking is the co-authors' reading aid, so it survives the
conversion as three macros -- \\new, \\old and \\del -- whose definitions in the preamble
reproduce the colours and which can be redefined to \\relax for a clean submission copy.

Equations come through omml2tex.py. Two legacy equations, (4) and (10), carry prose inside
the equation without marking it upright, so Word's XML cannot tell their text from their
mathematics; they are overridden explicitly rather than guessed at.
"""
import os, re, sys, collections, subprocess
import docx
from docx.oxml.ns import qn
from docx.table import Table
from docx.text.paragraph import Paragraph
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from omml2tex import to_tex

M = 'http://schemas.openxmlformats.org/officeDocument/2006/math'
NEW, OLD = '1F4E9E', 'FF0000'

OVERRIDE = {                      # equations whose OMML does not separate prose from maths
 '4':  r'm = \begin{cases} s_{\mathrm{resp}}\theta, & \text{if } \theta \text{ was reached},\\'
       r' y(L_{\mathrm{motion}}+D), & \text{otherwise},\end{cases}',
 '10': r'\tilde{y}(t) := \begin{cases} y(t), & \text{Model 1},\\'
       r' y(t)+s_{\mathrm{resp}}\beta, & \text{Model 2}.\end{cases}',
}

def esc(t):
    t = t.replace('\\', r'\textbackslash{}')
    for a, b in [('&', r'\&'), ('%', r'\%'), ('$', r'\$'), ('#', r'\#'), ('_', r'\_'),
                 ('{', r'\{'), ('}', r'\}'), ('~', r'\textasciitilde{}'),
                 ('^', r'\textasciicircum{}')]:
        t = t.replace(a, b)
    for a, b in [('—', '---'), ('–', '--'), ('“', '``'), ('”', "''"), ('‘', '`'),
                 ('’', "'"), ('×', r'$\times$'), ('≤', r'$\le$'), ('≥', r'$\ge$'),
                 ('−', '$-$'), ('≈', r'$\approx$'), ('→', r'$\to$'), ('∎', r'$\square$'),
                 ('²', r'\textsuperscript{2}'), ('³', r'\textsuperscript{3}'),
                 ('⁻¹', r'\textsuperscript{-1}'), ('α', r'$\alpha$'), ('β', r'$\beta$'),
                 ('γ', r'$\gamma$'), ('δ', r'$\delta$'), ('ε', r'$\varepsilon$'),
                 ('η', r'$\eta$'), ('θ', r'$\theta$'), ('κ', r'$\kappa$'),
                 ('λ', r'$\lambda$'), ('μ', r'$\mu$'), ('ν', r'$\nu$'), ('ρ', r'$\rho$'),
                 ('σ', r'$\sigma$'), ('τ', r'$\tau$'), ('φ', r'$\varphi$'),
                 ('ω', r'$\omega$'), ('Δ', r'$\Delta$'), ('Φ', r'$\Phi$'),
                 ('χ', r'$\chi$'), ('ψ', r'$\psi$'), ('π', r'$\pi$'), ('Σ', r'$\Sigma$'),
                 ('\u00a0', '~'), ('\u2009', r'\,'), ('…', r'\ldots{}'),
                 # the rest of what the manuscript actually contains, so that the output
                 # compiles under pdflatex rather than needing a Unicode engine
                 ('\u200b', ''), ('±', r'$\pm$'), ('¶', r'\P{}'), ('†', r'$\dagger$'),
                 ('‡', r'$\ddagger$'), ('§', r'\S{}'), ('°', r'$^{\circ}$'),
                 ('∈', r'$\in$'), ('€', r'\texteuro{}'), ('µ', r'$\mu$'),
                 ('‐', '-'), ('ỹ', r'$\tilde{y}$'),
                 ('ₒ', r'\textsubscript{o}'), ('ₓ', r'\textsubscript{x}'),
                 ('ₘ', r'\textsubscript{m}'), ('ₜ', r'\textsubscript{t}'),
                 ('ᵣ', r'\textsubscript{r}'), ('ᵢ', r'\textsubscript{i}'),
                 ('ₙ', r'\textsubscript{n}'), ('ₑ', r'\textsubscript{e}'),
                 ('ₛ', r'\textsubscript{s}'), ('ₚ', r'\textsubscript{p}'),
                 ('é', r"\'e"), ('è', r'\`e'), ('ê', r'\^e'), ('ï', r'\"i'),
                 ('ô', r'\^o'), ('à', r'\`a'), ('ç', r'\c{c}'), ('ü', r'\"u'),
                 ('ö', r'\"o'), ('ä', r'\"a'), ('ñ', r'\~n'), ('á', r"\'a")]:
        t = t.replace(a, b)
    return t

def run_tex(r):
    txt = ''.join(t.text or '' for t in r.findall(qn('w:t')))
    if not txt:
        return ''
    s = esc(txt)
    # Whitespace at the edges of a run is hoisted OUT of whatever macro wraps the run.
    # TeX drops a space at the start of a macro argument, so "$g$\new{ by the}" sets as
    # "gby the"; Word puts that space in the text run, not in the equation, so leaving it
    # inside the wrapper silently closes up every prose/equation boundary in the paper.
    lead = s[:len(s) - len(s.lstrip(' '))]
    trail = s[len(s.rstrip(' ')):]
    s = s.strip(' ')
    if not s:
        return lead + trail
    rPr = r.find(qn('w:rPr'))
    ital = bold = strike = False; col = None; va = None
    if rPr is not None:
        ital = rPr.find(qn('w:i')) is not None
        bold = rPr.find(qn('w:b')) is not None
        strike = rPr.find(qn('w:strike')) is not None
        c = rPr.find(qn('w:color'))
        if c is not None: col = (c.get(qn('w:val')) or '').upper()
        v = rPr.find(qn('w:vertAlign'))
        if v is not None: va = v.get(qn('w:val'))
    if va == 'superscript': s = r'\textsuperscript{%s}' % s
    elif va == 'subscript': s = r'\textsubscript{%s}' % s
    if ital: s = r'\emph{%s}' % s
    if bold: s = r'\textbf{%s}' % s
    if strike: s = r'\del{%s}' % s
    elif col == NEW: s = r'\new{%s}' % s
    elif col == OLD: s = r'\old{%s}' % s
    return lead + s + trail

def para_tex(p, eqcount):
    """LaTeX for one paragraph.

    Spacing around inline equations is decided on the RAW text, never on the wrapped
    LaTeX: a run that Word ended with a space becomes "\\new{... at }", whose last
    character is a brace, so a check on the output adds a second space. Word keeps the
    spaces in its own runs, so the rule is simply to insert one only when neither side
    already has it.
    """
    items = []          # (raw, tex) with raw = '' for display material
    for ch in p._element:
        if ch.tag == qn('w:r'):
            # A picture lives INSIDE the run (w:r/w:drawing), not beside it as a
            # paragraph child; testing for it at paragraph level silently drops
            # every figure, which is what an earlier version of this did.
            if ch.find(qn('w:drawing')) is not None or ch.find(qn('w:pict')) is not None:
                items.append(('', '%IMAGE%'))
                continue
            raw = ''.join(t.text or '' for t in ch.findall(qn('w:t')))
            items.append((raw, run_tex(ch)))
        elif ch.tag == qn('w:hyperlink'):
            inner = ''.join(run_tex(r) for r in ch.findall(qn('w:r')))
            raw = ''.join(t.text or '' for r in ch.findall(qn('w:r'))
                          for t in r.findall(qn('w:t')))
            tex = (r'\url{%s}' % re.sub(r'\\[a-zA-Z]+\{|\}', '', inner)) \
                  if ('@' in raw or '://' in raw or 'osf.io' in raw) else inner
            items.append((raw, tex))
        elif ch.tag == '{%s}oMath' % M:
            body, _ = to_tex(ch)
            items.append(('\x01', '$%s$' % body.strip()))
        elif ch.tag == '{%s}oMathPara' % M:
            for om in ch.findall('{%s}oMath' % M):
                body, num = to_tex(om)
                if num in OVERRIDE: body = OVERRIDE[num]
                if num:
                    eqcount.append(num)
                    items.append(('', '\n\\begin{equation}\n%s\n\\tag{%s}\n'
                                      '\\end{equation}\n' % (body.strip(), num)))
                else:
                    items.append(('', '\n\\[ %s \\]\n' % body.strip()))
        elif ch.tag == qn('w:drawing'):
            items.append(('', '%IMAGE%'))
    out = []
    for i, (raw, tex) in enumerate(items):
        if not tex:
            continue
        if out:
            praw = items[i-1][0] if i else ''
            need = (raw == '\x01' or praw == '\x01')
            has = (praw.endswith((' ', '\n')) or raw.startswith((' ', '\n'))
                   or praw == '' or raw == '')
            starts_punct = raw[:1] in ',.;:)]}?!' or (raw == '\x01' and False)
            if need and not has and not starts_punct:
                out.append(' ')
        out.append(tex)
    return ''.join(out)


HEAD = {'Heading 1': 'section', 'Heading 2': 'subsection',
        'Heading 3': 'subsubsection', 'Heading 4': 'paragraph'}

def table_tex(t, eqcount):
    """One tabular per Word table.

    Cells go through para_tex rather than through esc(cell.text): Table A1 carries two
    inline equations that plain-text extraction would drop silently.

    Column types are chosen from the content, because left-aligned columns take their
    natural width and overrun the margin. A table carrying prose (some column runs past
    forty characters) gets wrapping p-columns in proportion to its longest cells; that
    is Table A1, which is 418pt too wide otherwise. A table of short labels and figures
    would look absurd wrapped, so it keeps its natural columns and is scaled down --
    but only if it does not fit, which \\resizebox decides at typesetting time rather
    than from a character count here, since character widths are not uniform.
    """
    # Horizontal merges are read from w:gridSpan, not inferred from repeated text.
    # python-docx's row.cells repeats a merged cell once per grid column it covers, which
    # turns Table 2's single spanning "Period" header into seven copies of the word; and
    # the same table has genuinely repeated adjacent cells ("--"), so collapsing equal
    # neighbours instead would merge cells Word never merged.
    rows, plain, spans = [], [], []
    for tr in t._tbl.tr_lst:
        cells, raw, sp = [], [], []
        for tc in tr.tc_lst:
            cells.append(' '.join(para_tex(Paragraph(p, t), eqcount).strip()
                                  for p in tc.findall(qn('w:p'))).strip())
            raw.append(' '.join((p.text or '').strip()
                                for p in tc.iter(qn('w:t'))).strip())
            gs = tc.find(qn('w:tcPr'))
            gs = gs.find(qn('w:gridSpan')) if gs is not None else None
            sp.append(int(gs.get(qn('w:val'))) if gs is not None else 1)
        rows.append(cells); plain.append(raw); spans.append(sp)
    n = max(sum(s) for s in spans)
    for i in range(len(rows)):
        pad = n - sum(spans[i])
        rows[i] += [''] * pad; plain[i] += [''] * pad; spans[i] += [1] * pad
    # column widths: a spanning cell contributes to none of them
    wid = [1] * n
    for r, sp, pl in zip(rows, spans, plain):
        j = 0
        for k, s in enumerate(sp):
            if s == 1:
                wid[j] = max(wid[j], len(pl[k]))
            j += s
    wrap = max(wid) > 40
    if wrap:
        tot = float(sum(wid))
        spec = ''.join(r'>{\raggedright\arraybackslash}p{%.3f\linewidth}'
                       % max(0.93 * wid[j] / tot, 0.04) for j in range(n))
    else:
        spec = 'l' * n
    out = [r'\begin{center}\small']
    if not wrap:
        out.append(r'\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{%')
    def line(r, sp):
        return ' & '.join(c if s == 1 else r'\multicolumn{%d}{c}{%s}' % (s, c)
                          for c, s in zip(r, sp))
    out += [r'\begin{tabular}{%s}' % spec, r'\hline']
    out.append(line(rows[0], spans[0]) + r' \\ \hline')
    for r, sp in zip(rows[1:], spans[1:]):
        out.append(line(r, sp) + r' \\')
    out += [r'\hline', r'\end{tabular}']
    if not wrap:
        out.append('}')
    out.append(r'\end{center}')
    return '\n'.join(out)

def convert(path, out, figdir):
    d = docx.Document(path)
    body = d.element.body
    eqcount = []
    imgs = []
    # extract the images, in document order
    os.makedirs(figdir, exist_ok=True)
    seen = {}
    for i, blip in enumerate(body.findall('.//' + qn('a:blip'))):
        rid = blip.get(qn('r:embed'))
        if rid not in seen:
            part = d.part.related_parts[rid]
            ext = str(part.partname).rsplit('.', 1)[-1].lower()
            fn = 'fig%02d.%s' % (len(seen) + 1, ext)
            open(os.path.join(figdir, fn), 'wb').write(part.blob)
            if ext in ('tif', 'tiff', 'emf', 'wmf'):
                # pdflatex reads only pdf/png/jpg; Word happily embeds tiff.
                png = fn.rsplit('.', 1)[0] + '.png'
                subprocess.run(['convert', os.path.join(figdir, fn),
                                os.path.join(figdir, png)], check=True)
                os.remove(os.path.join(figdir, fn))   # keep only what the document uses
                fn = png
            seen[rid] = fn
        imgs.append(seen[rid])
    L = []
    it = iter(imgs)
    ip = it_ = 0
    for ch in body.iterchildren():
        if ch.tag == qn('w:p'):
            p = Paragraph(ch, d)
            st = p.style.name
            tex = para_tex(p, eqcount)
            while '%IMAGE%' in tex:
                fn = next(it)
                tex = tex.replace('%IMAGE%',
                    '\n\\begin{center}\\includegraphics[width=\\linewidth]{figures/%s}'
                    '\\end{center}\n' % fn, 1)
            if not tex.strip():
                continue
            if st in HEAD:
                plain = re.sub(r'\\(new|old|del|emph|textbf)\{', '', tex).replace('}', '')
                L.append('\n\\%s{%s}\n' % (HEAD[st], plain.strip()))
            elif st == 'Bibliography':
                # trailing newline -> blank line after the join, so each reference
                # is its own paragraph and can take a hanging indent
                L.append('%REF%' + tex.strip() + '\n')
            else:
                L.append(tex.strip() + '\n')
        elif ch.tag == qn('w:tbl'):
            L.append(table_tex(Table(ch, d), eqcount))
    text = '\n'.join(L)
    # The reference list Word carries is written out beside body.tex rather than into it.
    # checkbib.py reads it (it is the only place the APA text survives) and body.tex gets
    # \printbibliography instead, so the list is generated from the .bib by the same keys
    # the citations use. A reference nobody cites then disappears and a citation with no
    # reference becomes an error, which is the whole point of the exercise.
    refdir = os.path.dirname(os.path.abspath(out))
    m = re.search(r'(%REF%.*?)(?=\n\n\\section|\Z)', text, re.S)
    if m:
        verbatim = '\\begin{references}\n%s\n\\end{references}' % m.group(1).replace('%REF%', '')
        open(os.path.join(refdir, 'references_verbatim.tex'), 'w').write(verbatim + '\n')
        text = text[:m.start()] + '\\printbibliography[heading=none]\n' + text[m.end():]

    keyfile = os.path.join(refdir, 'citekeys.json')
    if os.path.exists(keyfile):
        import cites
        text, rw = cites.rewrite(text, cites.load_keys(keyfile))
        print('citations rewritten: %d   left as text: %d' % (len(rw.done), len(rw.skipped)))
        for s in rw.skipped:
            print('   NOT RECOGNISED: %s' % s)
    else:
        print('!! no citekeys.json yet -- citations left as plain text.')
        print('!! run checkbib.py, then run this again (see README).')
    open(out, 'w').write(text)
    print('%d paragraphs of body, %d numbered equations, %d images, %d tables'
          % (len(L), len(eqcount), len(seen), len(d.tables)))
    return eqcount

if __name__ == '__main__':
    # Default: the tracked-changes copy one directory up. Pass a path to override.
    #   python3 docx2tex.py [manuscript.docx]
    here = os.path.dirname(os.path.abspath(__file__))
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        here, '..', 'Achard_et_al_R1_v3_tracked-changes_2026-09-02.docx')
    convert(src, os.path.join(here, 'body.tex'), os.path.join(here, 'figures'))
