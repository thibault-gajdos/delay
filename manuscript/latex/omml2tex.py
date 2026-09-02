# -*- coding: utf-8 -*-
"""Office MathML to LaTeX, for the constructs this manuscript actually uses.

Anything unrecognised raises rather than being dropped silently: a converter that quietly
loses a subscript is worse than one that stops. The inverse of mathml.py, written against
the same set of structures plus the ones the pre-existing equations (1) to (10) use --
equation arrays, accents, and Word's '#' equation-number separator."""
import re
from lxml import etree

M = 'http://schemas.openxmlformats.org/officeDocument/2006/math'
def q(t): return '{%s}%s' % (M, t)

NUM = [None]        # equation number found inside a nested eqArr; read by to_tex

# Function names Word stores as ordinary italic runs. Left alone they set as products of
# variables -- "max" as m times a times x -- so they are recognised here by name.
FUNCS = {'max': r'\max', 'min': r'\min', 'exp': r'\exp', 'log': r'\log',
         'ln': r'\ln', 'sup': r'\sup', 'inf': r'\inf', 'det': r'\det',
         'sin': r'\sin', 'cos': r'\cos', 'tan': r'\tan'}

GREEK = {'α':r'\alpha','β':r'\beta','γ':r'\gamma','δ':r'\delta','ε':r'\varepsilon',
         'ζ':r'\zeta','η':r'\eta','θ':r'\theta','ι':r'\iota','κ':r'\kappa',
         'λ':r'\lambda','μ':r'\mu','ν':r'\nu','ξ':r'\xi','π':r'\pi','ρ':r'\rho',
         'σ':r'\sigma','τ':r'\tau','υ':r'\upsilon','φ':r'\varphi','χ':r'\chi',
         'ψ':r'\psi','ω':r'\omega','Δ':r'\Delta','Φ':r'\Phi','Σ':r'\Sigma',
         'Ω':r'\Omega','Γ':r'\Gamma','Λ':r'\Lambda','Θ':r'\Theta'}
OPS = {'−':'-','–':'-','—':'---','×':r'\times','·':r'\cdot','≤':r'\le','≥':r'\ge',
       '≠':r'\ne','≈':r'\approx','→':r'\to','∞':r'\infty','∈':r'\in','∫':r'\int',
       '∑':r'\sum','√':r'\sqrt','±':r'\pm','∼':r'\sim','↦':r'\mapsto','∂':r'\partial',
       '≲':r'\lesssim','≳':r'\gtrsim','∎':r'\qedhere','⇒':r'\Rightarrow','…':r'\ldots',
       '′':"'", '̇':'', '̅':'', '‰':r'\permil', '\u2009':r'\,', '\u00a0':'~',
       '⁻':'^-', '¹':'^1', '²':'^2', '³':'^3', '≔':':='}
ACC = {'̇': r'\dot', '̅': r'\bar', '̂': r'\hat', '̃': r'\tilde'}

def esc(s, upright=False):
    # every macro gets a trailing space: without it "a<=t<=b" compiles to \let\le, which
    # is a control sequence LaTeX knows and will not complain about until much later
    out = []
    for ch in s:
        if ch in GREEK: out.append(GREEK[ch] + ' ')
        elif ch in OPS: out.append(OPS[ch] + (' ' if OPS[ch].startswith('\\') else ''))
        elif ch in '#$%&_{}': out.append('\\' + ch)
        elif ch == '\\': out.append(r'\backslash ')
        else: out.append(ch)
    t = ''.join(out)
    t = re.sub(r'\s+', ' ', t)
    if not (upright and re.search('[A-Za-z]', s)):
        return t                                   # operators and digits need no wrapper
    # prose inside an equation -- "if theta was reached" -- is text, not an operator name
    if ' ' in s.strip() and len(re.findall('[A-Za-z]', s)) > 3:
        return r'\text{%s}' % t
    return r'\mathrm{%s}' % t

def conv(el):
    tag = el.tag.split('}')[1]
    if tag == 'r':
        upright = el.find(q('rPr')) is not None and el.find(q('rPr')).find(q('nor')) is not None
        txt = ''.join(t.text or '' for t in el.findall(q('t')))
        # A run of two or more capitals is an abbreviation, not a product of variables:
        # "SD" is the standard deviation, not S times D. Word marks most of them upright
        # but not all -- Table A1's SD_eps is one of the misses -- so the rule is applied
        # here rather than trusting m:nor.
        if not upright and re.fullmatch('[A-Z]{2,}', txt.strip()):
            upright = True
        if txt.strip() in FUNCS:
            return FUNCS[txt.strip()] + ' '
        return esc(txt, upright)
    if tag in ('oMath', 'e', 'num', 'den', 'sub', 'sup', 'deg', 'oMathPara', 'lim', 'fName'):
        return ''.join(conv(c) for c in el if c.tag.startswith('{' + M) and
                       c.tag.split('}')[1] not in ('rPr', 'ctrlPr', 'fPr', 'sSubPr',
                       'sSupPr', 'sSubSupPr', 'radPr', 'dPr', 'naryPr', 'accPr',
                       'eqArrPr', 'oMathParaPr', 'limLowPr', 'limUppPr', 'barPr',
                       'groupChrPr', 'mPr', 'funcPr'))
    if tag == 'f':
        typ = el.find(q('fPr'))
        lin = typ is not None and typ.find(q('type')) is not None and \
              typ.find(q('type')).get(q('val')) == 'lin'
        n, dd = conv(el.find(q('num'))), conv(el.find(q('den')))
        return (r'%s/%s' % (n, dd)) if lin else (r'\frac{%s}{%s}' % (n, dd))
    if tag == 'rad':
        pr = el.find(q('radPr'))
        hide = pr is not None and pr.find(q('degHide')) is not None
        deg = conv(el.find(q('deg'))) if el.find(q('deg')) is not None else ''
        body = conv(el.find(q('e')))
        return r'\sqrt{%s}' % body if (hide or not deg) else r'\sqrt[%s]{%s}' % (deg, body)
    if tag == 'sSub':  return '{%s}_{%s}' % (conv(el.find(q('e'))), conv(el.find(q('sub'))))
    if tag == 'sSup':  return '{%s}^{%s}' % (conv(el.find(q('e'))), conv(el.find(q('sup'))))
    if tag == 'sSubSup':
        return '{%s}_{%s}^{%s}' % (conv(el.find(q('e'))), conv(el.find(q('sub'))),
                                   conv(el.find(q('sup'))))
    if tag == 'd':
        pr = el.find(q('dPr')); b, e = '(', ')'
        if pr is not None:
            if pr.find(q('begChr')) is not None: b = pr.find(q('begChr')).get(q('val'))
            if pr.find(q('endChr')) is not None: e = pr.find(q('endChr')).get(q('val'))
        inner = ''.join(conv(c) for c in el.findall(q('e')))
        L = {'(': r'\left(', '[': r'\left[', '{': r'\left\{', '|': r'\left|', '': r'\left.'}
        R = {')': r'\right)', ']': r'\right]', '}': r'\right\}', '|': r'\right|', '': r'\right.'}
        return '%s%s%s' % (L.get(b, r'\left('), inner, R.get(e, r'\right)'))
    if tag == 'nary':
        pr = el.find(q('naryPr'))
        ch = pr.find(q('chr')).get(q('val')) if pr is not None and pr.find(q('chr')) is not None else '∫'
        op = {'∫': r'\int', '∑': r'\sum', '∏': r'\prod'}.get(ch, r'\int')
        lo = conv(el.find(q('sub'))) if el.find(q('sub')) is not None else ''
        hi = conv(el.find(q('sup'))) if el.find(q('sup')) is not None else ''
        body = conv(el.find(q('e'))) if el.find(q('e')) is not None else ''
        s = op
        if lo: s += '_{%s}' % lo
        if hi: s += '^{%s}' % hi
        return s + ' ' + body
    if tag == 'acc':
        pr = el.find(q('accPr'))
        ch = pr.find(q('chr')).get(q('val')) if pr is not None and pr.find(q('chr')) is not None else '̇'
        return '%s{%s}' % (ACC.get(ch, r'\dot'), conv(el.find(q('e'))))
    if tag == 'bar':
        return r'\overline{%s}' % conv(el.find(q('e')))
    if tag == 'eqArr':
        # Word's equation-number separator is "#", and in a multi-row array it appears
        # once per row -- equation (3) carries "#(3)" on its first line and a bare "#" on
        # its second. Both are markup, not content: the number is lifted out for \tag and
        # the rows are wrapped in aligned, since "\\" is not legal in bare display maths.
        rows = []
        for c in el.findall(q('e')):
            row = conv(c)
            i = row.find(r'\#')
            if i >= 0:
                # Everything from the "#" on is the number, however Word wrapped it:
                # a bare "22", a parenthesised "(3)", or a \left(22\right) delimiter.
                m = re.search(r'\d+', row[i:])
                if m:
                    NUM[0] = m.group(0)
                row = row[:i]
            row = row.strip()
            if row:
                rows.append(row)
        if len(rows) > 1:
            return r'\begin{aligned}%s\end{aligned}' % (' \\\\ '.join(rows))
        return rows[0] if rows else ''
    if tag == 'func':                      # a named function applied to an argument
        name = conv(el.find(q('fName'))) if el.find(q('fName')) is not None else ''
        return '%s %s' % (name, conv(el.find(q('e'))))
    if tag == 'fName':
        return ''.join(conv(c) for c in el)
    if tag == 'groupChr':
        return conv(el.find(q('e')))
    if tag in ('limLow', 'limUpp'):
        return conv(el.find(q('e')))
    if tag == 'm':      # matrix
        return ' \\\\ '.join(''.join(conv(c) for c in row.findall(q('e')))
                             for row in el.findall(q('mr')))
    if tag in ('rPr','ctrlPr','fPr','sSubPr','sSupPr','sSubSupPr','radPr','dPr','naryPr',
               'accPr','eqArrPr','oMathParaPr','barPr','groupChrPr','mPr','funcPr','t'):
        return ''
    raise ValueError('unhandled OMML element: m:' + tag)

TAGNUM = re.compile(r'\\#\s*\\?\(?(\d+)\)?\\?#?\s*$|#\s*\(?(\d+)\)?\s*#?\s*$')

def to_tex(el):
    """LaTeX for one m:oMath / m:oMathPara, and the equation number if it carries one."""
    NUM[0] = None                    # filled in by the eqArr branch, which sees it first
    s = conv(el)
    s = re.sub(r'\s+', ' ', s).strip()
    m = re.search(r'\\#\s*\\left\((\d+)\\right\)\s*$|\\#\s*\(?(\d+)\)?\s*$', s)
    num = NUM[0]
    if m:
        num = m.group(1) or m.group(2)
        s = s[:m.start()].strip()
    return s, num
