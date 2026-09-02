# -*- coding: utf-8 -*-
"""Turn the manuscript's plain-text APA citations into biblatex commands.

Word holds citations as literal text -- "(Kiani et al., 2008)" -- so nothing connects
them to the reference list. This rewrites them as \\parencite / \\textcite against the
Zotero keys, which is what makes the list come out of the .bib and makes a citation with
no reference (or a reference nobody cites) an error rather than a silent inconsistency.

Two things it deliberately does not do:

  - It leaves everything inside \\del{...} alone. Those are the passages proposed for
    deletion, and two of the references they cite are themselves struck, so they have no
    key to cite; striking through a live \\parencite would also mean biblatex counted a
    deleted citation as a use.
  - It does not guess. Every citation it does not recognise is reported by name, so the
    count of rewritten citations can be checked against the count in the source.

A citation that straddles a change-marking boundary -- "(Coles et al., 1985;} \\new{Servant
et al., 2021)" -- is emitted whole in the colour of its last fragment. A citation is one
unit and cannot be half one colour; the alternative is to leave it as text.
"""
import re, json, os

PRE = r'(?:e\.g\.|i\.e\.|see also|see|cf\.|but see|for example|for review)'
YEAR = r'(?:19|20)\d\d[a-z]?'
# The author phrase is not parsed. Author names run from "Wang" to "Pfurtscheller &
# Lopes da Silva" to "Aza\"is & Wschebor", and a regex precise enough to accept all of
# them accepts most of the sentence too. Instead everything before the year is taken as
# the phrase and looked up in the key map, which is the oracle: if the surname and year
# resolve to a Zotero key the reading was right, and if they do not the citation is
# reported rather than guessed at.
AUTHORS = r"(?:[A-Z\\]|van |von |de |da )[^;()]*?"

# The narrative form "Dendauw et al. (2024)" has no brackets to bound it, so here the
# phrase must be spelled out: a surname, optionally with a lowercase particle and a TeX
# accent, optionally followed by "et al." or a second surname. Reusing the loose pattern
# above would let "the voltage difference between each electrode pair was computed
# (2025)" read as an author.
SURNAME = r"(?:(?:van|von|de|del|della|da|di|dos|le|la)\s+)?" \
          r"[A-Z][A-Za-zÀ-ÿ'’\-]*(?:\\[\"'^`~][a-z][A-Za-zÀ-ÿ'’\-]*)?"
NARRATIVE = r"%s(?:\s+et\s+al\.|\s*\\?&\s*%s|\s+and\s+%s)?" % (SURNAME, SURNAME, SURNAME)


def load_keys(path=None):
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'citekeys.json')
    raw = json.load(open(path))
    return {tuple(k.split('|')): v for k, v in raw.items()}


def _norm_surname(s):
    """The first author's surname, spelled as the key map spells it: lower case, no TeX
    accent macros, no second author, no "et al."."""
    s = re.sub(r'\\[\"\'^`~]', '', s)                   # Aza\"is -> Azais
    s = re.sub(r'\s*(?:\\?&|\band\b)\s.*$', '', s)      # drop everything from "& ..."
    s = re.sub(r'\s+et\s+al\.?\s*$', '', s)
    return s.strip().rstrip(',').lower()


class Rewriter:
    def __init__(self, keys):
        self.keys = keys
        self.done = []          # (original, replacement)
        self.skipped = []       # citations left as text, with the reason

    def key(self, surname, year):
        k = self.keys.get((_norm_surname(surname), year))
        return k

    # ---------------------------------------------------------------- one citation
    def one(self, part):
        """"e.g., Kiani et al., 2008, Experiment 2" -> (prenote, [keys], postnote)."""
        part = part.strip()
        pre = []
        while True:                      # "e.g., 3.2--51.2\%:" is two prenotes in a row
            m = re.match(r'^(%s,?\s*)' % PRE, part) or re.match(r'^([^:,]{1,24}:\s*)', part)
            if not m:
                break
            pre.append(m.group(1).strip().rstrip(',:').strip())
            part = part[m.end():]
        pre = ', '.join(p for p in pre if p)
        m = re.match(r'^(%s),\s*(%s)((?:\s*,\s*%s)*)(.*)$' % (AUTHORS, YEAR, YEAR), part)
        if not m:
            # Table A1 writes "Dendauw et al. 2024" with no comma before the year
            m = re.match(r'^(%s)\s+(%s)()(.*)$' % (AUTHORS, YEAR), part)
            if not m:
                return None
        surname, y1, more, tail = m.groups()
        years = [y1] + re.findall(YEAR, more or '')
        keys = [self.key(surname, y) for y in years]
        if any(k is None for k in keys):
            return None
        post = tail.strip().lstrip(',').strip()
        return pre, keys, post

    # ---------------------------------------------------------------- emit
    def emit(self, parts, paren):
        """parts is a list of (pre, keys, post). biblatex takes optional arguments
        positionally: one bracket is a postnote, two are prenote then postnote."""
        flat = [(pre, k, post if k == keys[-1] else '')
                for pre, keys, post in parts for k in keys]
        cmd = 'parencite' if paren else 'textcite'
        if len(flat) == 1:
            pre, k, post = flat[0]
            if pre:
                return '\\%s[%s][%s]{%s}' % (cmd, pre, post, k)
            return ('\\%s[%s]{%s}' % (cmd, post, k)) if post else '\\%s{%s}' % (cmd, k)
        out = ['\\%ss' % cmd]
        for pre, k, post in flat:
            out.append(('[%s][%s]' % (pre, post)) if pre else
                       ('[%s]' % post if post else ''))
            out.append('{%s}' % k)
        return ''.join(out)

    # ---------------------------------------------------------------- the body
    def run(self, text):
        protected = _del_spans(text)

        def inside(i):
            return any(a <= i < b for a, b in protected)

        # (1) parenthetical groups, possibly straddling a change-marking boundary
        boundary = r'(?:\}\s*\\(?:new|old)\{)'
        pat = re.compile(r'\((?P<body>(?:[^()]|%s)*?%s(?:[^()]|%s)*?)\)'
                         % (boundary, YEAR, boundary))

        def sub_paren(m):
            if inside(m.start()):
                return m.group(0)
            raw = m.group('body')
            # the boundary tokens vanish: the citation is emitted in one colour, that of
            # the fragment that closes it
            flat = re.sub(boundary, ' ', raw)
            pieces = [p for p in flat.split(';') if p.strip()]
            parsed = [self.one(p) for p in pieces]
            if not parsed or any(p is None for p in parsed):
                if re.search(r'[A-Z][a-z]+[^;]*,\s*%s' % YEAR, flat):
                    self.skipped.append(raw.strip()[:80])
                return m.group(0)
            rep = self.emit(parsed, paren=True)
            if re.search(boundary, raw):
                rep = '}\\new{' + rep      # close the old group, reopen in the new colour
            self.done.append((m.group(0)[:70], rep))
            return rep

        text = pat.sub(sub_paren, text)

        # (2) narrative "Author et al. (2024)"
        def sub_narr(m):
            if inside(m.start()):
                return m.group(0)
            got = self.one('%s, %s' % (m.group(1), m.group(2)))
            if got is None:
                self.skipped.append(m.group(0)[:80])
                return m.group(0)
            rep = self.emit([got], paren=False)
            self.done.append((m.group(0), rep))
            return rep

        text = re.sub(r'(?<![A-Za-z])(%s)\s*\((%s)\)' % (NARRATIVE, YEAR), sub_narr, text)
        return text


def _del_spans(text):
    """Character ranges covered by \\del{...}, brace-matched."""
    spans = []
    for m in re.finditer(r'\\del\{', text):
        depth, i = 1, m.end()
        while i < len(text) and depth:
            if text[i] == '{':
                depth += 1
            elif text[i] == '}':
                depth -= 1
            i += 1
        spans.append((m.start(), i))
    return spans


def rewrite(text, keys=None):
    r = Rewriter(keys or load_keys())
    return r.run(text), r
