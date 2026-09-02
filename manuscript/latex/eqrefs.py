# -*- coding: utf-8 -*-
"""Turn the manuscript's literal equation numbers into \\eqref.

Word writes "(21)" as text, so the number in the prose and the number on the equation are
two unrelated strings that happen to agree. Insert an equation and they stop agreeing,
silently. This makes every reference symbolic: the equation carries \\label{eq:N} and the
prose carries \\eqref{eq:N}, so LaTeX prints whatever the equation actually got and a
reference to an equation that no longer exists becomes a warning instead of a wrong number.

The printed numbers still come from the docx through \\tag, so this copy and the Word copy
agree on what equation (21) is; docx2tex.WORD_NUMBERS switches that off and lets LaTeX
number the equations instead. Either way this module is unchanged, because the references
are symbolic by then.

Every parenthesised one- or two-digit number in the body is an equation reference; there
is nothing else of that shape in the paper. The check is still made against the set of
numbers that actually exist, and anything outside it is left alone and reported.
"""
import re

# The equation environments themselves, and their own \tag{N}, must not be touched.
BLOCK = re.compile(r'\\begin\{equation\}.*?\\end\{equation\}|\\\[.*?\\\]', re.S)


def rewrite(text, numbers):
    """Replace "(N)" by \\eqref{eq:N} outside display maths. `numbers` is the set of
    equation numbers that exist, as strings."""
    spans = [m.span() for m in BLOCK.finditer(text)]

    def inside(i):
        return any(a <= i < b for a, b in spans)

    done, unknown = [], set()

    def sub(m):
        if inside(m.start()):
            return m.group(0)
        n = m.group(1)
        if n not in numbers:
            unknown.add(n)
            return m.group(0)
        done.append(n)
        return r'\eqref{eq:%s}' % n

    out = re.sub(r'\((\d{1,2})\)', sub, text)

    def sub_bare(m):
        # "Equation 2 describes ..." carries no parentheses, so it takes \ref rather than
        # \eqref, which would add a pair the sentence never had.
        if inside(m.start()) or m.group(2) not in numbers:
            return m.group(0)
        done.append(m.group(2))
        return '%s~\\ref{eq:%s}' % (m.group(1), m.group(2))

    out = re.sub(r'(Equations?|Eqs?\.)\s+(\d{1,2})\b', sub_bare, out)
    return out, done, sorted(unknown)
