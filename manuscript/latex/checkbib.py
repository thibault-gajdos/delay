# -*- coding: utf-8 -*-
"""Check the manuscript's reference list against the Zotero BibTeX export.

The references in the Word copy are plain APA text, not \\cite keys, so nothing verifies
them: a reference can be in the paper and absent from Zotero, or in both with different
years, and neither Word nor LaTeX will say a word.

Two independent tests are run, because each fails differently. Matching by DOI is exact
but misses entries whose Zotero record has no DOI field -- 320 of the 1112 here. Matching
by first-author surname plus title words catches those, but can be fooled by two papers by
the same author on the same subject, so it demands that most content words agree. A
reference counts as present if either test finds it.

  python3 checkbib.py [body.tex] [bib file]

Writes zotero_missing_dois.txt: the DOIs of the references that are in neither, ready to
paste into Zotero's "Add Item by Identifier" box, which fetches the records itself rather
than trusting anything reconstructed from the APA text here.
"""
import re, sys, os, unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))
BODY = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'body.tex')
BIB = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser(
    '~/mnt/tools/bib_thib_zotero.bib')


def deaccent(s):
    return ''.join(c for c in unicodedata.normalize('NFKD', s)
                   if not unicodedata.combining(c))


def plain(s):
    s = re.sub(r'\\(url|emph|textbf|new|old|del)\{', '{', s)
    s = re.sub(r'\\[a-zA-Z]+\s*', ' ', s)
    s = s.replace('\\&', '&').replace('\\_', '_').replace('\\%', '%')
    s = re.sub(r'[{}]', '', s)
    return re.sub(r'\s+', ' ', s).strip()


def words(s):
    return set(re.findall(r'[a-z]{4,}', deaccent(s).lower()))


# ---------------------------------------------------------------- the manuscript
body = open(BODY, encoding='utf-8').read()
m = re.search(r'\\begin\{references\}(.*?)\\end\{references\}', body, re.S)
if not m:
    sys.exit('no references environment in ' + BODY)
refs = []
for line in m.group(1).split('\n'):
    if len(line.strip()) < 20:
        continue
    text = plain(line)
    doi = re.search(r'doi\.org/(\S+?)\.?$', text)
    refs.append({'text': text,
                 'struck': line.lstrip().startswith(r'\del{'),
                 'doi': doi.group(1).lower() if doi else None})

# ---------------------------------------------------------------- the .bib
entries = []
for chunk in re.split(r'\n@', '\n' + open(BIB, encoding='utf-8', errors='replace').read()):
    key = re.match(r'\w+\{([^,]+),', chunk)
    if not key:
        continue

    def field(name):
        """One field's value. Braces nest -- Zotero writes {{Using}} to protect case --
        so this scans with a brace counter; a non-greedy regex would stop at the first
        inner "}". Numeric fields (year = 2024) carry no braces at all."""
        f = re.search(r'\n\s*%s\s*=\s*' % name, chunk, re.I)
        if not f:
            return ''
        i = f.end()
        if i >= len(chunk):
            return ''
        if chunk[i] not in '{"':
            return plain(chunk[i:].split(',')[0].split('\n')[0])
        depth, j = 0, i
        while j < len(chunk):
            if chunk[j] == '{':
                depth += 1
            elif chunk[j] == '}':
                depth -= 1
                if depth == 0:
                    break
            j += 1
        return plain(chunk[i + 1:j])
    entries.append({'key': key.group(1).strip(), 'author': field('author'),
                    'title': field('title'), 'year': field('year') or field('date'),
                    'doi': field('doi').lower()})


def surname_of(author):
    """Surname of the first author: BibTeX writes "Last, First" or "First Last"."""
    first = author.split(' and ')[0].strip()
    if not first:
        return ''
    if ',' in first:
        return first.split(',')[0].strip()
    parts = first.split()
    return parts[-1] if parts else ''


by_doi = {e['doi']: e for e in entries if e['doi']}
by_surname = {}
for e in entries:
    s = deaccent(surname_of(e['author'])).lower().strip()
    by_surname.setdefault(s, []).append(e)
    if ' ' in s:                       # "van Ede" is also findable under "ede"
        by_surname.setdefault(s.split()[-1], []).append(e)

# ---------------------------------------------------------------- match
for r in refs:
    # APA gives "van Ede, F., & Nobre" -- the surname may be several words, so it is
    # everything up to the first comma, not the first word.
    sur = re.match(r"([A-Za-zÀ-ÿ''\- ]+),", r['text'])
    r['sur'] = deaccent(sur.group(1)).lower().strip() if sur else ''
    yr = re.search(r'\((\d{4})[a-z]?\)', r['text'])
    r['year'] = yr.group(1) if yr else ''
    title = r['text'][r['text'].find(').') + 2:][:170] if ').' in r['text'] else r['text']
    r['by_doi'] = by_doi.get(r['doi']) if r['doi'] else None
    r['by_name'] = None
    wt = words(title)
    for e in by_surname.get(r['sur'], []) or by_surname.get(
            r['sur'].split()[-1] if r['sur'] else '', []):
        we = words(e['title'])
        if not wt or not we:
            continue
        common = len(wt & we)
        # A title is the same title only when most of its content words agree; a bare
        # count matched "Estimating parameters of the diffusion model" to "Reexamining
        # the effects of speed" on three ordinary words.
        if common >= 3 and common / min(len(wt), len(we)) >= 0.6:
            r['by_name'] = e
            break

live = [r for r in refs if not r['struck']]
struck = [r for r in refs if r['struck']]
found = [r for r in live if r['by_doi'] or r['by_name']]
absent = [r for r in live if not (r['by_doi'] or r['by_name'])]
partial = [r for r in found if not r['by_doi'] and r['doi']]
noyear = [r for r in found if (r['by_name'] and not r['by_name']['year'])]

print('manuscript: %d references (%d live, %d struck)' % (len(refs), len(live), len(struck)))
print('Zotero export: %d entries, %d with a DOI\n' % (len(entries), len(by_doi)))
print('IN ZOTERO                    %d' % len(found))
print('NOT IN ZOTERO                %d' % len(absent))
print('  in Zotero but that record has no DOI  %d' % len(partial))
print('  in Zotero but that record has no year %d' % len(noyear))

print('\n--- not in Zotero ---')
for r in absent:
    print('  %-44s %s' % (r['doi'] or '(no DOI in the paper)', r['text'][:70]))
print('\n--- in Zotero, record incomplete ---')
for r in partial:
    print('  no DOI on %-42s %s' % (r['by_name']['key'][:42], r['text'][:52]))
for r in noyear:
    print('  no year on %-41s %s' % (r['by_name']['key'][:41], r['text'][:52]))
print('\n--- struck (proposed deletions) ---')
for r in struck:
    print('  %-14s %s' % ('in Zotero' if (r['by_doi'] or r['by_name']) else 'absent',
                          r['text'][:70]))

out = os.path.join(os.path.dirname(os.path.abspath(BODY)), 'zotero_missing_dois.txt')
dois = [r['doi'] for r in absent if r['doi']]
open(out, 'w').write('\n'.join(dois) + '\n')
print('\n%d DOIs written to %s' % (len(dois), out))
