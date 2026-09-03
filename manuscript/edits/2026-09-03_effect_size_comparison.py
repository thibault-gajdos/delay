#!/usr/bin/env python3
"""
manuscript/edits/2026-09-03_effect_size_comparison.py

Point 2 of the 2026-09-03 review pass: the comparison between a model's across-draw
effect-size range and the observed, between-participant value was stated mechanically
("a draw standing for one hypothetical participant") but never argued -- in particular,
nothing said why raising the simulated trial count per draw (4000/cell, far above a real
participant's ~156/cell) does not just manufacture spurious precision, the way pooling
4000-trial "replicates" of one fixed parameter setting would. It doesn't, because
admissible draws differ from each other in their SAMPLED PARAMETERS (Table A1), never in
trial count -- trial count only sharpens one draw's own estimate of its own value. This
inserts two sentences making that argument explicit, anchored to the one place in the
manuscript where it is checked numerically (the delay-period log-rate slope's
Monte-Carlo/between-draw variance decomposition, already in Results > Delay-locked,
sourced from scripts/pool_noise.py / results/noise_pooled.json), and stating plainly that
the scope of that check is one effect, not independently re-verified for the other eight
in Figure S5.

Splice point: inside the single run of paragraph 112 (the "What the evaluation reports is
a range of effect sizes..." paragraph) that already runs from "says what the architecture
can produce" to "...bound of the Appendix is drawn against." -- i.e. python-docx's
`.runs` boundaries here do NOT match the \\new{\\textbf{\\emph{...}}} chunking docx2tex.py
prints to LaTeX (that chunking splits one Word run into several \\new{} groups for its own
reasons); the actual XML has ONE run covering that whole span, formatted b=0 i=0 color
1F4E9E (Normal, plain, blue -- same convention as the point-1 paragraph, see that
edit script's note on docx2tex.py's bold-detection bug). Because formatting is uniform
across the whole run, the new sentences are spliced into that run's own text (no XML
surgery needed) rather than inserted as a separate run.

Usage: python3 add_effect_size_comparison.py <path-to-docx>   (edits in place)
"""
import sys
import docx
from docx.oxml.ns import qn

ANCHOR_BEFORE = "a draw standing for one hypothetical participant."
ANCHOR_AFTER = " The order matters for the delay-period log-rate slope"

NEW_TEXT = (
    " Trial count and parameter setting play different roles here: a draw's 4000 trials "
    "per cell sharpen the precision with which that setting's own value is estimated; "
    "they do not create additional hypothetical participants. Admissible draws differ "
    "from one another only in their sampled parameters, never in trial count, so the "
    "across-draw spread reported below is expected to reflect the parameter space rather "
    "than simulation noise, provided the latter is small next to it. This was checked "
    "directly for the delay-period log-rate slope (Results, Delay-locked): re-simulating "
    "admissible draws on several seeds and decomposing the variance attributes 19% of "
    "Model 1's total spread and 0.1% of Model 2's to Monte-Carlo noise, the rest "
    "genuinely between draws; the same trial count is used for every effect, though the "
    "decomposition was not repeated for each. Because the comparison below is with this "
    "across-draw range and not with a pooled point estimate, it asks whether any "
    "admissible setting reaches the observed value, not whether the sampled parameters "
    "reproduce the population distribution of individual differences --- which is also "
    "why an interval, not a mean and a standard error, is what is reported."
)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        "Achard_et_al_R1_v3_tracked-changes_2026-09-02.docx"
    d = docx.Document(path)

    hits = []
    for i, p in enumerate(d.paragraphs):
        if ANCHOR_BEFORE in p.text and ANCHOR_AFTER in p.text:
            hits.append(i)
    if len(hits) != 1:
        raise SystemExit(f"expected exactly one anchor paragraph, found {len(hits)}: {hits}")
    idx = hits[0]
    p = d.paragraphs[idx]

    target_runs = [r for r in p.runs if ANCHOR_BEFORE in r.text and ANCHOR_AFTER in r.text]
    if len(target_runs) != 1:
        raise SystemExit(f"expected exactly one run holding both anchors, found {len(target_runs)}")
    r = target_runs[0]

    rPr = r._r.find(qn("w:rPr"))
    b_el = rPr.find(qn("w:b"))
    i_el = rPr.find(qn("w:i"))
    bold = b_el is not None and b_el.get(qn("w:val")) not in ("0", "false")
    italic = i_el is not None and i_el.get(qn("w:val")) not in ("0", "false")
    color = r.font.color.rgb if r.font.color and r.font.color.type else None
    assert not bold and not italic, (bold, italic)
    assert str(color) == "1F4E9E", color

    old_text = r.text
    split_at = old_text.index(ANCHOR_BEFORE) + len(ANCHOR_BEFORE)
    assert old_text[split_at:split_at + len(ANCHOR_AFTER)] == ANCHOR_AFTER
    new_text = old_text[:split_at] + NEW_TEXT + old_text[split_at:]
    r.text = new_text

    d.save(path)

    # re-open and verify
    d2 = docx.Document(path)
    found = False
    for p2 in d2.paragraphs:
        if "Trial count and parameter setting play different roles here" in p2.text:
            found = True
            assert "population distribution of individual differences" in p2.text
            assert ANCHOR_BEFORE in p2.text and ANCHOR_AFTER.strip() in p2.text
            idx_before = p2.text.index(ANCHOR_BEFORE)
            idx_new = p2.text.index("Trial count and parameter setting")
            idx_after = p2.text.index(ANCHOR_AFTER.strip())
            assert idx_before < idx_new < idx_after, (idx_before, idx_new, idx_after)
    assert found
    print("OK: new sentences spliced in, ordering verified")


if __name__ == "__main__":
    main()
