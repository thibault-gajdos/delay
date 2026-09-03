#!/usr/bin/env python3
"""
manuscript/edits/2026-09-03_nine_effects_rationale.py

Point 3 of the 2026-09-03 review pass: "The effects reported are the contrasts that
carry an argument rather than the main effects, since a main effect that both
architectures produce with the same sign at every parameter setting adds nothing to the
derivation that establishes it" was flagged as incomprehensible, and the whole paragraph
asked to be redone.

Read closely against the rest of the manuscript, the sentence is not just badly worded --
it states a rule the paper does not actually follow. Three of the nine reported effects
(the duration benefit on accuracy, the duration effect on response time, and its
coherence interaction) ARE main effects that both architectures produce with the same
sign at every setting -- Behavior > Theoretical predictions says so explicitly ("Every
behavioral prediction above therefore holds in both models, with the same sign and the
same ordering") -- and they are reported anyway. So is a fourth (the coherence effect on
the stimulus-period PMA rate), which Stimulus-locked > Simulation validation calls "a
control" for exactly that reason. And the PMET/MET duration dissociation is ALSO a shared
prediction of both models (RS-locked > Theoretical predictions: "Both models predict...
this effect should be evident only in... PMET... whereas... MET should remain
unchanged") -- not a discriminator, but the claim that gives the response-time
decomposition its meaning.

scripts/effect_panel.py's own docstring makes the real, more nuanced rule explicit
("kept because their SIZE is the thing a derivation cannot give" for the RT pair; "a
control" for the stimulus-period coherence effect; the dissociation and B.2 named
separately) -- this rewrite adapts that reasoning into the manuscript's own prose,
correctly attributing each of the nine to the reason it is actually there, rather than
repeating the single (and, on five of the nine, false) "main effects are excluded" rule.
Sourced from: Behavior > Theoretical predictions and > Simulation validation (body.tex),
Stimulus-locked > Theoretical predictions and > Simulation validation, RS-locked >
Theoretical predictions, and scripts/effect_panel.py's docstring -- not invented.

Splice point: paragraph 112 (see 2026-09-03_effect_size_comparison.py for the paragraph's
other edit), runs 2-4, which together hold the confusing sentence and the itemised list
that follows it, up to and including "...whatever the sign of its prediction." Run 5
("Figure S5 shows all nine.") is untouched. Runs 2-4 are replaced by one run (their
combined content is superseded, not appended to), matching their own formatting (Normal,
plain, color 1F4E9E -- verified programmatically, not assumed, same as point 1 and 2).

Usage: python3 rewrite_nine_effects.py <path-to-docx>   (edits in place)
"""
import sys
import docx
from docx.oxml.ns import qn

RUN2_OLD_START = "The effects reported are the contrasts that carry an argument"
RUN4_OLD_END = "quantitatively inadequate on that effect whatever the sign of its prediction."

NEW_TEXT = (
    " The nine effects reported are the ones this paper's argument actually turns on, "
    "not an exhaustive catalogue of what six-cell contrasts could compute. A derivation "
    "gives a prediction's sign and shape, never its size, so an effect is worth reporting "
    "whenever the observed value's fit inside an architecture's producible range is not "
    "already settled --- whether or not the sign itself is. Three of the nine have a sign "
    "already settled, identically, for both architectures (Behavior, Theoretical "
    "predictions, above): the duration benefit on accuracy, the duration effect on "
    "response time, and its coherence interaction; kept because a shared sign says "
    "nothing about whether the observed SIZE is reached. A fourth, the coherence effect "
    "on the stimulus-period partial-activation rate, both models predict alike too, and "
    "is reported as a control, showing the two architectures describe that window "
    "equally well outside the one contrast on which they diverge (Stimulus-locked, "
    "above). Two more, the duration effects on premotor and on motor execution time, are "
    "likewise shared predictions of both models, reported together because together they "
    "are the claim that gives the response-time decomposition its meaning. The remaining "
    "three are where an architecture's sign is genuinely at stake: the duration effect on "
    "the stimulus-period rate, the one prediction outside the delay period on which the "
    "two models disagree; the coherence × duration interaction on the delay-period "
    "rate, the paper's acknowledged disagreement with the data; and the delay-period "
    "log-rate slope, on which the two architectures separate outright. A model whose "
    "range excludes the observed value is quantitatively inadequate on that effect "
    "whatever the sign of its prediction. "
)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        "Achard_et_al_R1_v3_tracked-changes_2026-09-02.docx"
    d = docx.Document(path)

    hits = [i for i, p in enumerate(d.paragraphs) if RUN2_OLD_START in p.text]
    if len(hits) != 1:
        raise SystemExit(f"expected exactly one anchor paragraph, found {len(hits)}: {hits}")
    idx = hits[0]
    p = d.paragraphs[idx]

    runs = p.runs
    i2 = next(j for j, r in enumerate(runs) if RUN2_OLD_START in r.text)
    i4 = next(j for j, r in enumerate(runs) if r.text.rstrip().endswith(RUN4_OLD_END))
    assert i4 == i2 + 2, (i2, i4, [r.text[:40] for r in runs])
    r2, r3, r4 = runs[i2], runs[i2 + 1], runs[i2 + 2]

    def fmt(r):
        rPr = r._r.find(qn("w:rPr"))
        b = rPr.find(qn("w:b"))
        it = rPr.find(qn("w:i"))
        bold = b is not None and b.get(qn("w:val")) not in ("0", "false")
        italic = it is not None and it.get(qn("w:val")) not in ("0", "false")
        color = r.font.color.rgb if r.font.color and r.font.color.type else None
        return bold, italic, str(color) if color else None

    f2, f3, f4 = fmt(r2), fmt(r3), fmt(r4)
    assert f2 == f3 == f4, (f2, f3, f4)
    assert f2 == (False, False, "1F4E9E"), f2

    # sanity: runs 2..4 concatenated should equal what we think we're replacing
    old_concat = r2.text + r3.text + r4.text
    assert old_concat.startswith(" " + RUN2_OLD_START) or old_concat.startswith(RUN2_OLD_START)
    assert old_concat.rstrip().endswith(RUN4_OLD_END)

    r2.text = NEW_TEXT
    # remove runs 3 and 4 from the XML (their content is now folded into r2)
    r3._r.getparent().remove(r3._r)
    r4._r.getparent().remove(r4._r)

    d.save(path)

    # re-open and verify
    d2 = docx.Document(path)
    found = [p2 for p2 in d2.paragraphs if "The nine effects reported are the ones" in p2.text]
    assert len(found) == 1
    p2 = found[0]
    for must in [
        "duration benefit on accuracy", "coherence interaction",
        "reported as a control", "gives the response-time decomposition its meaning",
        "delay-period log-rate slope, on which the two architectures separate outright",
        "A model whose range excludes the observed value is quantitatively inadequate",
        "Figure S5 shows all nine.",
    ]:
        assert must in p2.text, f"missing: {must}"
    assert "carry an argument rather than the main effects" not in p2.text
    print("OK: paragraph rewritten and verified, length", len(p2.text))


if __name__ == "__main__":
    main()
