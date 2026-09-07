# Observation audit: what was in the submitted version, what was added in this revision

Sources
- **ORIGINAL** = `Achard et al_submission_Jneuro.pdf` (the version submitted to J Neurosci),
  Tables 1, 2, S1–S4.
- **CURRENT** = `manuscript/latex/delay.tex`, Tables 3–6 and `tab:pp`.

The colour marking in `Achard_et_al_R1_v3_tracked-changes.docx` (blue = new) tracks
**rewritten text**, not new findings: every F-statistic in the Results is blue because
the section was rewritten. It cannot be used for this audit. The submitted PDF can.

---

## 1. Behaviour — unchanged

| Measure | Effects (ORIGINAL = CURRENT) |
|---|---|
| Accuracy | duration, coherence, duration × coherence |
| Mean RT | duration, coherence, duration × coherence |
| Mean PMET | duration, coherence, duration × coherence |
| Mean MET | none (a reported null) |

No change.

---

## 2. PMA rate, by window

`L` = level effect, `T` = time/ramp effect. ✓ significant, – not.

### Stimulus period (0–2 s from stimulus onset, 20 bins)

| Effect | ORIGINAL | CURRENT | status |
|---|---|---|---|
| duration (L) | 13.38, p=.001 ✓ | 13.381, p=.001 ✓ | same |
| coherence (L) | 5.25, p=.021 ✓ | 5.251, p=.021 ✓ | same |
| time (T) | 18.05, p<.001 ✓ *omnibus* | 30.869, p<.001 ✓ *linear trend* | **statistic changed** |
| duration × coherence | 0.312, p=.734 – | 0.312, p=.731 – | same |
| duration × time | 4.645, p<.001 ✓ *omnibus* | 9.532, p=.006 ✓ *linear* | **statistic changed** |
| coherence × time | 1.718, p=.005 ✓ *omnibus* | 6.145, p=.007 ✓ *linear* | **statistic changed** |
| 3-way | 1.058, p=.377 – | 0.729, p=.480 – | same |
| **time, quadratic** | — | 2.621, p=.120 – | **ADDED** |
| **dur × time, quadratic** | — | 3.555, p=.073 – | **ADDED** |
| **coh × time, quadratic** | — | 0.525, p=.586 – | **ADDED** |
| **3-way, quadratic** | — | 0.403, p=.625 – | **ADDED** |

### Delay period (1–3 s from stimulus offset, 20 bins)

| Effect | ORIGINAL | CURRENT | status |
|---|---|---|---|
| duration (L) | 33.19, p<.001 ✓ | 33.186, p<.001 ✓ | same |
| coherence (L) | 3.36, p=.044 ✓ | 3.361, p=.045 ✓ | same |
| duration × coherence (L) | 3.474, p=.040 ✓ | 3.474, p=.042 ✓ | same |
| time (T) | 31.42, p<.001 ✓ *omnibus* | 35.238, p<.001 ✓ *linear* | **statistic changed** |
| duration × time | 4.520, p<.001 ✓ *omnibus* | 11.898, p=.002 ✓ *linear* | **statistic changed** |
| coherence × time | 0.737, p=.879 – | 0.524, p=.567 – | same (null) |
| 3-way | 0.959, p=.542 – | 0.508, p=.580 – | same (null) |
| **time, quadratic** | — | 18.891, p<.001 ✓ | **ADDED, and significant** |
| **dur × time, quadratic** | — | 0.751, p=.396 – | **ADDED** |
| **coh × time, quadratic** | — | 1.855, p=.170 – | **ADDED** |
| **3-way, quadratic** | — | 0.946, p=.391 – | **ADDED** |
| **log-rate slope (point process)** | — | +0.929, F=102.95, p<.001 | **ADDED** |
| **slope × duration** | — | F=4.78, p=.040 | **ADDED** |
| **slope × coherence** | — | F=4.20, p=.038 | **ADDED** |
| **slope × dur × coh** | — | F=0.39, p=.663 | **ADDED** |

### RS period (0–0.2 s from RS onset, 2 bins)

| Effect | ORIGINAL | CURRENT | status |
|---|---|---|---|
| duration (L) | 36.56, p<.001 ✓ | 36.559, p<.001 ✓ | same |
| coherence (L) | 1.159, p=.315 – | 1.159, p=.315 – | same |
| time (T) | 21.96, p<.001 ✓ | 21.964, p<.001 ✓ | same |
| duration × time | 7.046, p=.015 ✓ | 7.046, p=.015 ✓ | same |
| others | – | – | same |

No change (with two bins the omnibus and the linear trend coincide).

### Response period (−0.2–0 s, 2 bins in ORIGINAL, 1 bin in CURRENT)

| Effect | ORIGINAL (2 bins) | CURRENT (1 bin) | status |
|---|---|---|---|
| duration (L) | 0.011, p=.919 – | 7.613, p=.012 ✓ | **REVERSED by the re-analysis** |
| coherence (L) | 0.724, p=.491 – | 1.845, p=.173 – | same |
| duration × coherence | 0.636, p=.534 – | 0.342, p=.697 – | same |
| time (T) | 38.31, p<.001 ✓ | — | **REMOVED** |
| duration × time | 5.013, p=.036 ✓ | — | **REMOVED** |

Reason for the change: the final 0.1 s bin is structurally zero for every participant and
condition, because a PMA is by definition an activation preceding the response-producing
burst. The two-bin analysis therefore compared an informative bin with a constant. This is
a methodological correction, not a new observation, but it **flips a duration effect from
null to significant** and **removes two significant effects**.

---

## 3. MBL and LRP — unchanged

| Measure | Window | ORIGINAL | CURRENT |
|---|---|---|---|
| MBL | stimulus | time ✓, coherence × time ✓ | same |
| MBL | delay | time ✓ | same |
| MBL | response | time ✓ | same |
| LRP | delay | time ✓ | same |
| LRP | RS | time ✓ | same |
| LRP | response | time ✓ | same |

---

## 4. Summary: what this revision added

1. **Quadratic (curvature) contrasts** in every window — 8 new rows, one of them
   significant (delay, p<.001).
2. **The point-process / log-rate analysis** of the delay period — 4 new rows, three
   significant. This is the whole of `tab:pp`, and with it the log-rate framing,
   the Kac–Rice appendix, and the S2/S3 predictions.
3. **The response-period re-analysis** on one bin — a correction, but it changes three
   entries.
4. **Trend contrasts replacing omnibus time effects** — same effects, different statistics.

Nothing else changed. The original observation set is 22 significant effects plus two
reported nulls (mean MET; the delay coherence × time interaction).
