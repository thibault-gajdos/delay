Achard, Gajdos Preuss & Servant — "Graded Motor-Preparatory Dynamics Preserve Perceptual Decisions over Short Delays" (Journal of Neuroscience, R1)

DRAFT — internal working copy, not yet formatted for submission.

## Note pour les auteurs — à lire avant envoi

Ce document est un brouillon complet de la lettre de réponse, point par point. Chaque réponse ci-dessous est rédigée comme si elle allait être envoyée aux rapporteurs. Il a été remis à jour le 8 septembre 2026 contre le texte courant du manuscrit (`manuscript/latex/delay.tex`), qui en est la seule source. Les points suivants ne sont PAS pleinement traités et méritent une décision avant l'envoi — ils sont aussi repérés en couleur directement à l'endroit concerné dans la lettre (« NOTE INTERNE »).

-   **C2.4 (renommer « second accumulation phase ») — volontairement NON traité**, sur votre décision explicite : vous jugez le rapporteur dans l'erreur sur ce point (il faut bien accumuler le signal de réponse perceptif, littérature à l'appui). La lettre de réponse ci-dessous défend donc la formulation actuelle plutôt que de la changer — à relire pour vérifier que l'argumentaire vous convient avant envoi.

-   **C2.7 (autres formats sensoriels) — traitement partiel** : le volet « animal work » est couvert (Cisek & Kalaska, 2005, décisions d'atteinte chez le singe). Le volet « autres formats sensoriels » (accumulation de preuves auditives, tactiles, etc.) n'est PAS traité — aucune référence identifiée dans le temps imparti.

-   **C1.7 (détail « α fonction du délai »)** — la demande principale (montrer le build-up pour différents délais) est traitée par la citation à Achard et al. (2025). Le point plus fin de la note — que cela fait du taux de préparation une fonction de la longueur du délai, à énoncer explicitement comme condition de portée — n'est pas formulé dans ces termes précis.

-   **Référence cassée dans la bibliographie** — l'entrée `MotorCorticalSignals` (NeuroImage 2022, doi 10.1016/j.neuroimage.2022.119667) est un import Zotero depuis « Elsevier Enhanced Reader » : ni auteurs ni année. Elle est citée six fois et s'affiche dans le PDF sous la forme `("Motor Cortical Signals Reflecting Decision Making and Action Preparation | Elsevier Enhanced Reader", n.d.)` — ce sont les sept dernières occurrences de « n.d. » du manuscrit. À réparer depuis votre Zotero avant tout envoi.

-   **Archive OSF** — vérifier que les exports corrigés du 19/08 (`PMA_rate_resp_JASP.csv` et `PMA_rate_rs_JASP.csv`) y sont bien en place : la réponse « Volunteered 3 » l'affirme au rapporteur.

### Ce qui a changé depuis la version du 5 septembre

-   **Le démarrage stationnaire.** La fluctuation de maintien part maintenant de son amplitude stationnaire, ε(τ) ~ N(0, SD_ε²), au lieu de zéro, et H(τ) = 0 est imposé par construction (la sortie du filtre est continue en τ) plutôt que supposé. Toutes les affirmations sur la période de délai deviennent de ce fait **exactes** au lieu d'asymptotiques : sous Model 1 le taux est rigoureusement constant sur la fenêtre, à tout κ, tout λ et tout niveau maintenu (Corollaire A.2). L'ancienne « frontière » — qui bornait la pente de Model 1 par un argument d'optimisation — a disparu du manuscrit et du code, remplacée par cet énoncé exact. Les réponses C1.8, C2.1 et C2.3 sont réécrites en conséquence.

-   **Statistiques sur le taux brut, non sur son logarithme.** Le manuscrit ne teste plus la pente du log-taux ni l'analyse de processus ponctuel non binnée : LEVEL, TREND et QUAD sont des contrastes polynomiaux orthogonaux sur le taux lui-même, dans les données comme dans les simulations. Le logarithme est indéfini dans 6 % des cellules participant × condition × bin de la fenêtre de délai, 30 % de celles de la fenêtre stimulus et 55 % de celles de la fenêtre RS. Rien n'est perdu : en divisant les contrastes du taux par son niveau on retrouve ν̇/ν et ν̈/ν, c'est-à-dire exactement ce que les contrastes sur le logarithme estimeraient. **Il n'y a donc plus de Table S6** (le complément de `_estim_body.tex` sur les trois estimandes tombe avec elle), et le paragraphe de C1.1 qui expliquait un désaccord apparent entre analyse binnée et non binnée est supprimé : il n'y a plus deux analyses.

-   **La courbure.** Le contraste quadratique est maintenant une prédiction à part entière : nulle sous Model 1, positive sous Model 2 dès que le niveau maintenu est à plus d'un écart-type de fluctuation de la porte. Observée à +0,210 s⁻³ (F(1,21) = 18,891), sans modulation par la cohérence ni par la durée. C'est un second test dans la même fenêtre, sur une statistique que Model 1 fixe elle aussi à zéro. Ajouté à C1.8.

-   **L'inversion (a, b) est retirée.** Le passage des Résultats qui remontait de niveau/pente/courbure aux deux quantités de la théorie (a = 1,38, b = 0,57 s⁻¹) a été supprimé, ainsi que le paragraphe correspondant de la Discussion : joli, mais il compliquait le message sans le renforcer. La lettre ne le mentionne plus.

-   **Chiffres à jour.** Tirages admissibles : **160** pour Model 1 et **145** pour Model 2 (et non 164 / 148). Effet du coherence × duration sur le taux de délai : sens observé à **3 %** des tirages de Model 1 et **5 %** de ceux de Model 2, sens opposé à **88 %** et **86 %** (et non 4 / 7 / 85 / 87). Effet de durée sur le niveau en fenêtre stimulus : Model 1 se trompe de signe à **96 %** de ses tirages (et non 95 %). Toutes les valeurs de la lettre viennent désormais de la campagne « chaude » (démarrage stationnaire) ; aucune de l'ancienne campagne froide.

-   **Un effet discriminant de plus.** La duration-dépendance de la rampe en fenêtre stimulus : Model 1 la produit à l'envers à 94 % de ses tirages et jamais dans le bon sens, Model 2 dans le bon sens à 37 %. La fenêtre stimulus exclut donc Model 1 deux fois, sur le niveau et sur la rampe. Ajouté à C2.1.

-   **Fuite : vérification rejouée.** Les nombres cités en C2.3 sont ceux de la campagne chaude (160 → 155 → 152 → 138 pour Model 1 ; 145 → 143 → 142 → 130 pour Model 2 ; « les six pentes positives » : 78 → 76 → 71 → 77 pour Model 2, **0 à tous les taux** pour Model 1). Contrairement à ce qu'indiquait la note du 5 septembre, l'argument **est** maintenant dans l'article : Remarque A.6 de l'annexe le démontre structurellement, et un paragraphe des Résultats le résume.

-   **Correction de la Figure 7 (8 septembre).** `scripts/mkfig.py` fixait $\sigma_\varepsilon = 1$ au lieu d'utiliser le SD_ε de chaque tirage, ce qui gonflait $\sigma_H$ d'un à deux ordres de grandeur, réduisait d'autant la distance standardisée à la porte, et sous-estimait l'erreur de l'approximation stationnaire de quatre ordres de grandeur. Corrigé : l'écart maximal sur la fenêtre est **6,8 × 10⁻³ %** (et non 1,4 × 10⁻³ %) sur les tirages admissibles, qui sont désormais les seuls tracés. La conclusion ne bouge pas — c'est toujours négligeable — mais la nouvelle **Corollaire A.2** dit pourquoi il fallait le calculer : l'erreur sur les moments est amplifiée par le carré de la distance standardisée à la porte, parce que celle-ci entre dans le taux par une gaussienne. Sur la boîte entière (tirages non admissibles compris) l'écart atteint le pour cent.

-   **Figures et annexe.** Sept figures principales : 1 dynamiques simulées, 2 design + essai unique, 3 comportement, 4 taux de PMA et EEG, 5 forme de la courbe de taux pendant le délai, 6 tailles d'effet (22 panneaux, groupés par fenêtre), 7 erreur de l'approximation stationnaire. L'annexe unique va de la **Proposition A.1** à la **Remarque A.7** (A.1 régime stationnaire, A.2 erreur de l'approximation stationnaire sur le taux, A.3 platitude sous Model 1, A.4 rampe et parabole sous Model 2, A.5 effet de durée en fenêtre stimulus, A.6 fuite, A.7 portée des approximations). Les tables supplémentaires sont S1 (comportement), S2 (taux de PMA), S3 (LRP), S4 (MBL), S5 (MBL, moitiés du stimulus).

Tout le reste des commentaires C1.1–C1.9, C2.1–C2.7, la note éditoriale et les points volontaires est traité intégralement dans le manuscrit, vérifié passage par passage contre le texte courant de `delay.tex`.

# Reviewer 1

## C1.1 — Slopes rather than bins.

> "Although there is a nice cartoon figure for the hypotheses, the authors use a repeated measures ANOVA for their tests. I would prefer tests on the different slopes in the EMG and EEG instead of the arbitrary 0.1s bins. [...] Ln 324: 'slope of this ramping activity' → that is not what you ANOVA tested."

**Response.** We agree that the earlier phrasing overstated what an omnibus ANOVA on a nineteen-degree-of-freedom time-bin factor establishes, and have gone beyond correcting the language. We now decompose the time-bin factor into orthogonal polynomial contrasts — the standard within-subjects contrast analysis — and report those as the primary statistics in the Results text for every measure and period where the decomposition applies (stimulus- and delay-period bins for partial-activation rate, mu–beta lateralization and the lateralized readiness potential); Tables S1–S4 report them as primary as well. The linear contrast is reported throughout as the ordinary-least-squares slope against time, in units of the measure per second, which is the same contrast rescaled so that the coefficient reads as a rate of change. Where a period has at least three bins the quadratic contrast is reported alongside it, for two reasons: it is where the two candidate models make opposite and exact predictions about the delay period (see C1.8), and mu–beta lateralization during the stimulus period is non-monotone, so a slope alone would misclassify it as null. Each contrast is tested against its own error term rather than against an error term pooled across components, which is what makes a sphericity correction unnecessary for the time-bin component; Greenhouse–Geisser correction is applied throughout to the three-level coherence factor. No signal was reprocessed: the contrasts are weighted sums of the bin means already computed. We should add that the omnibus time-bin test the comment objects to is not merely uninformative about slope but poorly behaved here in its own right — its Greenhouse–Geisser epsilon is 0.14 on 19 degrees of freedom. The analysis script that produces the tables is deposited with the data. Two periods need no substitution: the response-signal period has only two time bins, so the omnibus and linear-trend tests coincide exactly, and the response period is analyzed on a single bin, so no trend decomposition applies there.

For the partial-activation rate, the contrasts are taken on the rate itself and not on its logarithm, and this deserves a word, since the theory is naturally written on the logarithm — where the delay-period prediction is exactly a downward parabola of constant curvature (Appendix, Corollary A.4). The logarithm is undefined wherever a cell contains no activation, which is the case in 6% of the participant × condition × bin cells of the delay period, 30% of those of the stimulus period and 55% of those of the response-signal period; deleting or offsetting those cells would change the estimand by an amount that depends on the rate, which is the quantity under test. Nothing is lost by working on the rate: dividing its linear and quadratic contrasts by its level returns the relative slope and relative curvature that contrasts on the logarithm estimate, and the Results use exactly those two ratios to recover the two quantities of the delay-period theory without fitting anything (see C1.8).

Two consequences of the change of criterion are worth stating explicitly rather than leaving to be discovered, and they run in opposite directions. Mu–beta lateralization's stimulus-period coherence × time effect is significant on the omnibus test but null as a full-window linear trend (p = .596), because the measure rises over roughly the first second and then flattens; Results now says so, and an a priori first-half/second-half split (Table S5) locates the coherence-dependent trend in the first half. Conversely, mu–beta lateralization's stimulus-period duration × time effect was not significant on the omnibus test (p = .136) but is significant as a linear trend (p = .046) — a newly significant effect created by the change of criterion. Both are stated in the Results text and reported in full in Table S4.

Finally, one wording consequence for the Abstract. The delay-period claim is now scale-specific and rests on the contrast the comment asked for: across the delay the partial-activation rate has a strongly positive linear trend (+0.560 s⁻², F(1,21) = 35.238, p < .001) whose coherence-dependence — the diagnostic marker of continued evidence accumulation — is null (p = .567). The Abstract accordingly says that the rate "continued to ramp during the stimulus-free delay, but no longer at a coherence-dependent build-up rate."

*Manuscript location*: Materials and Methods ("Statistical Analysis"); Results, Stimulus-onset-locked and Stimulus-offset-locked subsections (trend and curvature statistics now primary); Tables S1–S4 (the same statistics primary there too); Table S5 (mu–beta lateralization, first-/second-half split); Abstract (delay-period claim).

## C1.2 — Clearer statement of the hypotheses and of what each measure does.

> "I would like a much clearer description of the hypotheses and how the different measures help unraveling. Now the EEG results are a bit of a hotch-potch."

**Response.** Results has been restructured by analysis window rather than by measure: Behavior, then Stimulus-onset-locked, Stimulus-offset-locked, RS-locked, and Response-locked, each split into "Theoretical predictions," "Empirical results" (partial-activation rate, then mu–beta lateralization, then the lateralized readiness potential), and, where available, "Simulation validation." Within every window, the model-specific predictions are stated explicitly before any data are shown, rather than after. For the delay period — the window carrying the paper's central discriminating prediction — "Theoretical predictions" derives what each model implies for the partial-activation rate (Model 1: a rate that is constant across the window at every parameter setting; Model 2: a rising, convex, coherence-independent ramp), the corresponding predictions for the two electroencephalographic measures are stated where those measures are reported, and "Simulation validation" reports the discrimination over the models' admissible parameter space (see C1.8 and C2.1).

*Manuscript location*: Results, the five window subsections, each split into "Theoretical predictions," "Empirical results," and "Simulation validation"; Stimulus-offset-locked > Simulation validation specifically for the model-discriminating claim.

## C1.3 — Muscle-activation onset procedure.

> "Ln 119: please be more explicit about the EMG onset procedure. [...] p. 10, ll. 157-159: PMA should be more precisely defined. Specifically, how is it decided that a muscle activation detected by EMG does not 'culminate in an overt response'?"

**Response.** Materials and Methods now gives the full procedure: candidate bursts were detected at a threshold of seven baseline standard deviations on the rectified signal, with threshold-defined periods separated by less than 25 ms merged into a single candidate burst, an adaptive minimum-sample criterion applied to reject noise-driven detections, and burst boundaries refined with the integrated-profile method (Liu & Liu, 2016; Santello & McDonagh, 1998). We also state the procedure's two consequences for how partial activations are counted: 27% of consecutive activation pairs in the exported event list remain separated by less than 25 ms despite the merge rule (this does not affect reported counts, since the refinement relocates but does not remove interval boundaries), and a partial muscle activation is defined as any detected interval other than the one immediately preceding the response — so 18.5% of trials contain an activation ending within 25 ms of the response-related burst, and in 9.3% of trials that activation would be absorbed into the response burst under strict application of the merge rule. We note that response-period results should be read with this possibility of segmentation in mind. The same detection procedure, with the same merge interval and the same minimum-sample rule, is applied to the simulated traces, so that the model and the muscle are measured alike.

*Manuscript location*: Materials and Methods, "EMG acquisition and preprocessing" (now three paragraphs), and "Evaluation over the parameter space" (the same rule applied to simulated traces).

## C1.4 — Coherence range used by others.

> "Ln 242: make explicit what range of coherences others have used."

**Response.** Added as a scope condition in the Discussion: our three-level coherence range (0%, 12%, 50%) is coarser than the five- or six-level logarithmic series more commonly used to characterize psychometric and chronometric functions in random-dot-motion studies (e.g., 3.2–51.2%: Roitman & Shadlen, 2002; 6.4–51.2%: Kiani et al., 2008); sufficient to establish the coherence-dependent effects reported here, though not their full functional form.

*Manuscript location*: Discussion, scope-conditions paragraph.

## C1.5 — Error bars and individual data.

> "Fig 3: SEM is not a good measure on data constrained between 0 and 1 [...] In modern research it is custom to also show individual participant datapoints."

**Response.** Figure 3 has been regenerated. Each participant's data now appear as faint connected points behind the group means in all four panels. For accuracy specifically — bounded in [0,1], where a symmetric standard error is not appropriate — the interval is now a 95% percentile bootstrap interval computed on the Cousineau–Morey-normalized values (10,000 resamples over participants, Morey-scaled), which is both a within-subjects interval, comparable to the ±1 within-subjects SEM shown for RT, PMET, and MET, and asymmetric where the data approach the bounds. The figure legend describes both intervals.

*Manuscript location*: Figure 3 and its legend.

## C1.6 — Reminder of MET and PMET.

> "Ln 308: please refer back to fig 2 to reiterate what MET and PMET are."

**Response.** A cross-reference to Figure 2 is given at first mention of mean PMET and MET in the Results.

*Manuscript location*: Results, RS-locked > Empirical results, opening paragraph.

## C1.7 — Different delays, to show build-up.

> "Ln 485 [...] If you really want to hammer in this point, you should show build-up for different delays (even within the same trial, where the RS can come at for example 1.5 or 3.0 sec)."

**Response.** Added as a second scope condition in the Discussion, citing a companion study (Achard et al., 2025) that directly manipulates delay length: the same delay-period build-up is preserved across 3–7 s delays when the response mapping is known in advance, but is largely absent when it is revealed only at the response signal.

*Manuscript location*: Discussion, scope-conditions paragraph.

> [NOTE INTERNE — ne pas envoyer telle quelle] Le point plus fin de la note (« cela fait de α une fonction du délai, à énoncer comme condition de portée ») n'est pas formulé dans ces termes précis — voir « Note pour les auteurs ».

## C1.8 — A quantitative model.

> "Ln 505-506: here you hint at a quantitative model, which I believe would strengthen your data better. Now you use ANOVAs, and thus significances to rule out hypotheses that are only expressed in words or cartoons, but are never quantified based on the data."

**Response.** We agree that the paper needed a quantitative model, and it now has one. It works in two stages, and we have tried to be careful about what each of them establishes.

The first is derivation. Every prediction the paper tests is derived from the model equations at the point where it is used, in the Results, with the heavier proofs collected in a single Appendix. That is what replaces the verbal and cartoon arguments the comment objects to: the claims are no longer stated in words and then adjudicated by an ANOVA, they follow from the architecture, and the ANOVAs test them. Three examples. Because the drift of the second accumulation is the same in every condition (Equation 9), coherence and duration can act on response time only through the level maintained at the end of the first stage. Because that level is capped at the internal threshold (Equation 4), the coherence effect on response time is larger after a short stimulus than after a long one. And because temporal preparation is anchored to stimulus offset while the response signal follows it at a fixed interval, stimulus duration cancels from the ramp: it contributes the same amount in all six coherence × duration cells, so no behavioral effect whatever can separate the two models.

The second stage is an evaluation over the models' parameter space, which supplies what a derivation cannot, namely the size of each effect. We do not fit. No parameter is estimated from these data, no parameter value is interpreted, and no comparison between the models rests on a goodness-of-fit statistic. Fitting would be the wrong instrument here for three reasons, which the Discussion states: the models are deliberately reduced, so several parameters are close to unidentifiable from behavior and a fitted value would describe the objective's shape rather than a participant; Model 2 contains Model 1 as the special case β = 0, so a better fit for Model 2 is guaranteed and therefore uninformative; and the predictions the paper tests are ordinal, so nothing about them turns on which setting is optimal. What we do instead is sample the parameter space — 1200 Latin-hypercube draws over the ranges of Table 1, 4000 trials per cell per model per draw — screen those draws against behavioral admissibility criteria stated in Materials and Methods (160 draws admissible for Model 1, 145 for Model 2), and report, for every effect, the range of magnitudes each model can produce over its admissible draws and the fraction of draws at which it produces the observed sign. The observed value and its between-participant interval are read against that range (Figure 6).

The result is that behavior does not separate the two models and the muscle does. Every behavioral effect falls inside both models' ranges, and the two ranges are barely distinguishable from one another — the expected consequence of the derivation above. The delay period separates them, and not by a matter of degree. Because the maintenance fluctuation starts at its stationary amplitude, Model 1's delay-period rate is exactly constant in time: every linear and quadratic contrast of it across the window is zero, at every value of the mean-reversion rate, of the filter gain and of the maintained level (Appendix, Corollary A.3). Over the whole sampled parameter domain the largest departure from flatness the closed-form rate returns is 3 × 10⁻¹¹ s⁻², against an observed linear trend of +0.560 s⁻². The simulations agree: over Model 1's 160 admissible draws the trend has a median of −0.007 s⁻², a 5th-to-95th-percentile range of −0.058 to +0.004, and a largest value of +0.007 — a seventy-fifth of the observation, and itself Monte-Carlo noise, since the Monte-Carlo standard deviation of one draw's trend is 0.0035 s⁻² and re-simulating the seven most promising draws at ten times the trial count leaves the largest at +0.0007 s⁻². Model 2's median is +0.054, its 5th-to-95th-percentile range −0.017 to +0.418, and the observation falls at the 96th percentile of that distribution. Seventy-eight of Model 2's 145 admissible draws produce all six condition trends significantly positive; none of Model 1's 160 do, at any parameter setting, and none at any leak rate either (see C2.3).

The shape of the ramp gives a second, independent test in the same window, and one on which the two models make opposite exact predictions. Under Model 1 a constant rate has zero curvature as well as zero trend. Under Model 2 the rate is convex wherever the maintained level stands more than one fluctuation standard deviation from the gate. The observed quadratic contrast is +0.210 s⁻³ (F(1,21) = 18.891, p < .001), and neither factor modulates it (all p > .17), as the theory requires, since the curvature depends on properties of the participant and not of the condition.


*Manuscript location*: Results, the "Theoretical predictions" subsection of each window (the derivations, stated where they are used); Appendix, Proposition A.1 to Remark A.7 (the proofs); Materials and Methods, "Evaluation over the parameter space"; Results, the "Simulation validation" subsections, Figures 5–7; Discussion, the two paragraphs on what the modeling does and does not establish.

## C1.9 — Statistics in the text.

> "In J Neurosci it is customary to present statistics in the results text instead of everything in a table."

**Response.** Exact test statistics (F, df, p, ηp²) are now reported inline in the Results text throughout, in addition to the summary tables.

*Manuscript location*: Results, throughout.

# Reviewer 2

## C2.1 — Simulate the model; show the predictions hold across parameters.

> "the predictions in Fig. 1 are qualitative rather than quantitative and seem to be based on simple drawings of hypothetical scenarios rather than simulations of actual models. This does not hold up to the standards of the field of decision neuroscience. The authors should implement and simulate their extended model and show that the main qualitative predictions hold across a range of parameters."

**Response.** Figure 1 no longer schematizes the hypotheses: it shows simulated dynamics of the two models, at thirty admissible parameter vectors rather than one, so that what it displays is the spread across the parameter space. On the second half of the comment — that the qualitative predictions should be shown to hold across a range of parameters — the manuscript now does both halves of what the reviewer asks, in the order that lets each carry its own weight.

Each prediction is first derived from the model equations, in the Results where it is used, with the heavier proofs in the Appendix. A derivation establishes an ordering at every parameter setting for which its premises hold, which is more than a demonstration that the ordering held at a sample of settings, and it makes visible which structural feature the prediction depends on.

The predictions are then evaluated over the parameter space, which supplies the size a derivation cannot. We draw 1200 parameter vectors by Latin-hypercube sampling over the ranges of Table 1, simulate 4000 trials per cell per model per draw, and retain the draws that meet behavioral admissibility criteria (160 for Model 1, 145 for Model 2). For every effect we report the range of magnitudes each model produces over its admissible draws, and the fraction of draws at which it produces the observed sign, against the observed value and its between-participant interval (Figure 6).

On four effects the evaluation is the substance of the answer rather than a confirmation, and it is worth saying which way each of them cuts.

Two of them are in the stimulus-onset window, and both exclude Model 1 there — before the delay period is reached at all. The duration effect on the level of the partial-activation rate is +0.057 ± 0.015 s⁻¹ in the data; Model 1 produces the opposite sign at 96% of its admissible draws and the observed sign at none. The duration-dependence of the ramp within that same window repeats the split exactly: observed at +0.060 s⁻², Model 1 gets it backwards at 94% of its draws and right at none, against 37% for Model 2. Both follow from the same structural fact — Model 1's only source of duration dependence in that window is the extra evidence and the extra variance the long stimulus has accumulated, which push the rate up after the long stimulus, whereas the data show it higher after the short one.

The third is the delay-period trend, which is where the models are furthest apart: Model 1's is exactly zero at every parameter setting (Corollary A.3), 78 of Model 2's 145 admissible draws produce all six condition trends significantly positive, and 0 of Model 1's 160 do. The curvature says the same thing a second time, on a statistic Model 1 also fixes at zero (see C1.8).

The fourth runs the other way, and we report it as such: on the coherence × duration interaction on delay-period rate both models are wrong. It is +0.107 ± 0.051 s⁻¹ in the data, and runs in the observed direction at 3% of Model 1's admissible draws and 5% of Model 2's, against the opposite direction at 88% and 86%. We report that as a disagreement between the model and the data rather than as a parameter choice some other setting would repair, and the Discussion says so. It also resolves a tension the original submission carried without remarking on it: coherence-related differences in the maintained level were said to shrink after a long stimulus in the Results, for PMET and RT, and to grow after a long stimulus in the Discussion, for delay-period rate. The saturation of the maintained level predicts the first. The second no longer appears: the manuscript now treats that interaction as a limitation of the model, and does not use it to argue for either version.

We have also tried to be exact about the scope of each derivation. Some results are exact (accuracy is exactly one half at zero coherence whatever the duration; the final response-locked bin is empty by definition; Model 1's delay-period rate is exactly constant). The saturation result is conditional, and we state the condition. And the delay-period derivation itself is stated with the three idealizations that separate the computed rate from the measured one — merging of crossings closer together than 25 ms, absorption at the response criterion, and exclusion of the response-producing burst from the count — all three of which act in the same direction, so that what is derived bounds from above what either model can be measured to do on the trend (Appendix, Remark A.7). After the response signal the criterion becomes an absorbing barrier and the Gaussian structure the delay-period derivation relies on is destroyed, so that window is constrained through the two components of response time rather than through a closed-form rate; we say so rather than implying more.

*Manuscript location*: Figure 1; Results, the "Theoretical predictions" and "Simulation validation" subsections of each window; Appendix; Materials and Methods, "Evaluation over the parameter space"; Figures 5–7.

## C2.2 — The specific ordering in Fig. 1.

> "the slopes of the accumulator dynamics for the two motion coherence conditions are unknown, only their ordering can be safely assumed. Thus, it is not clear whether a long low-coherence stimulus yields a larger cumulative signal than a short high-coherence stimulus, as indicated in Fig. 1."

**Response.** Figure 1 has been replaced: it now shows the two models' simulated trajectories over thirty admissible parameter vectors rather than a schematic ordering. But replacing a drawing with simulations only relocates the assumption into the sampled vectors, so the quantitative answer is better given directly. The reviewer's premise is right — the ordering between a long low-coherence stimulus and a short high-coherence one cannot be assumed — and the account does not need it. The maintained level is capped at the internal threshold θ (Equation 4): a trial whose first accumulation reaches θ holds θ for the rest of the delay whatever its coherence or duration, and a longer stimulus raises the proportion of such trials at every coherence. The two conditions the reviewer contrasts are therefore both near the cap, and their ordering is close to a tie. What the paper rests on is that saturation, which is itself testable: it predicts a coherence effect on response time larger after a short stimulus than after a long one, which the data show, and a coherence effect on delay-period partial-activation rate that shrinks after a long stimulus, which they do not (see C2.1). Nothing in the paper's conclusions depends on the ordering the reviewer questioned.

*Manuscript location*: Figure 1 and its legend; Results, Behavior > Theoretical predictions and Stimulus-offset-locked > Theoretical predictions (the saturation argument, derived from Equations 4 and 5).

## C2.3 — Leak.

> "This is even more true once a non-zero leak is introduced in the integrator, as in generalized diffusion models."

**Response.** The reviewer is right, and we have answered the point in two ways, because they answer different halves of it.

Structurally, adding a leak to the decision variable changes none of the derivations, since none of them uses leak-freeness. Three are worth spelling out. Discriminability at the end of a stimulus of length t is, under a leaky accumulator with leak γ, equal to (kc/σ)√(2/γ)·√(u/(2−u)) with u = 1 − e^(−γt); this increases in coherence, increases in duration through u, is exactly one half at c = 0 by symmetry, and reduces to kc√t/σ as γ → 0. So the accuracy predictions are leak-proof at any leak rate. The response-time predictions require only that the drift of the second accumulation be the same in every condition (Equation 9), so that coherence and duration act on response time through the maintained level alone, and that the maintained level be capped at the internal threshold (Equation 4); a leak alters the mapping from coherence to the maintained level but leaves both properties intact.

The delay-period results are untouched for a third and stronger reason, which the manuscript now states as an appendix result of its own (Remark A.6). A leak acts only before maintenance begins. It changes the law of the accumulated evidence, hence the law of the threshold-crossing time and of the maintained level, and therefore the distribution across trials of the standardized distance from that level to the gate. It does not appear in the maintenance dynamics: from maintenance onset the input to the filter is the maintained level plus the Ornstein–Uhlenbeck fluctuation, and neither depends on how the accumulation arrived. The fluctuation, its moments and the crossing rate are therefore exactly as derived, at any leak rate. Corollary A.3 is thus invariant: under Model 1 the delay-period rate is constant in time on every trial and at every parameter setting, for every γ. A leak can move the level of that rate; it cannot give it a time dependence, and so it is not a route by which Model 1 could acquire a delay-period ramp.

Numerically, we re-ran the parameter-space evaluation at leak rates of 0.25, 0.50 and 1.00 s⁻¹ on the decision variable, over exactly the draws admissible at zero leak and with the same parameter vectors and the same 4000 trials per cell. Nothing directional moves. Model 2 produces all six delay-period condition trends significantly positive at 76, 71 and 77 draws, against 78 at zero leak; Model 1 does so at none of its draws, at any leak rate, as at zero leak. The coherence and the duration effect on response time keep the observed sign at 136 to 138 of Model 2's 145 evaluated draws at every leak rate, against 136 at zero leak. What does erode, mildly, is admissibility itself: 160 of Model 1's draws and 145 of Model 2's are admissible at zero leak, and 138 and 130 at the largest leak. The draws that drop out are those with the smallest ratio of drift to threshold, kc/θ — median 5.4 among the fifteen Model 1 draws lost at γ = 1 s⁻¹, against 7.8 among those kept — which is what the mechanism predicts: a leak caps the accumulator at kc/γ, and once that asymptote approaches θ the threshold stops being reached. That ratio is the criterion for when a leak can matter at all, and over the admissible draws its smallest value at 50% coherence is 1.5 s⁻¹, above the largest rate we tested.

It is worth adding, as the reviewer implicitly does, that the models already contain a first-order filter between accumulator and muscle, which is a leak on the transmitted signal, and an Ornstein–Uhlenbeck process during the delay, which is mean-reversion of the maintenance fluctuation. Neither is what the comment asks for. The leak on the decision variable itself is the one addressed above.

*Manuscript location*: Appendix, Remark A.6 (the structural argument, with the simulation result stated in one sentence); Results, Stimulus-offset-locked > Simulation validation (the same, in the Results); Materials and Methods, Theory (Equations 1–8) and Results, the "Theoretical predictions" subsections, none of which uses leak-freeness. The per-leak-rate counts are reported here rather than in the manuscript; the code and the output are in the deposited archive (`scripts/leak_sweep.py`, `scripts/results/leak/`).

## C2.4 — What is accumulated in the second phase.

> "the response signal is unambiguous and presented at fixation, leaving no sensory evidence to be accumulated before it can be detected even in the presence of strong internal noise. This aspect should be unpacked."

**Response.** We considered renaming this stage but have kept "second accumulation phase." The reviewer is right that the response signal itself carries no directional evidence to resolve, but we intend the second phase as accumulation of the perceptual evidence that the signal is present and has been detected — a well-established construct in the one-choice/detection literature (e.g., Ratcliff & Van Dongen, 2011; Smith, 1995), supported also by electroencephalographic signatures of accumulation under near-ceiling target detection (O'Connell et al., 2012), and distinct from the motion-discrimination accumulation of the first phase. The model uses exactly two properties of this stage: the drive is condition-independent (Equation 9), which places all condition effects on PMET into the starting point, and its increments are independent, which is what makes MET invariant across conditions. Nothing in the model depends on the signal carrying directional evidence, so the reviewer's objection bears on the name we gave this stage, not on the mechanism itself, which would be unchanged under any label.

*Manuscript location*: Introduction; Materials and Methods, model description (response-signal accumulator).

> [NOTE INTERNE — ne pas envoyer telle quelle] Décision explicite de l'utilisateur de NE PAS renommer — la réponse ci-dessus défend la formulation actuelle. À relire pour vérifier que l'argumentaire vous convient avant envoi.

## C2.5 — Scholarship: motor and premotor circuits.

> "studies from the lab of Karel Svoboda have consistently identified a 'ramping mode' of non-selective drive of mouse (pre)motor cortex [...] both the Svoboda work in mice as well as work in monkey indicate that effector-selective motor preparatory activity and activity driving motor execution 'live' in different subspaces."

**Response.** The Introduction already related the model's architecture to this literature: a selective maintained level that carries the choice and does not by itself trigger movement, and a non-selective drive released at the response signal — corresponding, respectively, to the non-selective ramping mode reported in mouse motor circuits and to the output-null/output-potent subspace separation described in primate motor cortex (Inagaki et al., 2022; Churchland & Shenoy, 2024). We have now added a sentence in the Discussion making this correspondence explicit at the point where the EMG evidence for it is presented: the effector-selective partial activations sustained through the delay, followed by a non-selective, response-signal-triggered drive to response, are the empirical counterpart of that architecture.

*Manuscript location*: Introduction (model-architecture paragraph); Discussion (second paragraph).

## C2.6 — Writing: abbreviations and placement of the reasoning.

> "the text contains a large number of abbreviations, many of which are non-standard [...] much of the key reasoning underlying the hypotheses is folded into the EMG part of the Methods section."

**Response.** We have defined every heavily used abbreviation (RS, PMA, LRP, MBL) at its true first occurrence in the running text; several had previously been used well before being spelled out, including PMA, which was not introduced as an abbreviation anywhere in the body text despite being used from its first Results subheading onward. We have not removed these abbreviations, since RS, PMA, PMET, MET, MBL, and LRP each occur between 13 and 94 times and are standard shorthand in this literature; spelling them out throughout would add length without adding clarity. On placement of the reasoning: it has been moved out of Methods altogether. Results is organized by analysis window, and each window opens with a "Theoretical predictions" subsection that derives what the two models predict there from the model equations, immediately before the data that test it.

*Manuscript location*: Introduction ("a response-signal (RS) variant"); Materials and Methods, "EMG-derived measures" ("Partial muscle activations (PMAs)") and "Measurement of motor preparation EEG signals" ("(LRP)" and "(MBL)"); Results, the "Theoretical predictions" subsections.

## C2.7 — Broader evidence on selective motor preparation.

> "p. 21, l. 418 [...] I suggest to extend it by also incorporating the accumulation of sensory evidence in other formats and to animal work."

**Response.** We have added a citation to non-human-primate reaching decisions alongside the existing citations for evidence accumulation at the motor-preparatory level (Cisek & Kalaska, 2005).

*Manuscript location*: Discussion, second paragraph.

> [NOTE INTERNE — ne pas envoyer telle quelle] Traitement partiel — seul le volet « animal work » est couvert. Le volet « autres formats sensoriels » (accumulation de preuves auditives, tactiles, etc.) n'est PAS traité. Voir « Note pour les auteurs ».

# Editorial Note

> "in framing the dynamic regulation of action during working memory and decision making, the authors may also wish to consider other work by van Ede."

**Response.** Added a citation to van Ede & Nobre (2023) in the Introduction, in the context of accounts of working memory as action-oriented.

*Manuscript location*: Introduction, model-architecture paragraph.

# Additional changes made on our own initiative

In preparing this revision we identified three further issues in the original submission, none of which either reviewer raised, but which we think are better disclosed here than left for a reader to discover independently — the data and analysis code are public.

## Volunteered 1 — an unremarked significant effect.

**Response.** The main effect of duration on stimulus-period partial-activation rate was significant with a large effect size and appeared in the summary tables, but no sentence of the original Results or Discussion referred to it. We have added a sentence reporting it; the theory explains it with no additional assumption, as the average counterpart of the duration × time interaction already discussed. It has since become one of the two effects that exclude Model 1 in that window (see C2.1).

*Manuscript location*: Results, Stimulus-onset-locked > Empirical results > Partial muscle activation (PMA) rate.

## Volunteered 2 — a two-phase mu–beta lateralization buildup.

**Response.** As originally written, the claim that mu–beta lateralization shows a coherence-dependent buildup rate across the whole stimulus period was not fully accurate: the coherence-dependent trend is present only in approximately the first second of stimulation and is extinguished thereafter, as coherence-dependent trajectories converge on a common sustained level. We have restated the finding in two phases and added a split-half analysis confirming the pattern directly; this is evidence for accumulation toward an internal threshold, not against it.

*Manuscript location*: Results, Stimulus-onset-locked > Empirical results > Mu–beta lateralization; Tables S4 and S5.

## Volunteered 3 — a post-submission correction to two exported data files.

**Response.** After the original submission we found an error in the two exported electromyographic data files used for the response-signal and response periods, and all results reported for those periods reflect corrected exports. The error was in the export step, not in the recorded data or in the activation-detection procedure. Both files had been written with a 0.3-s window rather than the 0.2-s window specified in Materials and Methods, so each contained three 0.1-s bins per cell instead of two; and in the response-locked file the bins were referenced to the onset of the response-producing muscle burst rather than to the button press. The two bins actually analyzed in the original submission were consequently the second and third of those three, that is the interval from roughly 0.4 to 0.2 s before muscle-activation onset, rather than the 0.2 s immediately preceding the button press. The response-signal period was affected only by the window length and its results are essentially unchanged; the response period was affected by both, and its numbers change.

The consequence for the response period is that the main effect of duration on partial-activation rate, null as originally analyzed (F(1,21) = 0.011, p = .919), is significant in the corrected export (F(1,21) = 7.61, p = .012, ηp² = 0.266), in the same direction as the response-signal period. A second consequence is structural rather than numerical, and explains why the manuscript now analyzes this period on one bin rather than two. Once the bins are correctly referenced to the button press, the final bin (−0.1 to 0 s) is exactly zero for every participant and every condition, and necessarily so: a partial muscle activation is by definition an activation preceding the response-producing burst, and that burst begins on average about 160 ms before the press. The manuscript states this as the structural fact it is, without recounting the correction; we describe the correction here so that a reader who re-runs the analysis from the originally deposited export, and obtains F = 0.011, can see why. No other period or dependent measure was affected, and the corrected exports replace the originals in the archive deposited with the paper.

*Manuscript location*: Results, Response-locked > Empirical results > Partial muscle activation (PMA) rate; Materials and Methods, Statistical Analysis (final paragraph, the structurally empty final bin). The correction itself is documented here rather than in the manuscript.

> [NOTE INTERNE — ne pas envoyer telle quelle] Vérifier avant envoi que l'archive OSF a bien été mise à jour avec les exports corrigés du 19/08 (`PMA_rate_resp_JASP.csv` et `PMA_rate_rs_JASP.csv`) : la phrase ci-dessus l'affirme au rapporteur. Sur le Mac, les fichiers du 05/08 sont conservés dans `data/86cbw-osfstorage-archive/analysis_PMA/data/_superseded_aug5/`.
