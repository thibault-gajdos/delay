Achard, Gajdos Preuss & Servant — "Graded Motor-Preparatory Dynamics Preserve Perceptual Decisions over Short Delays" (Journal of Neuroscience, R1)

DRAFT — internal working copy, not yet formatted for submission.

## Note pour les auteurs — à lire avant envoi

Cette lettre est mise à jour contre `manuscript/latex/delay.tex` du 25 octobre 2026, qui intègre la v5 (Introduction, Methods) et la nouvelle évaluation de l'espace des paramètres. Les notes internes sont repérées dans le texte (« NOTE INTERNE »). Points à régler avant l'envoi :

-   **C2.4 (« second accumulation phase »)** : réponse inchangée sur le fond (on garde le nom), à relire.
-   **C2.7 (autres formats sensoriels)** : traitement partiel, comme avant.
-   **C1.7 (α fonction du délai)** : comme avant.
-   **Plages de r et g** : le texte dit que g est échantillonné au-dessus de l'estimation de Dendauw et al. (2024). Les valeurs de r et g de Dendauw et al. n'ont pas été vérifiées.
-   **Archive OSF** : vérifier que les exports corrigés du 19/08 (`PMA_rate_resp_JASP.csv`, `PMA_rate_rs_JASP.csv`) y sont en place (Volunteered 3).

### Ce qui a changé depuis la version du 8 septembre

-   **Texte de la v5** : Introduction et Methods réécrits ; PMET et MET deviennent premotor time (PMT) et motor time (MT) ; la section « EMG-derived measures » et le paragraphe sur la segmentation des bouffées disparaissent ; le cadre théorique est dans les Methods.
-   **Statistiques** : dans les fenêtres stimulus et délai, chaque profil est résumé par son niveau moyen, sa pente et sa courbure (coefficients des moindres carrés, en s⁻¹, s⁻² et s⁻³ pour les PMA), chacun analysé par une ANOVA durée × cohérence. Dans les fenêtres RS et réponse, les PMA et le MBL sont décrits sans test ; le LRP est testé sur son niveau et sa pente.
-   **Simulations** : aucun critère de détection EMG appliqué aux traces simulées ; 40 000 essais par condition ; critères d'admissibilité sans borne tirée des données (228 vecteurs admissibles pour le Modèle 1, 204 pour le Modèle 2) ; deux fuites indépendantes, grille 4 × 4.
-   **Règle de comparaison (v5)** : un vecteur reproduit un effet si son intervalle Monte Carlo est contenu dans l'intervalle de confiance empirique à 95 %. On rapporte minimum, maximum et effectifs ; les signes seulement pour les effets dont la direction n'est pas fixée par une dérivation.
-   **Results réorganisés** : « Accuracy, RT, PMT and MT », puis une sous-section par fenêtre ; Table 2 = comparaison modèles–données.
-   **Annexe** : Proposition A.5 remplacée (décomposition exacte ; signe établi numériquement) ; Remarque A.6 couvre les deux fuites.

# Reviewer 1

## C1.1 — Slopes rather than bins.

> "Although there is a nice cartoon figure for the hypotheses, the authors use a repeated measures ANOVA for their tests. I would prefer tests on the different slopes in the EMG and EEG instead of the arbitrary 0.1s bins. [...] Ln 324: 'slope of this ramping activity' → that is not what you ANOVA tested."

**Response.** We agree, and the analysis now tests slopes directly. For the stimulus-onset and delay windows, each participant's profile in each coherence × duration condition is summarized by three numbers: its mean level across bins, its linear slope, and its quadratic curvature. The slope and the curvature are least-squares coefficients on bin-center times expressed in seconds and centered within the window; for the partial-activation rate they are in s⁻² and s⁻³. Each of the three numbers is then analyzed by its own coherence × duration repeated-measures ANOVA, and the mean slope and mean curvature are tested against zero. No signal was reprocessed: these are weighted sums of the bin means already computed. The curvature is reported because it is where the two candidate models make opposite predictions about the delay window (see C1.8), and because mu–beta lateralization during the stimulus is non-monotone, so that a slope alone would misdescribe it. We add that the omnibus time-bin test the comment objects to was also poorly behaved in its own right here, with a Greenhouse–Geisser epsilon of 0.14 on 19 degrees of freedom.

For the partial-activation rate, the slope and curvature are taken on the rate itself, not on its logarithm. The logarithm is undefined wherever a cell contains no activation: 6% of the participant × condition × bin cells of the delay window, 30% of those of the stimulus-onset window. Deleting or offsetting those cells would change the estimand by an amount that depends on the rate, which is the quantity under test.

The response-signal and response windows contain two bins each, and are treated differently. The partial-activation rate and mu–beta lateralization are described there rather than tested: no partial activation can occur after the onset of the response-related burst, so a decline is expected by construction, and successive spectral estimates overlap too much for mu–beta lateralization to be resolved in 0.2 s. The lateralized readiness potential, whose temporal resolution is higher, is tested on its level and its slope.

The change of analysis has one consequence worth stating explicitly. For mu–beta lateralization in the stimulus-onset window, the coherence dependence is significant on the omnibus time-bin test but not on the full-window slope (p = .596), because the measure rises over roughly the first second and then flattens. The Results say so, and a split of the window at its midpoint (Table S5) locates the coherence-dependent slope in the first half. Conversely, the duration effect on the full-window slope of mu–beta lateralization is significant (p = .046) where the omnibus test was not (p = .136); it is reported in Table S4.

Across the delay, the partial-activation rate has a positive slope (0.560 s⁻² [0.364, 0.757], F(1,21) = 35.24, p < .001) whose coherence dependence, the marker of continued evidence accumulation, is null (p = .567). The Abstract accordingly says that the rate "continued to ramp during the stimulus-free delay, but no longer at a coherence-dependent rate."

*Manuscript location*: Materials and Methods, "Statistical analysis"; Results, the stimulus-onset, delay, response-signal and response windows; Tables S2–S5; Abstract.

## C1.2 — Clearer statement of the hypotheses and of what each measure does.

> "I would like a much clearer description of the hypotheses and how the different measures help unraveling. Now the EEG results are a bit of a hotch-potch."

**Response.** The Introduction now states the hypothesis and the two candidate models in its last two paragraphs, and says which measure tests what: partial muscle activations index motor preparation, Model 1 predicts a constant rate of them during the delay and Model 2 a rising one, and the two EEG measures are complementary cortical indices of effector-selective preparation. The Results are organized accordingly. A first subsection takes accuracy, response time, premotor time and motor time. Each following subsection takes one analysis window: stimulus-onset, delay, response-signal, and response. Each states what the two models predict there before any data are shown ("Theoretical predictions"), then reports the data, measure by measure (partial-activation rate, mu–beta lateralization, lateralized readiness potential), and, for accuracy, response time and the partial-activation rate in the stimulus-onset and delay windows, compares the observed effect magnitudes with those the models produce ("Model–data comparison"). Materials and Methods also state why the EEG measures are used qualitatively: the model specifies the neural drive to the muscles but no forward model from its internal variables to scalp potentials or spectral amplitude.

*Manuscript location*: Introduction, last two paragraphs; Materials and Methods, "Theoretical framework" (last paragraph) and "EEG acquisition and preprocessing"; Results, all subsections.

## C1.3 — Muscle-activation onset procedure.

> "Ln 119: please be more explicit about the EMG onset procedure. [...] p. 10, ll. 157-159: PMA should be more precisely defined. Specifically, how is it decided that a muscle activation detected by EMG does not 'culminate in an overt response'?"

**Response.** Materials and Methods now give the full procedure: candidate bursts were detected at seven baseline standard deviations on the rectified signal; threshold-defined periods separated by less than 25 ms were merged; an adaptive minimum-sample criterion rejected noise-driven detections; burst boundaries were refined with the integrated-profile method (Liu & Liu, 2016; Santello & McDonagh, 1998); and every trial was inspected visually. The classification is now stated explicitly. In each retained trial, the response-related muscle activation is the last burst in the responding hand before the button press; its onset splits response time into premotor time and motor time. Every other burst beginning between stimulus onset and the button press is a partial muscle activation. In the simulated traces, a muscle activation begins each time the neural drive becomes positive; the last one before the response is the response-related activation, every other one is a partial activation, and no further detection criterion is applied.

*Manuscript location*: Materials and Methods, "EMG acquisition and preprocessing" (second paragraph) and "Evaluation over the parameter space" (first paragraph).

## C1.4 — Coherence range used by others.

> "Ln 242: make explicit what range of coherences others have used."

**Response.** Added as a scope condition in the Discussion: our three coherence levels (0%, 12%, 50%) are coarser than the five- or six-level logarithmic series commonly used to characterize psychometric and chronometric functions in random-dot-motion studies (e.g., 3.2–51.2%: Roitman & Shadlen, 2002; 6.4–51.2%: Kiani et al., 2008). They are sufficient to establish the coherence-dependent effects reported here, though not their full functional form.

*Manuscript location*: Discussion, scope-conditions paragraph.

## C1.5 — Error bars and individual data.

> "Fig 3: SEM is not a good measure on data constrained between 0 and 1 [...] In modern research it is custom to also show individual participant datapoints."

**Response.** Figure 3 has been regenerated. Each participant's data now appear as faint connected points behind the group means in all four panels. For accuracy, which is bounded in [0, 1], the interval is now a 95% percentile bootstrap interval computed on the Cousineau–Morey-normalized values (10,000 resamples over participants, Morey-scaled). It is a within-subjects interval, comparable to the ±1 within-subjects SEM shown for response time, premotor time and motor time, and it is asymmetric where the data approach the bounds. The legend describes both intervals.

*Manuscript location*: Figure 3 and its legend.

## C1.6 — Reminder of MET and PMET.

> "Ln 308: please refer back to fig 2 to reiterate what MET and PMET are."

**Response.** The two components are now called premotor time (PMT) and motor time (MT). They are defined in Materials and Methods, illustrated in Figure 2A, and the Results refer back to Figure 2A where mean PMT and MT are first reported.

*Manuscript location*: Materials and Methods, "Theoretical framework" (last paragraph but one); Figure 2A; Results, "Accuracy, RT, PMT and MT" > Empirical results.

## C1.7 — Different delays, to show build-up.

> "Ln 485 [...] If you really want to hammer in this point, you should show build-up for different delays (even within the same trial, where the RS can come at for example 1.5 or 3.0 sec)."

**Response.** Added as a second scope condition in the Discussion, citing a companion study (Achard et al., 2025) that manipulates delay length directly: the same delay-period build-up is preserved across 3–7 s delays when the response mapping is known in advance, but is largely absent when it is revealed only at the response signal.

*Manuscript location*: Discussion, scope-conditions paragraph.

> [NOTE INTERNE — ne pas envoyer telle quelle] Le point plus fin (« α fonction du délai, condition de portée ») n'est pas formulé dans ces termes.

## C1.8 — A quantitative model.

> "Ln 505-506: here you hint at a quantitative model, which I believe would strengthen your data better. Now you use ANOVAs, and thus significances to rule out hypotheses that are only expressed in words or cartoons, but are never quantified based on the data."

**Response.** We agree, and the paper now has a quantitative model, used in two stages.

The first stage is derivation. Every prediction the paper tests is derived from the model equations where it is used, in the Results, with the longer proofs in a single Appendix. Three examples. Because the drift of the second accumulation is the same in every condition (Equation 9), coherence and duration can act on response time only through the motor-preparatory level reached at the end of the first stage. Because that level is capped at the internal threshold (Equation 4), the coherence effect on response time is larger after a short stimulus than after a long one. And because temporal preparation is anchored to stimulus offset while the response signal follows it at a fixed interval, stimulus duration cancels from the temporal-preparation term: it contributes the same amount in all six conditions, so no effect on accuracy or response time can separate the two models.

The second stage is an evaluation over the models' parameter space, which gives the size of each effect, which a derivation does not. We do not fit participant-level parameters: doing so would first require showing that the parameters can be recovered, through a dedicated identifiability analysis. We sample instead. We draw 1,200 parameter vectors by Latin-hypercube sampling over the ranges of Table 1, simulate 40,000 trials per condition per model per vector, and retain the vectors that meet admissibility criteria stated in Materials and Methods, none of which is derived from the data (228 vectors for Model 1, 204 for Model 2). For each effect, the observed group mean and its 95% confidence interval across participants are compared with the effect each admissible vector produces, with its Monte Carlo 95% interval. A vector reproduces the observed magnitude if its interval lies inside the observed one, does not if the two are disjoint, and is uncertain otherwise. A model can reproduce an effect if at least one admissible vector does. Table 2 and Figure 6 report, for each effect and model, the minimum and the maximum over the admissible vectors and the number of vectors in each category.

Accuracy, response time and premotor time do not separate the models: each effect is reproduced by both models, by similar numbers of vectors. The delay window separates them, and not by degree. Because the maintenance fluctuation starts at its stationary amplitude, Model 1's delay-window rate of partial activations is exactly constant in time, at every value of every parameter (Appendix, Corollary A.3); its slope and its curvature are therefore zero. In the simulations, Model 1's slope ranges from −0.090 to 0.007 s⁻² over its admissible vectors, and no vector reproduces the observed 0.560 s⁻² [0.364, 0.757]; the negative values come from the exclusion of the response-related burst from the count (Appendix, Remark A.7). Model 2's slope ranges from −0.029 to 1.064 s⁻², and 7 of its vectors reproduce the observation. Figure 5 shows the profiles: every admissible vector of Model 1 is flat, whereas Model 2's fan out.

The curvature gives a second test in the same window. Under Model 1 it is zero. Under Model 2 the mean rate of a condition is convex wherever the rate-weighted mean of the squared standardized distance to the gate exceeds 1 (Appendix, Corollary A.4). The observed curvature is 0.210 s⁻³ [0.109, 0.310] (F(1,21) = 18.89, p < .001), and neither duration nor coherence modulated it. No Model 1 vector reproduces it; 12 Model 2 vectors do.

*Manuscript location*: Results, the "Theoretical predictions" subsections (derivations); Appendix, Proposition A.1 to Remark A.7 (proofs); Materials and Methods, "Theoretical framework" (last paragraph), "Evaluation over the parameter space" and "Model–data comparison"; Results, the "Model–data comparison" subsections; Table 2; Figures 5–7.

## C1.9 — Statistics in the text.

> "In J Neurosci it is customary to present statistics in the results text instead of everything in a table."

**Response.** Test statistics (F, degrees of freedom, p, ηp²) are now reported in the Results text throughout, in addition to the supplementary tables.

*Manuscript location*: Results, throughout.

# Reviewer 2

## C2.1 — Simulate the model; show the predictions hold across parameters.

> "the predictions in Fig. 1 are qualitative rather than quantitative and seem to be based on simple drawings of hypothetical scenarios rather than simulations of actual models. This does not hold up to the standards of the field of decision neuroscience. The authors should implement and simulate their extended model and show that the main qualitative predictions hold across a range of parameters."

**Response.** Figure 1 no longer schematizes the hypotheses: it shows the simulated dynamics of the two models at thirty admissible parameter vectors, so that it displays the spread across the parameter space. On the second half of the comment, the manuscript now does two things.

Each prediction is first derived from the model equations, in the Results where it is used, with the longer proofs in the Appendix. A derivation establishes an ordering at every parameter setting for which its premises hold, and shows which structural feature the prediction depends on.

The predictions are then evaluated over the parameter space: 1,200 parameter vectors, 40,000 trials per condition, both models, with the comparison rule described in C1.8 (Table 2, Figure 6). Four results carry the answer.

The first two are in the stimulus-onset window, and both exclude Model 1 there. The partial-activation rate is higher after the short stimulus than after the long one (0.057 s⁻¹ [0.026, 0.088]). In Model 1, a trial that is still accumulating cannot cross the gate, because accumulation stops at the internal threshold, which lies below the gate (Appendix, Proposition A.5). The duration effect in this window is then the difference between two terms: the rate of the trials that reach the threshold, the level closest to the gate, during the further accumulation of the long stimulus, and the rate of the trials that the short stimulus left maintained below the threshold. The derivation does not fix which dominates. Numerically, the first does: Model 1 produces a higher rate after the long stimulus at 227 of its 228 admissible vectors, an uncertain sign at the remaining one, and never reproduces the observed effect. Model 2's temporal preparation begins at stimulus offset, inside this window only after the short stimulus, and opposes that effect; 10 of its vectors reproduce the observation, and the sign is the observed one at 76. The duration dependence of the slope in that window repeats the split (0.060 s⁻² [0.020, 0.100]): Model 1 produces the opposite sign at 222 vectors and never reproduces it; 21 Model 2 vectors reproduce it.

The third is the delay-window slope and curvature, where the models are furthest apart: Model 1's are exactly zero at every parameter setting (Corollary A.3), and only Model 2 reproduces the observed values (see C1.8). The duration dependence of the slope separates them the same way: no Model 1 vector reproduces it, 37 Model 2 vectors do.

The fourth runs the other way, and we report it as such. The coherence × duration interaction on the delay-window rate (0.107 s⁻¹ [0.001, 0.213]) has the opposite sign at most admissible vectors of both models (206 of 228 in Model 1, 188 of 204 in Model 2); only 10 vectors in each model reproduce it. We report this as a disagreement between the data and both models, and the Discussion says so. It also resolves a tension the original submission carried without remarking on it: coherence-related differences in the maintained level were said to shrink after a long stimulus in the Results, for premotor time and response time, and to grow after a long stimulus in the Discussion, for the delay-window rate. The saturation of the maintained level predicts the first, and the manuscript no longer uses the second to argue for either model.

The evaluation also exposes a smaller discrepancy, which we report. The duration effect on motor time is null in the data (−0.5 ms [−2.1, 1.1]). Both models reproduce it at some vectors (63 and 51), but at others they make motor time shorter after the short stimulus, by up to 82 and 90 ms. At these vectors, the response-related activation begins before the response signal on a larger proportion of trials after the long stimulus; on those trials, motor time includes the interval up to the response signal. In the data, the response-related activation almost never begins before the response signal.

We have also tried to be exact about the scope of each derivation. Some results are exact (accuracy is one half at zero coherence whatever the duration; Model 1's delay-window rate is constant). Others are conditional, and the condition is stated: the coherence × duration interaction on the delay-window rate is predicted to be negative only where the 0%-coherence, 2-s condition has mostly reached the threshold. The delay-window derivation is stated with the three idealizations that separate the computed rate from the measured one (Appendix, Remark A.7): merging of activity periods closer than 25 ms in the data, absorption at the response criterion, and exclusion of the response-related burst from the count. The first cannot create a slope under Model 1 and the last two lower the measured slope, so the flatness of Model 1 is not an artefact of them. After the response signal the criterion acts as an absorbing barrier, the closed-form rate no longer applies, and the partial-activation rate in that window is described rather than compared with the models.

*Manuscript location*: Figure 1; Results, the "Theoretical predictions" and "Model–data comparison" subsections; Appendix; Materials and Methods, "Evaluation over the parameter space" and "Model–data comparison"; Table 2; Figures 5–7.

## C2.2 — The specific ordering in Fig. 1.

> "the slopes of the accumulator dynamics for the two motion coherence conditions are unknown, only their ordering can be safely assumed. Thus, it is not clear whether a long low-coherence stimulus yields a larger cumulative signal than a short high-coherence stimulus, as indicated in Fig. 1."

**Response.** Figure 1 has been replaced by simulated trajectories at thirty admissible parameter vectors. The reviewer's premise is right: the ordering between a long low-coherence stimulus and a short high-coherence one cannot be assumed, and the account does not need it. The maintained level is capped at the internal threshold θ (Equation 4): a trial whose first accumulation reaches θ holds θ for the rest of the delay whatever its coherence or duration, and a longer stimulus raises the proportion of such trials at every coherence. What the paper uses is this saturation, which is testable: it predicts a coherence effect on response time and premotor time that is larger after a short stimulus than after a long one, which the data show. Nothing in the paper's conclusions depends on the ordering the reviewer questioned.

*Manuscript location*: Figure 1 and its legend; Results, "Accuracy, RT, PMT and MT" > Theoretical predictions; Delay window > Theoretical predictions.

## C2.3 — Leak.

> "This is even more true once a non-zero leak is introduced in the integrator, as in generalized diffusion models."

**Response.** The model has two accumulators, the motion-evidence accumulator and the response-signal accumulator, and we have added a leak to each, independently. We answer structurally and numerically.

Structurally, a leak changes none of the derivations, because none of them uses the absence of leak (Appendix, Remark A.6). For accuracy: with a leak γ on the motion-evidence accumulator, discriminability at the end of a stimulus of length t is (kc/σ)√(2/γ)·√(u/(2−u)) with u = 1 − e^(−γt); it increases with coherence and with duration, is one half at c = 0 by symmetry, and reduces to kc√t/σ as γ tends to 0. For the delay window: a leak on the motion-evidence accumulator acts only before maintenance begins. It changes the law of the maintained level, and therefore the distribution across trials of the distance from that level to the gate, but it does not enter the maintenance dynamics. The fluctuation, its moments and the crossing rate are therefore exactly as derived, and Model 1's delay-window rate remains constant in time: a leak can move its level, not give it a time dependence. A leak on the response-signal accumulator acts only after the response signal, after the delay window has ended. After the response signal, the conditions still act on response time only through the value of the signal reaching the gate at that time, whatever the leak.

Numerically, we varied the two leak rates independently over 0, 0.25, 0.5 and 1 s⁻¹ (integration time constants of 4, 2 and 1 s for the non-zero values), which gives 15 combinations besides the main evaluation. Each combination was run on the 304 parameter vectors admissible in either model without leak, with the same parameter values and 40,000 trials per condition. The number of admissible vectors ranges from 218 to 228 for Model 1 and from 194 to 204 for Model 2 across the grid, and no combination changes whether a model reproduces any effect of Table 2.

*Manuscript location*: Materials and Methods, "Evaluation over the parameter space" (last paragraph); Table 1 (γ_x, γ_q); Appendix, Remark A.6; Results, Delay window > Model–data comparison ("Leakage"); Discussion.

## C2.4 — What is accumulated in the second phase.

> "the response signal is unambiguous and presented at fixation, leaving no sensory evidence to be accumulated before it can be detected even in the presence of strong internal noise. This aspect should be unpacked."

**Response.** We considered renaming this stage but have kept the term. The reviewer is right that the response signal carries no directional evidence to resolve, but the second phase accumulates the perceptual evidence that the signal is present. This is an established construct in the one-choice and detection literature (e.g., Ratcliff & Van Dongen, 2011; Smith, 1995), supported by electroencephalographic signatures of accumulation under near-ceiling target detection (O'Connell et al., 2012), and it is distinct from the motion-discrimination accumulation of the first phase. The Introduction and the model description now say this explicitly. The model uses one property of this stage: its drift does not depend on the condition (Equation 9), which places every condition effect on premotor time in the starting point. Nothing in the model depends on the signal carrying directional evidence.

*Manuscript location*: Introduction, third paragraph; Materials and Methods, "Theoretical framework" > RS-driven response initiation.

> [NOTE INTERNE — ne pas envoyer telle quelle] Décision de ne pas renommer ; argumentaire à relire.

## C2.5 — Scholarship: motor and premotor circuits.

> "studies from the lab of Karel Svoboda have consistently identified a 'ramping mode' of non-selective drive of mouse (pre)motor cortex [...] both the Svoboda work in mice as well as work in monkey indicate that effector-selective motor preparatory activity and activity driving motor execution 'live' in different subspaces."

**Response.** The Discussion now relates the model's architecture to this literature, at the point where the EMG evidence for it is presented. The model has a selective maintained state that carries the choice without triggering movement, and a non-selective drive released at the response signal. These correspond to the output-null/output-potent separation described in primate motor cortex and to the subcortical control signals that reshape motor-planning dynamics and promote movement execution in mice (Inagaki et al., 2022; Churchland & Shenoy, 2024). The effector-selective partial activations sustained through the delay, followed by a response-signal-triggered drive to the response, are the empirical counterpart of that architecture.

*Manuscript location*: Discussion, second paragraph.

## C2.6 — Writing: abbreviations and placement of the reasoning.

> "the text contains a large number of abbreviations, many of which are non-standard [...] much of the key reasoning underlying the hypotheses is folded into the EMG part of the Methods section."

**Response.** We replaced the two non-standard abbreviations, PMET and MET, by premotor time (PMT) and motor time (MT), the terms used in the EMG literature. The remaining abbreviations (RS, RT, PMA, MBL, LRP, EMG, EEG) are each defined at their first occurrence. The reasoning has been moved out of the EMG part of the Methods. The model is now presented in its own section of Materials and Methods ("Theoretical framework"), and each Results subsection opens with "Theoretical predictions", which derives what the two models predict there, immediately before the data that test it.

*Manuscript location*: Introduction (RS, RT, EEG, EMG, PMA); Materials and Methods, "Theoretical framework" (PMT, MT) and "EEG acquisition and preprocessing" (MBL, LRP); Results, the "Theoretical predictions" subsections.

## C2.7 — Broader evidence on selective motor preparation.

> "p. 21, l. 418 [...] I suggest to extend it by also incorporating the accumulation of sensory evidence in other formats and to animal work."

**Response.** We have added citations to non-human-primate work alongside the existing citations for evidence accumulation at the motor-preparatory level (Cisek & Kalaska, 2005; Gold & Shadlen, 2000), and to human work showing that stretch-reflex gains and motor-cortical activity track evolving evidence (Selen et al., 2012; Donner et al., 2009).

*Manuscript location*: Materials and Methods, "Theoretical framework" (filter paragraph); Discussion, second paragraph.

> [NOTE INTERNE — ne pas envoyer telle quelle] Traitement partiel : le volet « autres formats sensoriels » n'est pas traité.

# Editorial Note

> "in framing the dynamic regulation of action during working memory and decision making, the authors may also wish to consider other work by van Ede."

**Response.** Added a citation to van Ede & Nobre (2023) in the Introduction, in the context of accounts of working memory as action-oriented.

*Manuscript location*: Introduction, third paragraph.

# Additional changes made on our own initiative

In preparing this revision we identified further issues, none of which the reviewers raised, which we think are better disclosed here.

## Volunteered 1 — an unremarked significant effect.

**Response.** The main effect of duration on the stimulus-onset partial-activation rate was significant with a large effect size and appeared in the summary tables, but no sentence of the original Results or Discussion referred to it. It is now reported, and it has become one of the two effects that exclude Model 1 in that window (see C2.1).

*Manuscript location*: Results, Stimulus-onset window > Empirical results > PMA rate.

## Volunteered 2 — a two-phase mu–beta lateralization buildup.

**Response.** As originally written, the claim that mu–beta lateralization shows a coherence-dependent buildup across the whole stimulus period was not accurate: the coherence-dependent slope is present only in approximately the first second and disappears thereafter, as the trajectories converge on a common sustained level. We now state the finding in two phases and report a split-half analysis (Table S5); this is consistent with accumulation toward an internal threshold.

*Manuscript location*: Results, Stimulus-onset window > Empirical results > MBL; Tables S4 and S5.

## Volunteered 3 — a post-submission correction to two exported data files.

**Response.** After the original submission we found an error in the two exported electromyographic data files used for the response-signal and response windows. The error was in the export step, not in the recorded data or in the detection procedure. Both files had been written with a 0.3-s window rather than the 0.2-s window specified in Materials and Methods, so each contained three 0.1-s bins per cell instead of two; and in the response-locked file the bins were referenced to the onset of the response-related burst rather than to the button press. All results for those windows now reflect the corrected exports. With the bins referenced to the button press, the final bin (−0.1 to 0 s) is exactly zero for every participant and condition, as it must be: a partial activation precedes the response-related burst by definition, and that burst begins on average about 160 ms before the press. This is one of the reasons why the partial-activation rate is now described rather than tested in the response-signal and response windows (see C1.1). For a reader who re-runs the original analysis: the main effect of duration in the response window, null with the original export (F(1,21) = 0.011, p = .919), is significant with the corrected one (F(1,21) = 7.61, p = .012). The corrected exports replace the originals in the deposited archive.

*Manuscript location*: Results, Response window > Empirical results > PMA rate. The correction itself is documented here.

> [NOTE INTERNE — ne pas envoyer telle quelle] Vérifier que l'archive OSF contient les exports corrigés du 19/08.

## Volunteered 4 — simulation procedure and admissibility.

**Response.** Three choices in the parameter-space evaluation were revised so that no part of it depends on the data it is compared with. First, no EMG detection criterion is applied to the simulated traces: every gate crossing is an activation (see C1.3). Second, the admissibility criteria contain no bound derived from the data: in every condition, the response-related activation may begin before the response signal on at most 10% of trials, responses may occur within 0.2 s of the signal on at most 10%, and failures to respond by the 4-s deadline on at most 5%; across conditions, the mean delay-window partial-activation rate must lie between 0.02 and 5 s⁻¹. Third, 40,000 trials are simulated per condition, so that no verdict of Table 2 depends on a vector whose Monte Carlo interval overlaps the observed one only partly.

*Manuscript location*: Materials and Methods, "Evaluation over the parameter space".
