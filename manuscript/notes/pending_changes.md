# Modifications à faire dans les manuscrits

Liste tenue à jour au fil de la révision, à appliquer en une fois.
Chaque modification sera soumise pour accord avant d'être faite.

Textes concernés :
- **v5** : `manuscript/Achard et al R1_v5.docx` (Mathieu Servant, 21 octobre 2026) ;
- **delay.tex** : `manuscript/latex/delay.tex` ;
- **lettre** : `manuscript/response_letter_2026-09-08.md` (et .docx).

Question ouverte, à trancher avant toute modification : quel texte fait référence,
la v5 ou delay.tex ?

## 1. Détection des activations dans les simulations (décidé le 22 octobre)

Aucun post-traitement EMG n'est appliqué aux simulations. Une activation commence à chaque
franchissement de la porte vers le haut : pas de fusion, pas de nombre minimal
d'échantillons, pas de raffinement. L'excursion qui atteint r est la bouffée de réponse.

- **v5, note de bas de page 1.** Elle dit que l'application des critères de détection EMG
  aux simulations « yielded similar results ». Cette affirmation n'a pas été vérifiée avec
  le pipeline actuel. Il faut soit la vérifier, soit retirer la note.
- **delay.tex, l. 389.** Décrit la fusion à 25 ms et le minimum adaptatif d'échantillons
  appliqués aux simulations. À remplacer par la règle ci-dessus.
- **delay.tex, Appendix, l. 988 (Remark).** Parle de l'intervalle fusionné et de la
  fusion à 25 ms. Ne s'applique plus aux simulations.
- **delay.tex, l. 207.** Porte sur la détection dans les données ; à relire pour vérifier
  qu'aucune phrase n'y suppose la même règle dans les simulations.
- `scripts/theory_checks.py` compte encore les grappes fusionnées : à aligner (code, pas
  manuscrit).

## 2. Nombre d'essais simulés (décidé le 23 octobre)

1200 jeux de paramètres (même hypercube), 40 000 essais par condition, par modèle et par
jeu.

- **v5** : « 4,000 trials were simulated per experimental condition and model » → 40,000.
- **delay.tex, l. 394** : « 4000 trials per condition per model per draw » → 40000.

## 3. Critère d'admissibilité (décidé le 24 octobre)

Aucune borne n'est calculée à partir des données. Critères, dans chaque condition :
au plus 10 % d'essais (sur tous les essais) dont la bouffée de réponse commence avant le
signal ; au plus 10 % de réponses dans les 0,2 s après le signal ; au plus 5 % d'absences
de réponse avant 4 s. Sur l'ensemble des conditions : taux moyen de PMA pendant le délai
entre 0,02 et 5 s⁻¹.

- **v5** (Evaluation over the parameter space) :
  - remplacer « responses occurring before the RS could occur on no more than 10% of
    trials » par le critère sur les bouffées de réponse commençant avant le signal ;
  - supprimer « Accuracy was additionally required to exceed 80% when motion coherence
    was 50% » ;
  - garder la plage 0,02–5 s⁻¹.
- **delay.tex, l. 394** : même critère ; supprimer la plage 0,04–3 s⁻¹, les fenêtres de
  premotor time (181–677 ms) et de motor time (46–288 ms), et la phrase qui les justifie
  par l'étendue observée.

## 4. Plages des paramètres

- **delay.tex, l. 391** : « They come from the design of the experiment, from published
  estimates, from ranges observed in the data, or from the structure of the model ».
  Aucune plage ne vient des données : supprimer « from ranges observed in the data ».
- **v5** : dit que k, λ et r couvrent « roughly a factor of two on either side » des
  estimations de Dendauw et al. (2024). Dans le code, λ ∈ [10, 50] est justifié par la
  constante de temps du filtre (100 à 20 ms), et r ∈ [0,01 ; 0,2] n'a pas de source
  indiquée. À vérifier contre les valeurs de Dendauw et al., puis corriger le texte ou les
  plages.

## 5. Fuite (décidé le 22 octobre)

Les fuites des deux accumulateurs varient indépendamment : γ_x (accumulateur du
mouvement) et γ_q (accumulateur du signal de réponse), grille complète
{0 ; 0,25 ; 0,5 ; 1} s⁻¹, sur les jeux admissibles à fuite nulle.

- **v5** : « we imposed the same leakage rate γ on both accumulators and examined
  γ = 0.25, 0.50, and 1.00 s⁻¹ » → grille indépendante.
- **delay.tex, l. 578** : parle d'une fuite sur la seule variable de décision, à trois
  valeurs → grille indépendante, nouveaux résultats.
- **lettre, C2.3** : table à refaire avec la grille.

## 6. Règle de comparaison modèle–données

- **v5** (Model–data comparison) : la règle (intervalle Monte Carlo contre intervalle de
  confiance empirique à 95 %) est celle de `scripts/classify.py`. Rien à changer sur la
  règle.
- **delay.tex, l. 398 et Results** : décrit l'ancienne règle (médiane et percentiles
  5–95, fraction des jeux du signe observé). À remplacer par celle de la v5.
- **delay.tex, l. 572** : l'argument sur l'écart-type Monte Carlo mesuré en
  re-simulant quarante jeux (`effect_sizes.py`, `ramp_check.py`) n'a plus lieu d'être
  avec la règle de la v5.
- **Analyse du signe.** Mathieu propose de la supprimer quand une dérivation analytique
  fixe le signe. À décider effet par effet (l'effet de la durée sur le taux pendant le
  stimulus n'a pas de dérivation de signe générale).

## 7. Points de la v5 à vérifier

- **Courbure.** La v5 la définit comme le coefficient de t² (s⁻³). Le script R rapporte un
  contraste orthogonal sans unité (`contr.poly`) ; les F sont identiques. Les valeurs
  citées doivent être les coefficients (Python, `targets.py`).
- **Pente.** Vérifié : régression sur les centres des bins en secondes, en Python et en R.
- **Fenêtres du signal de réponse et de la réponse.** La v5 traite les PMA de ces fenêtres
  de façon descriptive ; `run_stats.R` y fait des ANOVA, et delay.tex rapporte
  F(1,21) = 7,613 pour l'effet de la durée dans la fenêtre de réponse. À décider.

## 8. Chiffres et figures à mettre à jour après la relance

Tous les nombres issus des simulations, dans la v5, delay.tex et la lettre, notamment
(delay.tex) l. 425, 458, 460, 572, 574, 576, 578, 582, 588, 626, 628, et les figures 1,
2B, 6 (panneau des effets), 7 et 8.

- **Interaction cohérence × durée sur le taux pendant le délai.** delay.tex (l. 588, et
  Discussion l. 667) la présente comme « the paper's clearest disagreement with its own
  model ». Avec 40 000 essais et l'ancien critère, chacun des deux modèles la reproduit
  avec quelques jeux (5 et 8). À réécrire selon les résultats de la relance.
