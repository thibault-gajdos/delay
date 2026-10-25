# Changements depuis la v5 (`Achard et al R1_v5.docx`)

État au 25 octobre 2026. Source : `manuscript/latex/delay.tex`.

La v5 contenait le front matter, l'Introduction et les Methods. Ces parties sont reprises dans `delay.tex`, avec les changements ci-dessous. Les Results, la Discussion, les figures, les tables et l'Annexe viennent de `delay.tex` et ont été mis à jour.

## Termes

- « covert muscle activations » → « partial muscle activations (PMAs) ».
- « premotor time » et « motor time » abrégés en PMT et MT.
- Taux de retour de la fluctuation de maintien : κ_ε partout, Annexe comprise.
- Un seul nom par objet : « maintenance target » (m), « signal reaching the gate » (ỹ), « slope » (pente), « parameter vector ».
- Fenêtres : « stimulus-onset window », « delay window », « RS window », « response window ».

## Introduction

- Citation de van Ede & Nobre (2023) rétablie (demandée par l'éditeur ; retirée dans la v5).

## Methods

**Theoretical framework**
- L'équation (11), ỹ après le RS, est numérotée.
- Coquille corrigée : ε(τ) → ε(t) dans « the maintenance-noise process ε(t) is discontinued ».

**Evaluation over the parameter space**
- 4 000 → 40 000 essais par condition et par modèle.
- Note de bas de page 1 supprimée. Elle est remplacée par la règle appliquée aux simulations : une activation commence chaque fois que la commande neuronale devient positive ; la dernière avant la réponse est l'activation liée à la réponse ; toutes les autres sont des PMA ; aucun autre critère de détection.
- Admissibilité :
  - « responses before the RS ≤ 10 % » est remplacé par « the response-related muscle activation begins before the RS on ≤ 10 % of trials (this includes every response before the RS) » ;
  - critère d'exactitude (> 80 % à 50 % de cohérence) supprimé ;
  - nombre de vecteurs admissibles ajouté : 228 (Modèle 1) et 204 (Modèle 2) sur 1 200.
- Plages des paramètres :
  - k et λ couvrent environ un facteur deux de part et d'autre des estimations de Dendauw et al. (2024) ;
  - r et g sont échantillonnés sur des plages plus larges, g au-dessus de l'estimation de Dendauw et al.
- Ajout de deux phrases sur les figures :
  - la Figure 1 montre trente vecteurs admissibles dans les deux modèles ;
  - la Figure 2B est tracée au vecteur représentatif, avec la règle de choix.
- Fuite :
  - v5 : un même taux γ ∈ {0,25 ; 0,5 ; 1} s⁻¹ pour les deux accumulateurs ;
  - maintenant : γ_x et γ_q varient indépendamment sur {0 ; 0,25 ; 0,5 ; 1} s⁻¹, soit 15 combinaisons en plus de (0, 0), sur les 304 vecteurs admissibles dans au moins un modèle, avec les mêmes paramètres et le même nombre d'essais.

**Statistical analysis**
- Renvoi aux Tables S1–S4 ajouté.

**Model–data comparison**
- Les vingt effets comparés sont ceux du Table 2.
- Deux vecteurs du Modèle 1 sur 228 sont exclus des trois effets sur la pente et la courbure pendant le délai (trop peu de PMA dans une condition).

## Results (réorganisés)

- **Structure :**
  - une sous-section « Accuracy, RT, PMT and MT », puis une par fenêtre ;
  - chacune contient : Theoretical predictions, Empirical results et, pour les trois premières, Model–data comparison.
- **Règle de la v5 appliquée partout :**
  - minimum et maximum de chaque effet sur les vecteurs admissibles ;
  - nombre de vecteurs qui reproduisent, ne reproduisent pas, ou sont incertains ;
  - signes donnés seulement pour les effets dont la direction n'est pas fixée par une dérivation ;
  - le signe de l'effet de la durée sur le niveau pendant le délai est rapporté, faute de dérivation.
- **Prédictions corrigées :**
  - Effet de la durée pendant le stimulus : un essai encore en accumulation ne franchit pas la porte. La différence entre durées est la différence de deux termes, et son signe sous le Modèle 1 est établi numériquement (227 vecteurs sur 228).
  - L'effet de la cohérence sur la pente pendant le stimulus n'a pas de signe fixé ; seule la prédiction sur le niveau est maintenue.
  - Pendant le délai, la pente du log-taux est proportionnelle à la distance à la porte moyennée en pondérant chaque essai par son taux.
  - Condition de convexité sous le Modèle 2 : la moyenne de a² pondérée par le taux dépasse 1.
  - Interaction cohérence × durée pendant le délai : signe non fixé par la dérivation.
- **Erreur de l'approximation stationnaire :** 3 × 10⁻³ (Modèle 1) et 1,2 × 10⁻³ (Modèle 2) en relatif, sur la fenêtre (au moins 0,8 s après le début du maintien).
- **MT :** à de nombreux vecteurs, les modèles produisent un MT plus court après le stimulus court. L'écart vient d'essais où l'activation liée à la réponse commence avant le RS, plus fréquents après le stimulus long. Dans les données, ces essais sont presque absents.
- **Fenêtres RS et réponse :** PMA et MBL décrits sans test ; LRP testé sur son niveau et sa pente ; pas de comparaison avec les modèles.
- **Ajout :** l'effet de la durée sur la courbure du LRP pendant le stimulus est rapporté (p = .014). Le texte précédent disait « aucun effet significatif ».
- **MBL pendant le stimulus :** « built up toward stronger effector-selective preparation », au lieu de « rose ». L'estimation de la pente est négative.

## Discussion

- La correspondance avec Inagaki et al. (2022) et Churchland & Shenoy (2024) est maintenant dans la Discussion (elle était dans l'Introduction de `delay.tex`, absente de la v5).
- Effet de la durée pendant le stimulus : 227 vecteurs du Modèle 1 sur 228 donnent le signe opposé ; aucun ne reproduit l'effet.
- Interaction cohérence × durée pendant le délai : le signe est opposé pour la plupart des vecteurs des deux modèles (206 sur 228, 188 sur 204), et 10 vecteurs par modèle reproduisent l'effet. C'est un désaccord avec les deux modèles.
- Les deux fuites (γ_x, γ_q) remplacent la fuite unique.

## Annexe

- κ → κ_ε.
- Proposition A.1 : la fenêtre du délai commence au moins 0,8 s après le début du maintien.
- Canal non sélectionné : énoncé exact.
- Corollaire A.4 : ajout de l'énoncé au niveau de la condition. La vitesse standardisée est notée b (et non k, déjà pris par le paramètre de cohérence).
- Proposition A.5 remplacée par une décomposition exacte ; le signe est établi numériquement.
- Remarque A.6 : les deux fuites.
- Remarque A.7 : le premier point porte maintenant sur la différence entre simulations (chaque franchissement compte) et données (fusion à 25 ms).

## Figures et tables

- Table 1 : κ_ε ; lignes ajoutées pour γ_x, γ_q et ỹ.
- Table 2 (nouveau) : comparaison modèles–données, générée par `scripts/table_comparison.py`.
- Tables S1–S4 générées par les scripts. Pentes et courbures en unités physiques (s⁻², s⁻³ pour les PMA).
- Figures 1, 2B, 3, 5, 6 et 7 régénérées.
  - Figure 6 : 20 panneaux.
  - Figure 7 : fenêtre à partir de 0,8 s.
  - Légendes mises à jour.

## Bibliographie

- Ajouts : Gold & Shadlen (2000) ; Servant, White, Montagnini & Burle (2015).
- Rogge, Jocham & Ullsperger (2022) : entrée corrigée (auteurs et année manquaient).

## Front matter

- 7 figures (+ 1 supplémentaire), 2 tables (+ 5 supplémentaires).
- 51 pages.
- Mots : résumé 248, introduction 644, discussion 1302, citations comprises.
