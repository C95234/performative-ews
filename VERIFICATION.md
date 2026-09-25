# Vérification indépendante

Ce document rapporte une vérification indépendante du modèle et des scripts
de ce dépôt, menée par un second système (Claude Code) sur demande de
l'auteur — **pas une relecture par un chercheur extérieur au domaine**, qui
reste à faire (voir README, section Limites). L'objectif : écrire des
ré-implémentations *from scratch* à partir des seules équations publiées
(jamais en copiant le code fourni), avec des graines aléatoires disjointes,
pour vérifier que les résultats numériques sont réels et non un artefact
d'un script particulier — dans l'esprit de la revue qui a déjà été faite
sur la contribution [`ewstools`](https://github.com/ThomasMBury/ewstools/pull/482)
associée à ce projet (comparaison à une implémentation indépendante,
recherche de bistabilité cachée, tests de sensibilité).

Tout le code de vérification est dans [`verification/`](verification/).
Aucun fichier de ce dossier n'a été utilisé pour produire les résultats du
document lui-même — c'est un contrôle a posteriori, distinct.

## 1. Équation de point fixe (`verification/check_fixed_point.py`)

Ré-implémentation indépendante de $\lambda_{\text{eff}} = \lambda(t) + \kappa \cdot A(\sigma^2/2|\lambda_{\text{eff}}|)$, résolue par un scan fin + bissection (pas la bissection simple du script fourni).

- **Tableaux du document reproduits exactement** : $\kappa=-0{,}1$ donne $\lambda_{\text{eff}}=-0{,}1994/-0{,}1100/-0{,}1010$ pour $\lambda(t)=-0{,}100/-0{,}010/-0{,}001$ (document : $-0{,}199/-0{,}110/-0{,}101$) ; $\kappa=+0{,}05$ confirme la disparition du point fixe entre $\lambda(t)=-0{,}100$ et $-0{,}050$.
- **Pas de bistabilité cachée** : scan à 21 points entre $\lambda=-0{,}06$ et $-0{,}05$ pour $\kappa=+0{,}05$ — une seule racine partout, disparition nette à $\lambda=-0{,}0500$, jamais deux solutions simultanées.
- **Constat supplémentaire, non signalé dans le document** : $\lambda_c(\kappa) \approx -\kappa$ quasiment exactement pour les petits $\kappa$ (ex. $\kappa=0{,}05 \Rightarrow \lambda_c \approx -0{,}0496$ ; $\kappa=0{,}10 \Rightarrow \lambda_c \approx -0{,}0996$).

## 2. Simulation stochastique complète (`verification/check_simulation.py`)

Ré-implémentation indépendante de la SDE (Euler-Maruyama), variance glissante recalculée à **chaque pas** (le script fourni ne la recalcule que tous les 20 pas), graines `seed0=10000+` (jamais utilisées par les scripts fournis, qui utilisent `seed=42`).

| $\kappa$ | Document (seed=42) | Vérification indépendante (seed≠42) |
|---|---|---|
| 0 (classique) | 63,6% | 61,6% |
| −0,05 | 23,2% | 30,8% |
| −0,10 | 11,6% | 9,6% |
| −0,20 | 2,6% | 1,4% |
| +0,02 | 77,0% | 77,2% |
| +0,05 | 89,2% | 87,2% |
| +0,08 | 95,8% | 95,8% |

Monotonie stricte confirmée sur les deux jeux de graines. Les écarts (quelques points de %) sont cohérents avec le bruit d'échantillonnage Monte-Carlo attendu sur 500 réalisations, pas un signe de divergence entre implémentations.

## 3. Extension Kuramoto (`verification/check_kuramoto.py`)

Même exercice sur le modèle de synchronisation, graines `seed0=20000+`.

| $\kappa$ | Document (seed=42) | Vérification indépendante |
|---|---|---|
| 0 | 100% / t=141,7 | 100% / t=144,9 |
| −0,5 | 100% / t=273,0 | 100% / t=273,6 |
| −1,0 | 99,0% / t=402,1 | 99,2% / t=403,0 |
| −2,0 | 0% | 0% |
| +0,5 | 100% / t=50,5 | 100% / t=50,8 |
| +1,0 | 100% / t=22,5 | 100% / t=23,1 |
| +2,0 | 100% / t=12,2 | 100% / t=12,6 |

Accord très net, y compris sur le plancher de stabilité à $\kappa=-2{,}0$ (0/500 des deux côtés).

**Test de sensibilité supplémentaire** (seuil de bascule $r_{\text{tip}} \in \{0{,}6; 0{,}7; 0{,}8\}$, non fait dans le document) : à $r_{\text{tip}}=0{,}6$, $\kappa=-2{,}0$ bascule tout de même 21,5% du temps — ce qui **confirme** la mise en garde déjà présente dans le document ("0/30 bascule... pas que la bascule est impossible dans l'absolu").

## 4. Validité de l'approximation quasi-stationnaire (AQS) — nouveau, non fait dans le document

Le document suppose $V_t \approx \sigma^2/(2|\lambda_{\text{eff}}|)$ à chaque instant (séparation d'échelles de temps). Vérifié ici directement : simulation de **20 000 réalisations en parallèle**, variance d'ENSEMBLE à travers les réalisations à chaque instant (un estimateur non biaisé de $\mathrm{Var}(y_t)$, plutôt que la fenêtre glissante d'une seule trajectoire, qui mélange biais de mesure et validité de l'AQS), comparée à la prédiction quasi-stationnaire.

**Sans rétroaction ($\kappa=0$, `check_qsa.py`)** : l'AQS tient à moins de 3% d'écart de $t=20$ à $t=400$ ($\lambda$ de $-0{,}96$ à $-0{,}20$), puis décroche progressivement : 6,5% à $\lambda=-0{,}12$, 11% à $\lambda=-0{,}08$, **27,4% à $\lambda=-0{,}04$**. La variance réelle est systématiquement **inférieure** à la prédiction quasi-stationnaire dans cette zone (elle n'a pas le temps de "rattraper" l'équilibre théorique qui croît de plus en plus vite).

**Avec rétroaction (`check_qsa_feedback.py`)**, variance d'ensemble sur la dynamique réelle (chaque trajectoire réagit à sa propre fenêtre glissante) comparée à la prédiction du point fixe :

| | $\kappa=-0{,}10$ (auto-invalidant) | $\kappa=+0{,}05$ (auto-réalisateur) |
|---|---|---|
| Écart à $\lambda_{\text{bare}}=-0{,}28$ | 2,6% | 2,0% |
| Écart à $\lambda_{\text{bare}}=-0{,}20$ | 8,4% | 18,3% |
| Écart à $\lambda_{\text{bare}}=-0{,}16$ | 18,2% | 20,9% |
| Écart à $\lambda_{\text{bare}}=-0{,}12$ | 29,7% | 27,0% |
| Écart à $\lambda_{\text{bare}}=-0{,}08$ | 31,0% | 47,5% |

Ici, contrairement au cas sans rétroaction, la variance réelle **dépasse** la prédiction pour $\kappa=-0{,}10$ (le plateau performatif s'installe avec un dépassement transitoire avant de se stabiliser). Pour $\kappa=+0{,}05$, le décrochage est plus rapide et plus sévère : à $\lambda_{\text{bare}}=-0{,}04$ (où le point fixe a déjà officiellement disparu, $\lambda_c \approx -0{,}0496$), la variance d'ensemble a déjà explosé (21,0, contre 12,5 pour le cas sans rétroaction au même point).

### Ce que ça implique, honnêtement

- **La validation numérique principale du document (§5, simulation stochastique complète) n'est pas affectée** : elle a toujours été faite par simulation complète, pas par l'AQS elle-même — c'est exactement pour ça que le document ne se contente pas de l'analyse du point fixe. Les sections 2 et 3 ci-dessus confirment cette validation indépendamment.
- **Ce qui doit être nuancé** : l'analyse du point fixe (section "Analyse du point fixe" du document, le résultat en deux régimes avec un $\lambda_c(\kappa)$ précis) est une **bonne image qualitative**, mais **quantitativement imprécise dans le dernier quart de l'approche vers $\lambda=0$** — exactement la zone où se joue le phénomène qu'elle décrit. Le décrochage est plus rapide et plus sévère côté auto-réalisateur ($\kappa>0$) que côté auto-invalidant ($\kappa<0$), une asymétrie qui n'était pas documentée.
- Recommandation pour une future version : présenter $\lambda_c(\kappa)$ comme un ordre de grandeur utile, pas une prédiction quantitative exacte du moment de bascule — la simulation complète reste la seule source fiable pour un chiffre précis, ce que le document fait déjà pour ses résultats principaux.

## Ce qui n'a pas été vérifié ici

- La revendication de nouveauté bibliographique (aucun outil de vérification par le code ne peut trancher ça).
- Une relecture par un chercheur extérieur au domaine — cette vérification reste interne au processus de production du document, pas une évaluation par les pairs.
