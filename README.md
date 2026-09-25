# Signaux précurseurs performatifs

**Auteur :** Clément Marin — septembre 2026
**Statut : premier brouillon, non relu par un tiers extérieur au domaine, cohérence mathématique en cours de vérification.** Ceci n'est pas présenté comme une théorie aboutie de la réflexivité sociale, ni comme une découverte, mais comme un premier modèle jouet, dérivé et vérifié numériquement, destiné à servir de point de départ à une discussion ouverte.

## En une page

Les indicateurs de ralentissement critique (variance, autocorrélation) supposent un paramètre de contrôle $\lambda(t)$ purement exogène — le système observé ne « lit » jamais l'indicateur calculé sur lui-même. Cette hypothèse est raisonnable pour un lac ou un plasma ; elle est fausse par construction pour un système social observé par un outil de détection : si l'alerte est publiée, les acteurs peuvent y réagir, et cette réaction change la trajectoire même que l'indicateur cherchait à décrire.

Ce document pose un premier modèle formel reliant deux champs établis séparément — les signaux précurseurs de bascule et la **prédiction performative** (Perdomo, Zrnic, Mendler-Dünner & Hardt, 2020) — qu'une revue bibliographique (recherche multi-sources, puis recoupement direct des bibliographies des deux synthèses les plus récentes de chaque champ) n'a trouvés combinés formellement nulle part ailleurs à ce jour. Cette revue s'appuie sur une recherche documentaire, pas sur une base de citations exhaustive — l'absence de résultat trouvé n'est pas une preuve d'absence définitive.

Le postulat : un détecteur qui publie une alerte $a_t = A(V_t)$ (fonction sigmoïde de la variance glissante $V_t$) modifie en retour le paramètre de stabilité effectif du système, $\lambda_{\text{eff}} = \lambda(t) + \kappa \cdot A(V_t)$. Selon le signe de $\kappa$ (coefficient de rétroaction performative), l'équation de point fixe qui en résulte se comporte de deux façons radicalement différentes : pour $\kappa<0$ (réponse auto-invalidante), $\lambda_{\text{eff}}$ reste borné loin de zéro — la variance plafonne au lieu de diverger ; pour $\kappa>0$ (réponse auto-réalisatrice), le point fixe cesse d'exister avant que le système « nu » n'atteigne la bifurcation — la bascule survient plus tôt.

Ce résultat est vérifié par simulation stochastique complète (Euler-Maruyama, variance glissante calculée sur une fenêtre finie plutôt que par une formule fermée) : 500 réalisations indépendantes par valeur de $\kappa$ sur un modèle de bifurcation nœud-col, et une extension au modèle de synchronisation de Kuramoto pour vérifier que le mécanisme ne dépend pas du modèle sous-jacent particulier.

## Limites assumées

- Un seul mécanisme de rétroaction est modélisé (réponse linéaire du paramètre de stabilité à une alerte sigmoïde). Le masquage stratégique (loi de Goodhart) et l'aléa moral (comportement plus risqué en l'*absence* d'alerte) ne sont pas modélisés ici.
- Le détecteur ne se recalibre pas — pas de boucle de ré-entraînement au sens de la « stabilité performative » de Perdomo et al. Ce modèle teste un $\kappa$ fixé, pas la convergence d'un processus d'apprentissage.
- Validé sur deux modèles simulés (nœud-col, Kuramoto), à 500-1000 réalisations, avec vérification de sensibilité aux paramètres de la fonction d'alerte et au seuil de détection de bascule — mais **aucun test sur données réelles**, **aucune relecture par un chercheur extérieur au domaine à ce stade**, et une vérification bibliographique par recoupement de synthèses plutôt que par une base de citations exhaustive.
- Premier modèle jouet, pas une théorie de la réflexivité sociale.

Le détail complet (dérivations, tableaux de résultats, code annoté) est dans [`performative_ews.pdf`](performative_ews.pdf) / [`performative_ews.tex`](performative_ews.tex). Les scripts de validation (`performative_ews_validation_500_1000.py`, `performative_ews_kuramoto_extension.py`) reproduisent exactement les chiffres cités dans le document.

## Précurseurs directs

Aucun des travaux ci-dessous ne combine formellement les deux champs — chacun se situe d'un seul côté du pont que ce document tente d'esquisser :

- **Perdomo, Zrnic, Mendler-Dünner & Hardt (2020)**, *ICML* — définit la prédiction performative ; $\kappa$ ici joue le rôle de la sensibilité $\varepsilon$ de leur cadre, transposée au paramètre de stabilité d'un processus stochastique plutôt qu'à une distribution de caractéristiques.
- **Hardt & Mendler-Dünner (2025)**, *Statistical Science* — synthèse la plus récente du champ, utilisée pour le recoupement bibliographique.
- **Diekert, Heyen, Nesje & Shayegh (2025)**, *J. R. Soc. Interface* — le lien formel le plus proche : un système d'alerte précoce change les décisions d'un planificateur, avec un comportement parfois plus risqué en l'absence d'alerte. Un seul décideur rationnel, pas de population réactive, pas de statistique de ralentissement critique.
- **Bauch, Sigdel, Pharaon & Anand (2016)**, *PNAS* — la rétroaction humaine peut atténuer le signal précurseur dans un modèle couplé forêt-opinion.
- **Diks, Hommes & Wang (2019)**, *Empirical Economics* — notent, en une phrase, qu'une crise prédite peut déclencher une réaction de panique et ainsi survenir plus tôt, sans le formaliser.
- **Sornette (2003)**, *Physics Reports* — liste, de façon informelle, les scénarios où une prédiction publique de krach s'auto-réalise ou s'auto-invalide.
- **Drehmann & Juselius (2014)**, BIS / *Int. J. Forecasting* — notent qu'un indicateur d'alerte utilisé en politique publique devient sujet à la critique de Lucas.
- **Scheffer et al. (2009)**, *Nature* — synthèse de référence des signaux précurseurs classiques (sans rétroaction).

## Prochaines étapes possibles

Listées dans le document (section finale) : formaliser le mécanisme de masquage stratégique comme second modèle ; étudier la convergence d'un détecteur qui recalibre son seuil à partir de données déjà influencées par ses propres alertes passées ; tester sur un jeu de données réel où la réflexivité est plausible (marchés financiers). La question de l'hébergement (extension d'un projet existant, ou projet et interface séparés) reste ouverte.

## Contact

Retours et critiques bienvenus — via une issue sur ce dépôt.
