# Métaheuristiques — Set Packing Problem (SPP)

**Nantes Université · Master Informatique · 2026–2027**  
**Projet :** exercices d'implémentation (EI1, puis EI2) en Julia  
**Dépôt :** [BenlaroussiMiskarEI](https://github.com/Elias8674/BenlaroussiMiskarEI)

> **État du projet :** EI1 implémenté et testé sur **5 instances** à ce stade. L'expérimentation demandée sur **au moins 10 instances variées** et la finalisation du rapport restent à effectuer avant le dépôt. EI2 : section réservée en fin de document.

## 1. Objectif du projet

Résoudre le **Set Packing Problem pondéré (WSPP)** : sélectionner des colonnes de valeur positive ou nulle, sans que deux colonnes sélectionnées utilisent une même ressource (une même ligne de contrainte), afin de maximiser la somme de leurs valeurs.

Le modèle mathématique est :

\[
\max z = \sum_{j=1}^{n} c_j x_j
\]

sous contraintes :

\[
\sum_{j=1}^{n} a_{ij}x_j \leq 1 \quad \forall i=1,\ldots,m, \qquad x_j\in\{0,1\}.
\]

- `C` : vecteur des valeurs `c_j` des `n` colonnes ;
- `A` : matrice d'incidence binaire de taille `m × n` ;
- `x[j] = 1` : colonne `j` sélectionnée ; `x[j] = 0` : non sélectionnée ;
- chaque ligne de `A` impose qu'au plus une des colonnes portant un `1` sur cette ligne soit choisie.

**Exemple :** si une ligne contient `[1, 1, 0]`, choisir simultanément les colonnes 1 et 2 donnerait `1 + 1 > 1` et violerait la contrainte.

## 2. Organisation du dépôt

```text
BenlaroussiMiskarEI/
├── livrableEI1.jl       # Point d'entrée de l'EI1
├── src/
│   ├── loadSPP.jl       # Lecture des instances SPP
│   ├── setSPP.jl        # Utilitaires fournis pour le SPP
│   └── getfname.jl      # Utilitaires liés aux noms de fichiers
├── dat/                 # Instances .dat
├── res/                 # Résultats expérimentaux (CSV)
├── doc/                 # Rapport LaTeX et PDF
└── README.md
```

Les fichiers d'instances sont au format attendu par `loadSPP.jl` (références SPP de type OR-Library). Le chargeur fournit `C, A = loadSPP(fname)` ; il faut conserver le format compatible avec ce chargeur. La signification exacte des champs du fichier `.dat` dépend de son implémentation.

### Instances utilisées lors des premiers essais

- `didactic.dat` : 7 contraintes, 9 colonnes ;
- `pb_100rnd0100.dat` : 500 contraintes, 100 colonnes ;
- `pb_100rnd0700.dat` : 100 contraintes, 100 colonnes ;
- `pb_200rnd0100.dat` : 1 000 contraintes, 200 colonnes ;
- `pb_500rnd0100.dat` : 2 500 contraintes, 500 colonnes.

**Important :** ces **5** instances constituent des tests préliminaires ; elles ne suffisent pas à satisfaire l'exigence des **10 instances variées** de l'EI1. D'autres fichiers peuvent déjà être présents dans `dat/`, mais ils ne figurent pas dans les résultats présentés ici.

## 3. Méthodes implémentées pour l'EI1

### 3.1 Construction gloutonne

L'heuristique construit une solution réalisable à partir d'une solution vide.

1. Initialiser `solution = falses(n)` (aucune colonne sélectionnée) et `mask = trues(n)` (toutes les colonnes disponibles).
2. Pour chaque colonne disponible, calculer une **utilité** : `C[j] / d[j]`, où `d[j]` est le nombre de `1` dans la colonne `j` de `A` ; une colonne sans contrainte reçoit, par convention dans ce code, l'utilité `C[j]`.
3. Choisir la colonne disponible de plus grande utilité (`argmax`).
4. L'ajouter à `solution`, actualiser le score et rendre indisponibles dans `mask` toutes les colonnes partageant au moins une contrainte avec elle.
5. Répéter jusqu'à épuisement des candidats utiles.

**Différence fondamentale :** `solution[j]` indique si la colonne est *choisie* ; `mask[j]` indique si elle est *encore autorisée*. Une colonne peut être non choisie et néanmoins interdite à cause d'un conflit.

Ce critère est **heuristique** : il privilégie un bon rapport valeur/contraintes mobilisées, mais ne garantit pas l'optimum global. Les variantes du code doivent éviter de sélectionner un candidat interdit (par exemple en lui attribuant une utilité `-Inf`) et traiter les cas limites.

### 3.2 Recherche locale — échange 1-1

La recherche locale part de la solution gloutonne et examine des voisins obtenus par un **échange 1-1** :

- retirer une colonne actuellement sélectionnée ;
- ajouter une colonne actuellement non sélectionnée ;
- vérifier la faisabilité `A * x .<= 1` ;
- accepter le voisin si son score est **strictement supérieur**.

La fonction de descente adopte la stratégie **première amélioration** (*first improvement*) : elle retourne dès qu'elle trouve un échange admissible améliorant. La recherche locale relance cette descente jusqu'à ce qu'aucun échange 1-1 améliorant ne soit trouvé.

On obtient alors un **optimum local relativement au voisinage 1-1**, et non nécessairement l'optimum global. Un mouvement d'ajout seul ou un échange 1-2 pourrait parfois produire une amélioration inaccessible au seul 1-1.

### 3.3 Résolution exacte — JuMP + GLPK

- **JuMP** est l'interface Julia permettant de décrire le modèle mathématique : variables binaires, fonction objectif et contraintes.
- **GLPK** (*GNU Linear Programming Kit*) est le solveur qui traite le modèle de programmation linéaire en nombres entiers.
- Le statut `OPTIMAL` signifie que le solveur a **prouvé** l'optimalité de la valeur trouvée.
- Le statut `TIME_LIMIT` signifie qu'il a atteint la limite de temps **sans preuve d'optimalité**. Une valeur affichée dans ce cas est une solution réalisable trouvée pendant la recherche, **pas forcément l'optimum**.

La comparaison avec GLPK permet d'évaluer la qualité des heuristiques. Si un optimum positif `z*` est prouvé, l'écart relatif de la recherche locale se calcule par :

\[
\operatorname{gap}(\%)=100\frac{z^*-z_{\mathrm{local}}}{z^*}.
\]

Sans optimum prouvé, le champ `gap_pct` est laissé à `missing` ; on ne doit pas calculer un « écart à l'optimum » à partir d'une valeur `TIME_LIMIT`.

## 4. Installation et exécution

### Prérequis

- Julia 1.x ;
- packages Julia `JuMP`, `GLPK` et leurs dépendances ;
- les fichiers utilitaires dans `src/` et les instances dans `dat/`.

Installer les dépendances dans le REPL Julia (une seule fois dans l'environnement concerné) :

```julia
import Pkg
Pkg.add(["JuMP", "GLPK"])
```

Le code utilise aussi des bibliothèques standard de Julia telles que `LinearAlgebra` et `Printf`.

Depuis la racine du projet, démarrer Julia, puis :

```julia
include("livrableEI1.jl")
resoudreSPP("dat/didactic.dat")
```

La fonction affiche les scores, solutions, temps et statut du solveur ; elle retourne également une structure contenant les résultats.

### Expérimenter sur plusieurs instances

Pour choisir explicitement les fichiers :

```julia
fichiers = [
    "dat/didactic.dat",
    "dat/pb_100rnd0100.dat",
    "dat/pb_100rnd0700.dat",
    "dat/pb_200rnd0100.dat",
    "dat/pb_500rnd0100.dat"
]

experimentationSPP(fichiers)
```

La version de travail de `experimentationSPP` doit parcourir **`for fichier in fichiers`** et non `fichiers[1:10]` pour accepter les tests de moins de 10 fichiers. La contrainte pédagogique des 10 instances s'applique à l'**expérimentation finale**, pas à chaque essai de développement.

Pour l'expérimentation finale, fournir **au moins 10 chemins** valides et variés. Vérifier la version effectivement enregistrée du code avant l'exécution : certaines versions intermédiaires imposaient 10 fichiers et tronquaient la liste.

Paramètres optionnels :

```julia
# Test rapide des seules heuristiques
experimentationSPP(fichiers; exact=false)

# Accorder jusqu'à 180 secondes à GLPK par instance
experimentationSPP(fichiers; limite_secondes=180.0)
```

**Sortie :** `res/resultats_EI1.csv` (séparateur `;`). Une nouvelle expérimentation peut écraser ce fichier : sauvegarder les résultats à conserver.

## 5. Premiers résultats expérimentaux (5 instances)

Les résultats ci-dessous proviennent d'une exécution locale de développement. Ils sont **indicatifs** et ne remplacent pas les expériences finales reproductibles.

| Instance | Glouton | Recherche locale | GLPK (meilleure valeur trouvée) | Statut GLPK | Gap local à l'optimum |
|---|---:|---:|---:|---|---:|
| `didactic.dat` | 30 | 30 | 30 | `OPTIMAL` | 0,00 % |
| `pb_100rnd0100.dat` | 342 | 343 | 372 | `OPTIMAL` | 7,80 % |
| `pb_100rnd0700.dat` | 495 | 495 | 503 | `OPTIMAL` | 1,59 % |
| `pb_200rnd0100.dat` | 351 | 365 | 395 | `TIME_LIMIT` | Non disponible |
| `pb_500rnd0100.dat` | 285 | 285 | 262 | `TIME_LIMIT` | Non disponible |

**Lecture :** la recherche locale améliore notamment 342 → 343 et 351 → 365 ; elle n'améliore pas les trois autres solutions sur ces essais. Pour `pb_500rnd0100.dat`, GLPK n'a trouvé que 262 dans le temps imparti alors que l'heuristique atteint 285 : **ce n'est pas une contradiction**, car GLPK n'a pas prouvé l'optimum.

### Mesures de temps

Le code de travail utilise `time_ns()` pour mesurer des **durées écoulées** (en secondes), et non une mesure stricte du **temps CPU**. `t_local_s` correspond à la phase de recherche locale seule ; le temps total de l'approche glouton + recherche locale est `t_glouton_s + t_local_s`.

Pour une comparaison rigoureuse : utiliser la même machine et le même environnement, faire une exécution de préchauffage pour réduire l'effet de la compilation JIT de Julia, répéter les mesures et préciser la méthode de chronométrage. Si l'énoncé exige spécifiquement `CPUt`, distinguer explicitement cette mesure du temps écoulé.

## 6. Préparation du dépôt au professeur — checklist EI1

- [ ] **Ajouter au moins 5 autres instances variées** pour atteindre **10 instances distinctes minimum** ; ne pas compter deux exécutions d'un même fichier comme deux instances.
- [ ] Vérifier que les 10 fichiers sont bien lisibles par `loadSPP` et couvrent des tailles / caractéristiques variées.
- [ ] Relancer `experimentationSPP` sur les 10 instances et conserver le CSV final dans `res/`.
- [ ] Vérifier pour chaque instance la faisabilité des solutions, les scores et l'absence d'erreurs d'exécution.
- [ ] Comparer glouton, recherche locale et GLPK ; différencier les cas `OPTIMAL` et `TIME_LIMIT`.
- [ ] Documenter la limite de temps GLPK, les temps mesurés et les conditions expérimentales ; relever le temps CPU si demandé.
- [ ] Compléter le rapport LaTeX dans `doc/` : modèle SPP, algorithmes, protocole, tableaux, analyse des écarts, limites et conclusion ; produire le PDF.
- [ ] Vérifier les commandes imposées `include("livrableEI1.jl")`, `resoudreSPP(fname)` et `experimentationSPP()` depuis la **racine** du projet.
- [ ] Tester sur **Linux** (environnement du CIE), y compris avec une instance inconnue mais au format attendu.
- [ ] Vérifier la présence des dossiers `src/`, `dat/`, `res/`, `doc/`, du README et des dépendances documentées.
- [ ] Créer l'archive ZIP demandée, avec la bonne arborescence et le bon nom, puis vérifier son contenu après extraction.
- [ ] Respecter les consignes du cours relatives à l'usage d'outils d'IA et aux déclarations éventuelles dans le rapport.

**Source des instances de référence indiquée dans l'énoncé :** http://www.emse.fr/~delorme/SetPackingFr.html (format OR-Library). Le choix des 10 instances et les résultats doivent être justifiés dans le rapport.

> **État avant remise :** ne pas présenter l'EI1 comme expérimentalement finalisé tant que la campagne des 10 instances, la validation Linux et le rapport ne sont pas terminés. Selon les consignes, le compte rendu PDF fait partie du livrable à remettre au terme de l'EI3 ; vérifier les modalités et dates exactes sur Madoc.

---

## 7. EI2 — GRASP et composant additionnel (à compléter)

**Statut :** à développer. Cette section sera mise à jour lors de l'exercice d'implémentation n°2.

### 7.1 Objectif

Étendre la construction gloutonne de l'EI1 à une **métaheuristique GRASP** (*Greedy Randomized Adaptive Search Procedure*) : construire plusieurs solutions par sélection gloutonne randomisée, puis les améliorer par recherche locale et conserver la meilleure.

### 7.2 Composant additionnel retenu

**À renseigner :** VND / Path Relinking / ReactiveGRASP / Destroy-and-Repair / autre composant autorisé par l'énoncé.

- Justification du choix : **à compléter**.
- Voisinages, structures de données et paramètres : **à compléter**.
- Réglage de `α` et protocole d'étude : **à compléter**.

### 7.3 Expérimentations EI2

**À compléter après implémentation :**

- comparer **EI1**, **GRASP seul** et **GRASP + composant additionnel** ;
- utiliser les **mêmes 10 instances variées** que pour l'EI1 ;
- définir un **budget de calcul commun** (par exemple 60 s) ;
- répéter les exécutions (par exemple 10 fois sur une même machine) ;
- rapporter `z_min`, `z_max`, `z_moy`, les temps et les paramètres utilisés ;
- analyser les résultats et recommander un réglage des paramètres.

### 7.4 Commandes EI2 et résultats

**À compléter :** noms des fichiers, fonctions publiques, exemples de commandes Julia, emplacement des résultats et lien vers le rapport EI2.

---

*Projet pédagogique — Nantes Université, module Métaheuristiques, année 2026–2027.*
