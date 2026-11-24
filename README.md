# Métaheuristiques — Set Packing Problem (SPP)

Projet réalisé dans le cadre du module **Métaheuristiques** — Master Informatique, Nantes Université (2026–2027).

Ce dépôt regroupe des implémentations en **Julia** de méthodes de résolution du *Set Packing Problem* (SPP), ainsi que des instances et des résultats expérimentaux.

## 1. Problème étudié

Le **Set Packing Problem pondéré** consiste à sélectionner un ensemble de colonnes afin de maximiser leur valeur totale, tout en respectant des contraintes d'incompatibilité.

La formulation mathématique est :

$$
\max \; z = \sum_{j=1}^{n} c_j x_j
$$

sous les contraintes :

$$
\sum_{j=1}^{n} a_{ij}x_j \leq 1 \quad \forall i \in \{1,\ldots,m\},
\qquad x_j \in \{0,1\}.
$$

Avec :

- $c_j$ : valeur associée à la colonne $j$ ;
- $a_{ij}$ : coefficient de la matrice binaire des contraintes $A$ ;
- $x_j$ : variable de décision indiquant si la colonne $j$ est sélectionnée ;
- $m$ et $n$ : nombres de contraintes et de colonnes.

Deux colonnes contenant un `1` sur une même ligne de $A$ ne peuvent pas être sélectionnées simultanément.

## 2. Organisation du projet

```text
BenlaroussiMiskarEI/
├── livrableEI1.jl       # Algorithmes et expérimentations EI1
├── src/
│   ├── loadSPP.jl       # Chargement des instances
│   ├── setSPP.jl        # Utilitaires SPP
│   └── getfname.jl      # Gestion des noms de fichiers
├── dat/                 # Instances au format .dat
├── res/                 # Résultats des expérimentations
├── doc/                 # Documentation et rapport
└── README.md
```

## 3. Exercice d'implémentation 1 (EI1)

### 3.1 Construction gloutonne

L'algorithme glouton construit progressivement une solution réalisable. Il démarre avec une solution vide et choisit à chaque itération une colonne disponible selon le critère d'utilité :

$$
u_j = \frac{c_j}{d_j},
$$

où $d_j$ est le nombre de contraintes auxquelles participe la colonne $j$. Pour une colonne sans contrainte, l'implémentation utilise directement sa valeur $c_j$.

Après chaque sélection, les colonnes incompatibles sont retirées de l'ensemble des candidats. Deux vecteurs booléens sont utilisés :

- `solution` indique les colonnes sélectionnées ;
- `mask` indique les colonnes encore disponibles.

La construction s'arrête lorsqu'il n'existe plus de candidat utile. Cette stratégie est rapide, mais ne garantit pas l'optimalité globale.

### 3.2 Recherche locale : échange 1-1

La recherche locale démarre à partir de la solution gloutonne et explore des voisins par **échange 1-1** : retirer une colonne sélectionnée et en ajouter une autre.

Chaque candidat est accepté uniquement s'il respecte les contraintes $Ax \leq \mathbf{1}$ et améliore strictement la fonction objectif.

La méthode suit une stratégie de **première amélioration** (*first improvement*). Les échanges sont répétés jusqu'à ce qu'aucun voisin améliorant ne soit trouvé. Le résultat est un optimum local pour le voisinage considéré.

### 3.3 Résolution exacte : JuMP et GLPK

Une résolution de référence est réalisée avec :

- **JuMP**, pour modéliser les variables binaires, l'objectif et les contraintes ;
- **GLPK** (*GNU Linear Programming Kit*), pour résoudre le programme linéaire en nombres entiers.

Le statut `OPTIMAL` indique que l'optimalité a été prouvée. Le statut `TIME_LIMIT` indique que le solveur s'est arrêté à la limite de temps sans preuve d'optimalité ; la valeur retournée, lorsqu'elle existe, est alors la meilleure solution réalisable trouvée par le solveur.

Lorsque l'optimum $z^*>0$ est prouvé, l'écart relatif de l'heuristique est calculé par :

$$
\operatorname{Gap}(\%) = 100\,\frac{z^* - z_{\mathrm{heuristique}}}{z^*}.
$$

## 4. Installation et exécution

### Prérequis

- Julia 1.x ;
- Packages `JuMP` et `GLPK`.

Installation des packages dans le REPL Julia :

```julia
import Pkg
Pkg.add(["JuMP", "GLPK"])
```

Depuis la racine du dépôt, charger le programme :

```julia
include("livrableEI1.jl")
```

Résoudre une instance :

```julia
resoudreSPP("dat/didactic.dat")
```

Exécuter plusieurs instances :

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

Lancer les expériences sur les instances détectées dans `dat/` :

```julia
experimentationSPP()
```

Options disponibles :

```julia
# Exécution sans solveur exact
experimentationSPP(fichiers; exact=false)

# Limite de temps de 180 secondes pour GLPK, par instance
experimentationSPP(fichiers; limite_secondes=180.0)
```

Les résultats sont enregistrés dans **`res/resultats_EI1.csv`** (séparateur `;`).

## 5. Résultats expérimentaux — EI1

Résultats obtenus lors des premiers essais sur cinq instances :

| Instance | Glouton | Recherche locale | GLPK | Statut GLPK | Gap local |
|---|---:|---:|---:|---|---:|
| `didactic.dat` | 30 | 30 | 30 | `OPTIMAL` | 0,00 % |
| `pb_100rnd0100.dat` | 342 | 343 | 372 | `OPTIMAL` | 7,80 % |
| `pb_100rnd0700.dat` | 495 | 495 | 503 | `OPTIMAL` | 1,59 % |
| `pb_200rnd0100.dat` | 351 | 365 | 395 | `TIME_LIMIT` | — |
| `pb_500rnd0100.dat` | 285 | 285 | 262 | `TIME_LIMIT` | — |

La recherche locale améliore certaines solutions gloutonnes, mais peut rester bloquée dans un optimum local. Pour les deux dernières instances, les valeurs GLPK ne constituent pas des optima prouvés.

Les temps enregistrés sont des **durées écoulées** en secondes. Le temps de recherche locale est mesuré séparément de celui de la construction gloutonne.

## 6. Exercice d'implémentation 2 (EI2)

*Section réservée aux développements de l'EI2.*

### GRASP

À documenter après implémentation : construction gloutonne randomisée, liste restreinte de candidats (RCL), recherche locale et critères d'arrêt.

### Composant additionnel

À documenter : méthode retenue, principe, paramètres et intégration avec GRASP.

### Expérimentations et résultats

À compléter : protocole expérimental, comparaison avec l'EI1, tableaux de résultats et analyse.

---

**Technologies :** Julia · JuMP · GLPK
