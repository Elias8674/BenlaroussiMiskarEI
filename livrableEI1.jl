# EI1 - Set Packing Problem (SPP)
# Construction gloutonne, recherche locale et comparaison exacte.

using JuMP
using GLPK
using LinearAlgebra
using Printf
const MOI = JuMP.MOI

# Chargement des instances fourni avec le projet.
include(joinpath(@__DIR__, "src", "loadSPP.jl"))

# Vérifie les dimensions et les données de l'instance.
function verifier_instance(C, A)
    m, n = size(A)
    length(C) == n || error("Le nombre de coûts ne correspond pas aux colonnes de A.")
    all(x -> x == 0 || x == 1, A) || error("La matrice A doit être binaire.")
    all(isfinite, C) || error("Les coûts doivent être finis.")
    return m, n
end

# Calcule le gain total d'une solution.
score_solution(C, x) = dot(C, Int.(x))

# Vérifie que chaque contrainte contient au plus une colonne choisie.
function est_realisable(A, x)
    size(A, 2) == length(x) || return false
    return all(A * Int.(x) .<= 1)
end

# Mesure le temps écoulé en secondes (temps réel, pas CPU strict).
function mesurer(f)
    t0 = time_ns()
    resultat = f()
    return resultat, (time_ns() - t0) / 1e9
end

# Construit une solution en choisissant la meilleure utilité admissible.
function greedy_algorithm(C, A)
    _, n = verifier_instance(C, A)
    solution = falses(n)
    disponibles = trues(n)

    while any(disponibles)
        utilites = fill(-Inf, n)

        for j in 1:n
            if disponibles[j] && C[j] > 0
                nb_contraintes = count(==(1), A[:, j])
                utilites[j] = nb_contraintes == 0 ? Float64(C[j]) : C[j] / nb_contraintes
            end
        end

        j = argmax(utilites)
        isfinite(utilites[j]) || break  # Aucun candidat rentable.

        solution[j] = true

        # Une colonne devient interdite si elle partage une contrainte avec j.
        lignes = A[:, j] .== 1
        if any(lignes)
            conflits = vec(any(A[lignes, :] .== 1, dims=1))
            disponibles .&= .!conflits
        end
        disponibles[j] = false  # Ne pas reprendre la colonne choisie.
    end

    @assert est_realisable(A, solution)
    return solution, score_solution(C, solution)
end

# Recherche le premier échange 1-1 qui améliore le score.
function simple_descente(C, A, solution)
    n = length(solution)
    score_actuel = score_solution(C, solution)

    for i in 1:n
        solution[i] || continue  # Retirer une colonne choisie.

        for j in 1:n
            solution[j] && continue  # Ajouter une colonne non choisie.

            candidat = copy(solution)
            candidat[i] = false
            candidat[j] = true

            if est_realisable(A, candidat)
                nouveau_score = score_solution(C, candidat)
                if nouveau_score > score_actuel
                    return candidat, nouveau_score, true
                end
            end
        end
    end

    return copy(solution), score_actuel, false
end

# Répète les échanges jusqu'à absence d'amélioration.
function recherche_locale(C, A, solution_initiale)
    solution = copy(solution_initiale)
    est_realisable(A, solution) || error("La solution initiale est irréalisable.")

    while true
        candidat, _, ameliore = simple_descente(C, A, solution)
        ameliore || break
        solution = candidat
    end

    return solution, score_solution(C, solution)
end

# Calcule l'optimum avec JuMP et GLPK.
function resolution_exacte(C, A; limite_secondes=60.0)
    m, n = verifier_instance(C, A)
    modele = Model(GLPK.Optimizer)
    set_silent(modele)
    set_time_limit_sec(modele, limite_secondes)

    @variable(modele, x[1:n], Bin)
    @objective(modele, Max, sum(C[j] * x[j] for j in 1:n))
    @constraint(modele, [i in 1:m], sum(A[i, j] * x[j] for j in 1:n) <= 1)

    (_, temps) = mesurer(() -> optimize!(modele))
    statut = termination_status(modele)
    optimum_prouve = statut == MOI.OPTIMAL
    valeur = has_values(modele) ? objective_value(modele) : missing

    return (score=valeur, temps=temps, statut=statut, optimal=optimum_prouve)
end

# Charge et résout une instance par les trois méthodes.
function resoudreSPP(fname; exact=true, limite_secondes=60.0, afficher=true)
    C, A = loadSPP(fname)
    m, n = verifier_instance(C, A)

    ((x_glouton, z_glouton), t_glouton) = mesurer(() -> greedy_algorithm(C, A))
    ((x_local, z_local), t_local) = mesurer(() -> recherche_locale(C, A, x_glouton))

    exact_resultat = exact ? resolution_exacte(C, A; limite_secondes=limite_secondes) : nothing

    if afficher
        println("\nInstance : ", fname, " (", m, " contraintes, ", n, " colonnes)")
        @printf("Glouton          : z = %.3f | temps = %.6f s\n", z_glouton, t_glouton)
        @printf("Recherche locale : z = %.3f | temps = %.6f s\n", z_local, t_local)
        println("Solution gloutonne : ", Int.(x_glouton))
        println("Solution locale    : ", Int.(x_local))
        if exact_resultat !== nothing
            println("Solveur exact : ", exact_resultat.statut)
            println("Meilleure valeur du solveur : ", exact_resultat.score)
            @printf("Temps solveur : %.6f s\n", exact_resultat.temps)
        end
    end

    return (instance=string(fname), m=m, n=n,
            glouton=(solution=x_glouton, score=z_glouton, temps=t_glouton),
            recherche_locale=(solution=x_local, score=z_local, temps=t_local),
            exact=exact_resultat)
end

# Trouve les fichiers d'instances du dossier dat (sauf fichiers cachés).
function lister_instances(dossier=joinpath(@__DIR__, "dat"))
    isdir(dossier) || error("Dossier dat introuvable : $dossier")
    fichiers = sort(filter(f -> isfile(joinpath(dossier, f)) && !startswith(f, "."), readdir(dossier)))
    return [joinpath(dossier, f) for f in fichiers]
end

# Exécute 10 instances et enregistre les résultats dans res/.
function experimentationSPP(fichiers=lister_instances(); exact=true, limite_secondes=60.0)
    length(fichiers) >= 5 || error("Il faut au moins 10 instances dans dat/ (ou fournir 10 chemins).")
    mkpath(joinpath(@__DIR__, "res"))
    sortie = joinpath(@__DIR__, "res", "resultats_EI1.csv")
    resultats = []

    open(sortie, "w") do io
        println(io, "instance;m;n;z_glouton;t_glouton_s;z_local;t_local_s;z_solveur;t_solveur_s;statut;optimal_prouve;gap_pct")
        for fichier in fichiers
            r = resoudreSPP(fichier; exact=exact, limite_secondes=limite_secondes, afficher=false)
            push!(resultats, r)
            e = r.exact
            z_exact = e === nothing ? missing : e.score
            t_exact = e === nothing ? missing : e.temps
            statut = e === nothing ? "non_calcule" : string(e.statut)
            optimal = e !== nothing && e.optimal
            # Écart uniquement si l'optimum positif est prouvé.
            gap = optimal && z_exact > 0 ? 100 * (z_exact - r.recherche_locale.score) / z_exact : missing
            valeurs = (basename(fichier), r.m, r.n, r.glouton.score, r.glouton.temps,
                       r.recherche_locale.score, r.recherche_locale.temps, z_exact, t_exact, statut, optimal, gap)
            println(io, join(string.(valeurs), ";"))
            println("Terminé : ", basename(fichier), " | z local = ", r.recherche_locale.score)
        end
    end

    println("\nRésultats enregistrés dans : ", sortie)
    return resultats
end
