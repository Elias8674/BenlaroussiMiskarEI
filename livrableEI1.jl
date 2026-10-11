# =========================================================================== #
# Compliant julia 1.x

# Using the following packages
using JuMP, GLPK
using LinearAlgebra

include("src/loadSPP.jl")
include("src/setSPP.jl")
include("src/getfname.jl")

# =========================================================================== #

function resoudreSPP(fname)
    println("\nAlgo glouton avec recherche local...")
    C, A = loadSPP(fname)
    m, n = size(A)

    function calcule_utilité(C, A, I) # vecteur de valeur, matrices des contraintes, et notre x
        nb_contrainte = count(==(1), A[:, I])
        valeur = C[I]
        
        if nb_contrainte == 0
            return valeur
        end

        return (valeur/nb_contrainte)
    end

    function greedy_algorithm(C, A)
        m, n = size(A)
        stop = false
        mask = trues(n)
        #solution = Int[]
        score_solution = 0

        solution = falses(n)


        while !stop

            utilités = Vector{Float64}(undef, n)

            for i=1:n 
                if mask[i] == false
                    utilités[i] = 0
                else
                    utilités[i] = calcule_utilité(C, A, i)
                end
            end

            # récupérer le meilleurs résultats de notre tableau d'utilités (on mets le i dans notre tableau de solution)
            # ajouter le score a noter score_solution
            meilleur_i = argmax(utilités)
        #push!(solution, meilleur_i)
            solution[meilleur_i] = true
            score_solution = score_solution + C[meilleur_i]


            # update le mask de x encore disponible
            B = A[A[:, meilleur_i] .== 1, :] # création d'une matrice ou on ne conserve que les lignes de contrainte concerné par notre meilleur_i
            
            if size(B, 1) > 0 # on vérifie que le nombre de lignes dans notre matrice est égal a 0
                mask_contrainte = vec(any(B .== 1, dims=1)) # création d'un mask qui représente nos contraintes
                mask = mask .& .!mask_contrainte # ne garde les élements de mask que ne sont pas dans mask_contrainte
            end

            mask[meilleur_i] = false


            # on verifie que tout le mask n'est pas a 0, sinon on arrête
            if sum(mask) == 0 
                stop = true
            end
        
        end

        return solution, score_solution

    end




    function check_validité(A, solution)
        comptes = A * Int.(solution)
        return all(comptes .<= 1)
    end


    function calcule_score(solution, C)
        return dot(C, solution)
    end


    function simple_descente(score_solution, solution, n, C, A)

        new_solution = copy(solution)

        for i=1:n
            if solution[i] == true
                new_solution[i] = false
                    for b=1:n
                        if solution[b] == false
                            new_solution[b] = true
                            # vérifié la validité 
                            if check_validité(A, new_solution)
                                new_score = calcule_score(new_solution, C)

                                if new_score > score_solution
                                    return [new_solution, new_score]
                                end
                            end
                            new_solution[b] = false
                        end
                    end
                    new_solution[i] = true
            end
        end
        

        return []
    end


    function recherche_local(C, A, score_solution, solution)

        stop = false
        m, n = size(A)

        score_améliorer = copy(score_solution)
        solution_améliorer = copy(solution)


        while !stop
            resultat_simple_descente = simple_descente(score_améliorer, solution_améliorer, n, C, A)
            if resultat_simple_descente == []
                stop = true
            else
                solution_améliorer = resultat_simple_descente[1]
                score_améliorer = resultat_simple_descente[2]    
            end
        end

        return solution_améliorer, score_améliorer
    end    




    solution, score = greedy_algorithm(C, A)
    solution_améliorer, score_améliorer = recherche_local(C, A, score, solution)

    println("résultats algorithme glouton")
    println("z = ", score)
    println("x = ", solution)

    println("résultats recherche local")
    println("z = ", score_améliorer)
    println("x = ", solution_améliorer)

end

# =========================================================================== #