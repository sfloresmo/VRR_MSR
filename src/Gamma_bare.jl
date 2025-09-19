##################################################
# Definition of the array containing the bare vertex
# Here, we try to be economic on memory consumption
# so that we only store the values for the non-trivial
# interaction triangles, as given by TAB_TRG
# In practice, the array TAB_GAMMA_BARE therefore
# onlyse Gamma_bare[(eps1,l1,t1),(eps2,l2,t2),(eps3,l3,t3)] for
# + (l1,l2,l3) in TAB_TRG
# + (eps1,eps2,eps3) in TAB_TRG_EPS of (l1,l2,l3)
# + t1 = t2 = t3
##################################################
# Initialising the array with NaN
# This ensures that the code will fail
# if we have not initialised some of the values
const TAB_GAMMA_BARE = [[NaN
                         for i_eps = 1:TAB_TRG_NB_EPS[i_trg]]
                         for i_trg = 1:NB_TRG]
##################################################
# Simple exclusion rule that returns true
# if (eps1,eps2,eps3) is not isomorphic to (+,+,-)
# In that case, we know for sure that Gamma_bare = 0.0
##################################################
function exclusion_Gamma_bare_eps(eps1::Int64,
                                  eps2::Int64,
                                  eps3::Int64) :: Bool
    #####
    # Case (+,+,-)
    if  (eps1 == 1) &&
        (eps2 == 1) &&
        (eps3 == 2)
        return false
    end
    #####
    # Case (+,-,+)
    if  (eps1 == 1) &&
        (eps2 == 2) &&
        (eps3 == 1)
        return false
    end
    #####
    # Case (-,+,+)
    if  (eps1 == 2) &&
        (eps2 == 1) &&
        (eps3 == 1)
        return false
    end
    #####
    # Otherwise, the result is necessarily zero
    return true
end
##################################################
# Analytical expression of
# Gamma_bare[(+,l1,t),(+,l2,t),(-,l3,t)]
##################################################
function compute_Gamma_bare_ppm(l1::Int64,
                                l2::Int64,
                                l3::Int64) :: Float64
    #####
    # Checking that it is a triangle
    if exclusion_trg(l1,l2,l3)
        return 0.0
    end
    #####
    # Determining the coupling coefficient
    J_1, J_2 = get_Jl(l1), get_Jl(l2)
    Delta_Jl = J_2 - J_1
    #####
    EL_123 = get_EL(l1,l2,l3) # Elsasser coefficient
    #####
    # Value of the coefficient
    # ATTENTION, to the prefactor in 1/DT^2
    # coming from the two Dirac deltas
    res = (1/DT^2) * EL_123 * Delta_Jl
    #####
    return res # Output
end
##################################################
# Analytical expression of
# Gamma_bare[(eps1,l1,t),(eps2,l2,t),(eps3,l3,t)]
# Here, we use antisymmetry, so that our goal is to
# fall back on the evaluation from compute_Gamma_bare_ppm
##################################################
function compute_Gamma_bare(eps1::Int64,l1::Int64,
                            eps2::Int64,l2::Int64,
                            eps3::Int64,l3::Int64) :: Float64
    #####
    # Checking that it is a triangle
    if exclusion_trg(l1,l2,l3)
        return 0.0
    end
    #####
    # Checking that (eps1,eps2,eps3) ~ (+,+,-)
    if exclusion_Gamma_bare_eps(eps1,eps2,eps3)
        return 0.0
    end
    #####
    # Using permutations to transform to the case
    # (eps1,eps2,eps3) = (+,+,-)
    #####
    # (+,+,-)  -->  No permutations
    if  (eps1 == 1) &&
        (eps2 == 1) &&
        (eps3 == 2)
        return compute_Gamma_bare_ppm(l1,l2,l3)
    end
    #####
    # (+,-,+)  --> (-1) * (+,+,-)  with  2 <-> 3
    if  (eps1 == 1) &&
        (eps2 == 2) &&
        (eps3 == 1)
        return -1 * compute_Gamma_bare_ppm(l1,l3,l2)
    end
    #####
    # (-,+,+)  -->  (-1) * (+,+,-)  with  1 <-> 3
    if  (eps1 == 2) &&
        (eps2 == 1) &&
        (eps3 == 1)
        return -1 * compute_Gamma_bare_ppm(l3,l2,l1)
    end
    #####
    # Otherwise, we return a bad value
    # because we reached a case we should not have
    # This is used as a sanity check
    return NaN
end
##################################################
# Computing the values of Gamma_bare
# for the considered triangles and values of eps
##################################################
function TAB_GAMMA_BARE!() :: Nothing
    #####
    for i_trg = 1:NB_TRG # Loop over the triangles
        #####
        trg = TAB_TRG[i_trg] # Current triangle
        l1, l2, l3 = trg[1], trg[2], trg[3] # Current harmonics
        #####
        tab_eps = TAB_TRG_EPS[i_trg] # Array of the eps-triplets to consider
        nb_eps = TAB_TRG_NB_EPS[i_trg] # Number of eps-triplets to consider
        #####
        for i_eps = 1:nb_eps # Loop over the eps-triplets
            #####
            eps = tab_eps[i_eps] # Current eps-triplet
            eps1, eps2, eps3 = eps[1], eps[2], eps[3] # Current epsilons
            #####
            # We bother computing this value only if
            # (eps1,eps2,eps3) ~ (+,+,-)
            if !exclusion_Gamma_bare_eps(eps1,eps2,eps3)
                #####
                # Computing the value of Gamma
                val = compute_Gamma_bare(eps1,l1,
                                         eps2,l2,
                                         eps3,l3)
                #####
            else # Otherwise, Gamma_bare = 0.0
                #####
                val = 0.0
                #####
            end
            #####
            # Filling in the array
            TAB_GAMMA_BARE[i_trg][i_eps] = val
            #####
        end
    end
    #####
    return nothing
end
##################################################
TAB_GAMMA_BARE!() # Pre-computing the value of TAB_GAMMA_BARE
##################################################
# Wrapped function that returns Gamma_bare
# by reading from TAB_GAMMA_BARE
# We use (i)   the full antisymmetry
#        (ii)  only triangular couplings matter
#        (iii) the fact that it is non-zero only for
#              (eps1,eps2,eps3) ~ (+,+,-)
##################################################
function get_Gamma_bare(eps1::Int64,l1::Int64,t1::Int64,
                        eps2::Int64,l2::Int64,t2::Int64,
                        eps3::Int64,l3::Int64,t3::Int64) :: Float64
    #####
    # Checking that all the times are equal
    if  !(t1 == t2 == t3)
        return 0.0
    end
    #####
    # Checking the range in harmonics
    if !(0 <= l1 <= LMAX) ||
       !(0 <= l2 <= LMAX) ||
       !(0 <= l3 <= LMAX) 
       return 0.0
    end
    #####
    # Checking that the harmonics is triangular
    if exclusion_trg(l1,l2,l3)
        return 0.0
    end
    #####
    # Checking that the eps-triplet is isomorphic to (+,+,-)
    if exclusion_Gamma_bare_eps(eps1,eps2,eps3)
        return 0.0
    end
    #####
    # We now have a non-trivial value to compute
    # that is available in TAB_GAMMA_BARE
    #####
    # Extracting the array with the needed information
    # to implement the symmetry
    tab_sym_loc = TAB_TRG_EPS_ALL[eps1,l1,eps2,l2,eps3,l3] # Array containing the needed information
    #####
    # Reading the information
    ind_trg   = tab_sym_loc[1]
    ind_eps   = tab_sym_loc[2]
    signature = tab_sym_loc[3]
    #####
    # Value pre-stored (without the signature of the permutation)
    res_stored = TAB_GAMMA_BARE[ind_trg][ind_eps]
    #####
    # Final result
    res = signature * res_stored
    #####
    return res # Output
end