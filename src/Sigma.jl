##################################################
# Array containing the first column of the matrix
# Sigma[(-,l,t),(+,l,0)] with
# + 0 <= l <= LMAX
# + 0 <= t <= NSTEPS
##################################################
const TAB_SIGMA = OffsetArray([OffsetArray([0.0 for t=0:get_Nsteps_G(l)],0:get_Nsteps_G(l)) for l=0:LMAX],0:LMAX) # Self-energy
##################################################
# Computation of the self-energy
# ATTENTION, we make sure that Sigma[l,t=0] = 0
# so that the response function always starts at 1 for t=0
# [Ideally, this should have been satisfied exactly.]
# @IMPROVE -- Parallelisation to be improved
##################################################
function TAB_SIGMA!() :: Nothing
    #####
    # Imposing the value for t=0
    # ATTENTION, VERY IMPORTANT
    # I checked that removing this condition
    # does not affect the prediction [except in a minor fashion for the starting point]
    for l=0:LMAX
        TAB_SIGMA[l][0] = 0.0
    end
    #####
    @parallel for l=0:LMAX
        for t=1:get_Nsteps_G(l) # ATTENTION, we start at t=1
            #####
            TAB_SIGMA[l][t] = compute_Sigma(l,t)
            #####
        end
    end
    #####
    return nothing
end
#################################################
# Computing Sigma
# ATTENTION, here @noinline is essential
##################################################
function compute_Sigma(l::Int64,t::Int64) :: Float64
    #####
    # Fixing the outside indices
    eps1 = 2 # eps_1 = -
    l1 = l
    t1 = t
    #####
    eps1p = 1 # eps_1p = +
    l1p = l
    t1p = 0
    #####
    s = 0.0 # Initialising the result
    #####
    for l2=0:LMAX
    for l3=0:LMAX
        #####
        # @IMPROVE -- Actually, these two tests are identical, since l1p=l1
        #             Though, to avoid confusion, it is better to leave as is
        @noinline exclusion_trg(l1, l3,l2) && continue # Triangular exclusion for Lambda_bare_132
        @noinline exclusion_trg(l1p,l2,l3) && continue # Triangular exclusion for Lambda_1p23
        #####
        W2 = get_W2(l1,l2,l3) # Value of the contraction diagram
        #####
    for eps2 =1:2
    for eps3 =1:2
        #####
    for t2= @noinline get_range_trg_t_trg_t(l1,l2,l3,t1,l1p,l2,l3,t1p)
    for t3= @noinline get_range_trg_tt_trg_tt(l1,l2,l3,t1,t2,l1p,l2,l3,t1p,t2) # Constraint from (t1,t2,t1p)
        #####
        Lambda_bare_132 = get_Lambda_bare(eps1,l1,t1,
                                          eps3,l3,t3,
                                          eps2,l2,t2)
        (Lambda_bare_132 == 0.0) && continue
        #####
        Lambda_1p23 = get_Lambda(eps1p,l1p,t1p,
                                 eps2,l2,t2,
                                 eps3,l3,t3)
        (Lambda_1p23 == 0.0) && continue
        #####
        s += W2 * Lambda_bare_132 * Lambda_1p23
        #####
    end
    end
    end
    end
    end
    end
    #####
    # We have performed 2 sums over time,
    # hence we must multiply the final result by DT^2
    # to comply with Riemann sum formula
    s *= DT^(2)
    #####
    # Adding prefactor (-1/2)
    return - 0.5 * s # Output
end

