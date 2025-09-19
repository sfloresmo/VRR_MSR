include("W3.jl") # To compute the diagram W3
##################################################
# Computing the values of Theta_123
# as given by the W3 diagram
# ATTENTION, here @noinline is essential,
# otherwise the compiler tries to optimise the bounds
# in the loops and can get stuck when optimising/compiling the code
# ATTENTION, we use the keyword "trustme" when evaluating Lambda,
# to turn off the sanity checks and accelerate the evaluation
##################################################
function compute_Theta(eps1::Int64,l1::Int64,t1::Int64,
                       eps2::Int64,l2::Int64,t2::Int64,
                       eps3::Int64,l3::Int64,t3::Int64) :: Float64
    ##############
    s = 0.0 # Initialising the result
        #####
    for l4=0:LMAX
    for l5=0:LMAX
        #####
        @noinline exclusion_trg(l1,l5,l4) && continue # Triangular exclusion for Lambda_154
        #####
    for l6=0:LMAX
        #####
        @noinline exclusion_trg(l2,l4,l6) && continue # Triangular exclusion for Lambda_246
        @noinline exclusion_trg(l3,l6,l5) && continue # Triangular exclusion for Lambda_365
        #####
        W3 = get_W3(l1,l2,l3,l4,l5,l6)
        (W3 == 0.0) && continue
        #####
    for eps4 =1:2
        #####
        @noinline exclusion_Tensor2(eps2,l2,eps4,l4,l6) && continue # Exclusion rule for Lambda
        #####
    for eps5 =1:2
        #####
        @noinline exclusion_Tensor2(eps1,l1,eps5,l5,l4) && continue # Exclusion rule for Lambda
        #####
    for eps6 =1:2
        #####
        @noinline exclusion_Tensor2(eps3,l3,eps6,l6,l5) && continue # Exclusion rule for Lambda
        #####
    for t4 =@noinline get_range_trg_t_trg_t(l1,l4,l5,t1,l2,l4,l6,t2)
    for t5 =@noinline get_range_trg_tt_trg_t(l1,l4,l5,t1,t4,l3,l5,l6,t3)
        #####
        Lambda_154 = get_Lambda(eps1,l1,t1,
                                eps5,l5,t5,
                                eps4,l4,t4;
                                trustme=true) # We turn off the sanity checks
        (Lambda_154 == 0.0) && continue
        #####
    for t6 =@noinline get_range_trg_tt_trg_tt(l2,l4,l6,t2,t4,l3,l5,l6,t3,t5)
        #####
        Lambda_246 = get_Lambda(eps2,l2,t2,
                                eps4,l4,t4,
                                eps6,l6,t6;
                                trustme=true) # We turn off the sanity checks
        (Lambda_246 == 0.0) && continue
        #####
        Lambda_365 = get_Lambda(eps3,l3,t3,
                                eps6,l6,t6,
                                eps5,l5,t5;
                                trustme=true) # We turn off the sanity checks
        (Lambda_365 == 0.0) && continue
        #####
        s += W3 * Lambda_154 * Lambda_246 * Lambda_365
        #####
    end
    end
    end
    end
    end
    end
    end
    end
    end
    #####
    # We have performed 3 sums over time,
    # hence we must multiply the final result by DT^3
    # to comply with Riemann sum formula
    s *= DT^(3)
    #####
    return s # Output
end