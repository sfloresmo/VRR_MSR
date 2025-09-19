##################################################
# Computing the values of G_123 defined as
# G_123 = \sum_{eps1p,t1p;eps2p,t2p;eps3p,t3p} G_11p G_22p G_33p Gamma_1p2p3p
# ATTENTION, there is no sum over l
# ATTENTION, here @noinline is essential
##################################################
# Computing G_123
##################################################
function compute_G_3(eps1::Int64,l1::Int64,t1::Int64,
                     eps2::Int64,l2::Int64,t2::Int64,
                     eps3::Int64,l3::Int64,t3::Int64) :: Float64
    ##############
    s = 0.0 # Initialising the result
    #####
    for eps1p = 1:2
        #####
        @noinline exclusion_G_eps(eps1,eps1p) && continue # Exclusion on eps for G_11p
        #####
    for eps2p = 1:2
        #####
        @noinline exclusion_G_eps(eps2,eps2p) && continue # Exclusion on eps for G_22p
        #####
    for eps3p = 1:2
        #####
        @noinline exclusion_G_eps(eps3,eps3p) && continue # Exclusion on eps for G_33p
        #####
    for t1p = @noinline get_range_G_t(l1,t1)
        #####
        @noinline exclusion_G_eps_t(eps1,t1,eps1p,t1p) && continue # Exclusion on eps,t for G_11p
        #####
        G_11p = get_G(l1,eps1,t1,eps1p,t1p)
        (G_11p == 0.0) && continue
        #####
    for t2p = @noinline get_range_G_t(l2,t2)
        #####
        @noinline exclusion_G_eps_t(eps2,t2,eps2p,t2p) && continue # Exclusion on eps,t for G_22p
        #####
        G_22p = get_G(l2,eps2,t2,eps2p,t2p)
        (G_22p == 0.0) && continue
        #####
    for t3p = @noinline get_range_G_t(l3,t3)
        #####
        @noinline exclusion_G_eps_t(eps3,t3,eps3p,t3p) && continue # Exclusion on eps,t for G_33p
        #####
        G_33p = get_G(l3,eps3,t3,eps3p,t3p)
        (G_33p == 0.0) && continue
        #####
        Gamma_1p2p3p = get_Gamma(eps1p,l1,t1p,
                                 eps2p,l2,t2p,
                                 eps3p,l3,t3p)
        (Gamma_1p2p3p == 0.0) && continue
        #####
        s += G_11p * G_22p * G_33p * Gamma_1p2p3p
        #####
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
    s *= DT^3
    #####
    return s # Output
end
