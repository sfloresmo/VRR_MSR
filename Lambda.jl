##################################################
# Computing the values of Lambda_bare_123 and Lambda_123 defined as
# Lambda_bare_123 = \sum_{eps3p,t3p} Gamma_bare_123p G_3p3
# Lambda_123 = \sum_{eps3p,t3p} Gamma_123p G_3p3
# [Here, G_3p3 is fully symmetric]
# ATTENTION, there is no sum over l3p
# ATTENTION, here @noinline is essential
##################################################
# Computing Lambda_123 and Lambda_bare_123
# We write a generic function to prevent code copy-paste
##################################################
function compute_Lambda_generic(eps1::Int64,l1::Int64,t1::Int64,
                                eps2::Int64,l2::Int64,t2::Int64,
                                eps3::Int64,l3::Int64,t3::Int64,
                                get_Gamma_generic::F) where {F <: Function}
    ##############
    s = 0.0 # Initialising the result
    #####
    for eps3p = 1:2
        #####
        @noinline exclusion_G_eps(eps3p,eps3) && continue # Exclusion on eps for G_3p3
        #####
    for t3p = @noinline get_range_G_t_trg_tt(l3,t3,l1,l2,l3,t1,t2)
        #####
        @noinline exclusion_G_eps_t(eps3p,t3p,eps3,t3) && continue # Exclusion on eps,t for G_3p3
        #####
        G_3p3 = get_G(l3,eps3p,t3p,eps3,t3)
        (G_3p3 == 0.0) && continue
        #####
        Gamma_generic_123p = get_Gamma_generic(eps1 ,l1,t1,
                                               eps2 ,l2,t2,
                                               eps3p,l3,t3p) # We have l3p=l3, through G_3p3
        (Gamma_generic_123p == 0.0) && continue
        #####
        s += Gamma_generic_123p * G_3p3
        #####
    end
    end
    #####
    # We have performed 1 sum over time,
    # hence we must multiply the final result by DT
    # to comply with Riemann sum formula
    s *= DT
    #####
    return s # Output
end
##################################################
# To compute Lambda_bare
##################################################
function compute_Lambda_bare(eps1::Int64,l1::Int64,t1::Int64,
                             eps2::Int64,l2::Int64,t2::Int64,
                             eps3::Int64,l3::Int64,t3::Int64)
    #####
    compute_Lambda_generic(eps1,l1,t1,
                           eps2,l2,t2,
                           eps3,l3,t3,
                           get_Gamma_bare)
    #####
end
##################################################
# To compute Lambda
##################################################
function compute_Lambda(eps1::Int64,l1::Int64,t1::Int64,
                        eps2::Int64,l2::Int64,t2::Int64,
                        eps3::Int64,l3::Int64,t3::Int64)
    #####
    compute_Lambda_generic(eps1,l1,t1,
                           eps2,l2,t2,
                           eps3,l3,t3,
                           get_Gamma)
    #####
end