##################################################
# Computing the value of Gamma for ORDER=2
# Attention, we suppose that Theta is filled
##################################################
function compute_Gamma_2(eps1::Int64,l1::Int64,t1::Int64,
                         eps2::Int64,l2::Int64,t2::Int64,
                         eps3::Int64,l3::Int64,t3::Int64) :: Float64
    #####
    Gamma_bare_123 = get_Gamma_bare(eps1,l1,t1,
                                    eps2,l2,t2,
                                    eps3,l3,t3)
    Theta_123 = get_Theta(eps1,l1,t1,
                          eps2,l2,t2,
                          eps3,l3,t3)
    #####
    return Gamma_bare_123 + Theta_123
end
