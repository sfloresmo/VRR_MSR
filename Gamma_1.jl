##################################################
# Computing the value of Gamma for ORDER=1
# In that case, it is simply the bare vertex
##################################################
function compute_Gamma_1(eps1::Int64,l1::Int64,t1::Int64,
                         eps2::Int64,l2::Int64,t2::Int64,
                         eps3::Int64,l3::Int64,t3::Int64) :: Float64
    #####
    return get_Gamma_bare(eps1,l1,t1,
                          eps2,l2,t2,
                          eps3,l3,t3)
end