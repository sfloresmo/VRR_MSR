##################################################
# Definition of the dressed propagator
# We recall the convention for the spin-indices
# +  <-->  1
# -  <-->  2
# @IMPROVE -- The implementation here is particularly lazy
# It can be accelerated
# @IMPROVE -- If abs(t1-t2) > get_Nsteps_G(l),
# we could already return 0
# @IMPROVE -- If l > LMAX,
# we could already return 0
##################################################
function get_G(l::Int64,
               eps1::Int64,t1::Int64,
               eps2::Int64,t2::Int64) :: Float64
    #####
    # Coefficient (eps1,eps2) = (+,+)
    if  (eps1 == 1) &&
        (eps2 == 1)
        return C0 * get_R(l,abs(t1-t2)) # ATTENTION, not to forget to multiply by C0
    end
    #####
    # Coefficient (eps1,eps2) = (+,-)
    if  (eps1 == 1) &&
        (eps2 == 2)
        return get_R(l,t1-t2)
    end
    #####
    # Coefficient (eps1,eps2) = (-,+)
    if  (eps1 == 2) &&
        (eps2 == 1)
        return get_R(l,t2-t1)
    end
    #####
    # Coefficient (eps1,eps2) = (-,-)
    if  (eps1 == 2) &&
        (eps2 == 2)
        return 0.0
    end
    #####
end
##################################################
# Returns true if G[eps1,eps2] is necessarily equal to 0
# based simply on the values of (eps1,eps2)
# This is the case if (eps1,eps2) = (-,-) = (2,2)
##################################################
function exclusion_G_eps(eps1::Int64,
                         eps2::Int64) :: Bool
    return (eps1 == 2) && (eps2 == 2)
end
##################################################
# Returns true if G[eps1,t1,eps2,t2] is necessarily equal to 0
# based simply on the values of (eps1,t1,eps2,t2)
# @IMPROVE -- Should we also redo the test on (eps1,eps2)=(-,-)?
##################################################
function exclusion_G_eps_t(eps1::Int64,t1::Int64,
                           eps2::Int64,t2::Int64) :: Bool
    #####
    # Case (eps1,eps2) = (+,-)
    if  (eps1 == 1) &&
        (eps2 == 2) &&
        (t2 > t1)
        return true
    end
    #####
    # Case (eps1,eps2) = (-,+)
    if  (eps1 == 2) &&
        (eps2 == 1) &&
        (t1 > t2)
        return true
    end
    #####
    # Otherwise, we cannot perform any exclusions
    return false
end