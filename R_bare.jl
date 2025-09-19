##################################################
# Definition of the bare propagator
# We recall the convention for the spin-indices
# +  <-->  1
# -  <-->  2
##################################################
const TAB_R_BARE     = OffsetArray([OffsetArray([0.0 for t=0:get_Nsteps_G(l)],0:get_Nsteps_G(l)) for l=0:LMAX],0:LMAX) # Bare response function
const TAB_INVR_BARE  = OffsetArray([OffsetArray([0.0 for t=0:get_Nsteps_G(l)],0:get_Nsteps_G(l)) for l=0:LMAX],0:LMAX) # Inverse bare response function
##################################################
# Filling in the bare response function
# R_bare = [1,1,1,...,1]
##################################################
function TAB_R_BARE!() :: Nothing
    #####
    for l=0:LMAX
        for t=0:get_Nsteps_G(l)
            TAB_R_BARE[l][t] = 1.0
        end
    end
    #####
    return nothing
end
##################################################
# Filling in the inverse of the bare response function
# R^-1 = 1/DT^2 [1,-1,0,0, ...]
##################################################
function TAB_INVR_BARE!() :: Nothing
    #####
    for l=0:LMAX
        TAB_INVR_BARE[l][0] =   1 / (DT^(2))
        TAB_INVR_BARE[l][1] = - 1 / (DT^(2))
        #####
        for t=2:get_Nsteps_G(l)
            TAB_INVR_BARE[l][t] = 0.0
        end
    end
    #####
    return nothing
end
##################################################
TAB_R_BARE!() # Pre-computing the bare response function
TAB_INVR_BARE!() # Pre-computing the inverse bare response function