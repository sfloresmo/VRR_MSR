##################################################
# Definition of the response function
# In practice, we store single array containing
# R[l][t] for
# 0 <= l <= LMAX
# 0 <= t <= NSTEPS(l) * DT
##################################################
const TAB_R     = OffsetArray([OffsetArray([0.0 for t=0:get_Nsteps_G(l)],0:get_Nsteps_G(l)) for l=0:LMAX],0:LMAX) # Response function
const TAB_R_NEW = OffsetArray([OffsetArray([0.0 for t=0:get_Nsteps_G(l)],0:get_Nsteps_G(l)) for l=0:LMAX],0:LMAX) # Updated response function
const TAB_INVR  = OffsetArray([OffsetArray([0.0 for t=0:get_Nsteps_G(l)],0:get_Nsteps_G(l)) for l=0:LMAX],0:LMAX) # Inverse of the response function
##################################################
# Function returning the value of R
##################################################
function get_R(l::Int64,t::Int64) :: Float64
    #####
    if !(0 <= l <= LMAX) # We are not within the considered harmonics range
        return 0.0
    end
    #####
    # Imposing causality
    # ATTENTION, we suppose that R[t=0] = 1
    if (t < 0)
        return 0.0
    end
    #####
    if (t > get_Nsteps_G(l)) # We are not within the considered time range for this harmonics
        return 0.0
    end
    #####
    # We are then safe to return the value from the array
    return TAB_R[l][t]
end
##################################################
# Updates TAB_INVR for a given value of Sigma
# This is Dyson equation
# ATTENTION, we assume that Sigma is up-to-date
##################################################
function TAB_INVR!() :: Nothing
    for l=0:LMAX
        for t=0:get_Nsteps_G(l)
            TAB_INVR[l][t] = TAB_INVR_BARE[l][t] - TAB_SIGMA[l][t]
        end
    end
    #####
    return nothing
end
##################################################
# Updates the value of TAB_R_NEW by
# 1/ solving the Dyson equation for R^-1
# 2/ inverting R^-1
##################################################
function TAB_R_NEW!() :: Nothing
    #####
    TAB_INVR!() # Computing R^-1
    INV_TAB_INVR!() # Inverting R^-1
    #####
    return nothing
end
##################################################
# Copying the content of TAB_R_NEW into TAB_R
##################################################
function TAB_R_COPY!() :: Nothing
    for l=0:LMAX
        for t=0:get_Nsteps_G(l)
            #####
            TAB_R[l][t] = TAB_R_NEW[l][t]
            #####
        end
    end
    #####
    return nothing
end