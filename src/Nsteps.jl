##################################################
# Definition of the decay time of an harmonics
# as obtained from the (exact) second-order
# derivative of the correlation function at the initial time
# @ATTENTION -- This function is well-defined only for l >= 2
# @ATTENTION -- We have added a bit of depth in the sum
#               over (l1,l2) to reach correct values of T_l
##################################################
function get_Tl(l::Int64) :: Float64
    #####
    J_l = get_Jl(l)
    #####
    res = 0.0 # Initialising the result
    #####
    for l1=0:(LMAX + 1) # ATTENTION, we take some margins in the sum
        #####
        J_l1 = get_Jl(l1)
        #####
        for l2=0:(LMAX + 1) # ATTENTION, we take some margins in the sum
            #####
            J_l2 = get_Jl(l2)
            #####
            EL_ll1l2 = get_EL(l,l1,l2)
            #####
            res += (J_l1 - J_l2) * (J_l1 - J_l) * (EL_ll1l2)^(2)
            #####
        end
    end
    #####
    res *= C0 / (2 * l + 1)
    #####
    return 1 / (sqrt(res))
end
##################################################
const TAB_TL = OffsetArray(zeros(Float64,LMAX+1),0:LMAX) # Stocking the coherence time of all the harmonics
##################################################
# Filling in the array of coherence time
# @ATTENTION -- For l=0,1, we pick the maximum value of Tl
#               over all harmonics 2 <= l <= LMAX
##################################################
function TAB_TL!() :: Nothing
    #####
    max_val = - Inf # Initialising the determination of the maximum
    #####
    for l=2:LMAX
        Tl = get_Tl(l)
        TAB_TL[l] = Tl
        #####
        # We determine the maximum coherence time encountered
        # In practice, it is given by l=2
        if (max_val < Tl)
            max_val = Tl
        end
    end
    #####
    # Filling the harmonics l=0,1
    TAB_TL[0] = max_val
    TAB_TL[1] = max_val
    #####
    return nothing
end
##################################################
TAB_TL!() # Computing all the coherence times
#####
const DT = TC / TC_OVER_DT # Discretisation timestep
const NSTEPS_MAX = TC_OVER_DT * TMAX_OVER_TC # Maximum number of steps used
const T2 = TAB_TL[2] # Coherence time for the harmonics l=2, i.e. the longest coherence time
##################################################
const TAB_NSTEPS = OffsetArray(zeros(Int64,LMAX+1),0:LMAX) # Storing the number of steps for each harmonics
const TAB_T_OVER_TC = [i * DT / TC for i = 0:NSTEPS_MAX] # Table of t/Tc
##################################################
# Once we have determined the coherence time,
# we may now determine the maximum number of steps
# associated with each harmonics
# @IMPROVE -- For l=2, it would be better not to compute
#             the ratio T2 / T2, to prevent not reaching exactly 1
##################################################
function TAB_NSTEPS!() :: Nothing
    #####
    for l=0:LMAX
        #####
        if CUT # Different harmonics have different depths in time
            Tl = TAB_TL[l] # Coherence time for the current harmonics
            nsteps = ceil(Int64, NSTEPS_MAX * Tl / T2) # We set the depth of l=0,1,2 to NSTEPS_MAX
        else # All harmonics have the same depth in time
            nsteps = NSTEPS_MAX # We use the maximum number of steps
        end
        #####
        TAB_NSTEPS[l] = nsteps # Filling in the array
        #####
    end
    #####
    # Some sanity checks, for the specification of NSTEPS
    if (TAB_NSTEPS[2] != NSTEPS_MAX)
        error("TAB_NSTEPS[2] != NSTEPS_MAX")
    end
    #####
    return nothing
end
##################################################
TAB_NSTEPS!() # Computing the depth of the signals
##################################################
# Wrapped function that returns the number of steps
# for a given harmonics
# This corresponds to the NSTEPS for R[l]
##################################################
function get_Nsteps_G(l::Int64) :: Int64
    return TAB_NSTEPS[l]
end
##################################################
# Wrapped function that returns the number of steps
# for a given triangle
# In that case, we pick the max NSTEPS
# from the three harmonics involved
##################################################
function get_Nsteps_trg(l1::Int64,l2::Int64,l3::Int64) :: Int64
    return max(TAB_NSTEPS[l1], TAB_NSTEPS[l2], TAB_NSTEPS[l3])
end