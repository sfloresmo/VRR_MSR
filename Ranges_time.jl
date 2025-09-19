##################################################
# Helper function that returns the summation range
# wrt time for an index that is constrained by other time indices
# The arguments are the times that are already known.
# For a propagator R[l], the undetermined time must be
# within +/- get_Nsteps_G(l) of every of the already known times
# For an interaction vertex [Gamma, Lambda](l1,l2,l3),
# the undetermined time must be within
# +/- get_Nsteps_trg(l1,l2,l3) of every of the already know times
##################################################
# One constraint from a G function
##################################################
function get_range_G_t(l1::Int64,
                       t1::Int64) :: UnitRange{Int64}
    #####
    nsteps = get_Nsteps_G(l1) # Depth of the current G function
    #####
    return t1-nsteps : t1+nsteps
    #####
end
##################################################
# One constraint from a G functions
# and two constraints from a triangle
##################################################
function get_range_G_t_trg_tt(l1_a::Int64,
                              t1_a::Int64,
                              l1_b::Int64,l2_b::Int64,l3_b::Int64,
                              t1_b::Int64,t2_b::Int64) :: UnitRange{Int64}
    #####
    range_a = get_range_G_t(l1_a,t1_a)
    range_b = get_range_trg_tt(l1_b,l2_b,l3_b,t1_b,t2_b)
    #####
    return max(range_a.start,range_b.start) : min(range_a.stop,range_b.stop) # Combining the constraints
    #####
end
##################################################
# Constraints from ONE interaction triangle
# with one or two times prescribed
# @IMPROVE -- It would be better to refactor the code
##################################################
function get_range_trg_t(l1::Int64,l2::Int64,l3::Int64,
                         t1::Int64) :: UnitRange{Int64}
    #####
    nsteps = get_Nsteps_trg(l1,l2,l3) # Depth of the current triangle
    #####
    return t1-nsteps : t1+nsteps
end
#####
function get_range_trg_tt(l1::Int64,l2::Int64,l3::Int64,
                          t1::Int64,t2::Int64) :: UnitRange{Int64}
    #####
    nsteps = get_Nsteps_trg(l1,l2,l3) # Depth of the current triangle
    #####
    return max(t1,t2)-nsteps : min(t1,t2)+nsteps
end
##################################################
# Constraints from TWO interaction triangles
# with various number of times prescribed
# @IMPROVE -- It would be better to refactor the code
##################################################
function get_range_trg_t_trg_t(l1_a::Int64,l2_a::Int64,l3_a::Int64,
                               t1_a::Int64,
                               l1_b::Int64,l2_b::Int64,l3_b::Int64,
                               t1_b::Int64) :: UnitRange{Int64}
    #####
    # Depth of the triangles
    range_a = get_range_trg_t(l1_a,l2_a,l3_a,t1_a)
    range_b = get_range_trg_t(l1_b,l2_b,l3_b,t1_b)
    #####
    return max(range_a.start,range_b.start) : min(range_a.stop,range_b.stop) # Combining the constraints
    #####
end
#####
function get_range_trg_tt_trg_t(l1_a::Int64,l2_a::Int64,l3_a::Int64,
                                t1_a::Int64,t2_a::Int64,
                                l1_b::Int64,l2_b::Int64,l3_b::Int64,
                                t1_b::Int64) :: UnitRange{Int64}
    #####
    # Depth of the triangles
    range_a = get_range_trg_tt(l1_a,l2_a,l3_a,t1_a,t2_a)
    range_b = get_range_trg_t(l1_b,l2_b,l3_b,t1_b)
    #####
    return max(range_a.start,range_b.start) : min(range_a.stop,range_b.stop) # Combining the constraints
    #####
end
#####
function get_range_trg_t_trg_tt(l1_a::Int64,l2_a::Int64,l3_a::Int64,
                                t1_a::Int64,
                                l1_b::Int64,l2_b::Int64,l3_b::Int64,
                                t1_b::Int64,t2_b::Int64) :: UnitRange{Int64}
    #####
    # Using symmetry to reduce the copy of code
    return get_range_trg_tt_trg_t(l1_b,l2_b,l3_b,
                                  t1_b,t2_b,
                                  l1_a,l2_a,l3_a,
                                  t1_a) 
    #####
end
#####
function get_range_trg_tt_trg_tt(l1_a::Int64,l2_a::Int64,l3_a::Int64,
                                 t1_a::Int64,t2_a::Int64,
                                 l1_b::Int64,l2_b::Int64,l3_b::Int64,
                                 t1_b::Int64,t2_b::Int64) :: UnitRange{Int64}
    #####
    # Depth of the triangles
    range_a = get_range_trg_tt(l1_a,l2_a,l3_a,t1_a,t2_a)
    range_b = get_range_trg_tt(l1_b,l2_b,l3_b,t1_b,t2_b)
    #####
    return max(range_a.start,range_b.start) : min(range_a.stop,range_b.stop) # Combining the constraints
    #####
end