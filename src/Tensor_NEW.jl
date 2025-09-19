###################################################
# Definition of the array containing the tensors Gamma, Theta and G_3
# Here, we try to be economic on memory consumption
# As a consequence, we only store the values for the non-trivial
# interaction triangles, as given by TAB_TRG.
# Then, for a given triangle, using antisymmetry,
# we store only a limited number of eps-triplets,
# as given by TAB_TRG_EPS.
# Finally, given that Gamma and Theta are stationary functions of time,
# we may always evaluate the function with
# 0 <= t1, t2, t3 <= get_Nsteps_trg(l1,l2,l3)
# and t1=0 or t2=0 or t3=0
# As a consequence, for a given (l1,l2,l3)
# and a given (eps1,eps2,eps3),
# we will need to store only 3 arrays of size (Nsteps+1)^2
# depending on whether t1=0, t2=0, or t3=0
# This allows for a (very) significant reduction of the memory imprint
##################################################
# Initialising the array with NaN
# This ensures that the code will fail
# if we have not correctly initialised some values
# With the present convention
# TAB[i_trg][i_eps][1][ta,tb]
#  -->  (t1,t2,t3) = (0,ta,tb)
# TAB[i_trg][i_eps][2][ta,tb]
#  -->  (t1,t2,t3) = (ta,0,tb)
# TAB[i_trg][i_eps][3][ta,tb]
#  -->  (t1,t2,t3) = (ta,tb,0)
const TAB_GAMMA = [[[OffsetArray([NaN for ta=0:TAB_TRG_NSTEPS[i_trg],
                                          tb=0:TAB_TRG_NSTEPS[i_trg]],
                                          0:TAB_TRG_NSTEPS[i_trg],0:TAB_TRG_NSTEPS[i_trg])
                                      for i_t = 1:3]
                                      for i_eps = 1:TAB_TRG_NB_EPS[i_trg]]
                                      for i_trg = 1:NB_TRG]
#####
const TAB_THETA = [[[OffsetArray([NaN for ta=0:TAB_TRG_NSTEPS[i_trg],
                                          tb=0:TAB_TRG_NSTEPS[i_trg]],
                                          0:TAB_TRG_NSTEPS[i_trg],0:TAB_TRG_NSTEPS[i_trg])
                                      for i_t = 1:3]
                                      for i_eps = 1:TAB_TRG_NB_EPS[i_trg]]
                                      for i_trg = 1:NB_TRG]
#####
const TAB_G_3 = [[[OffsetArray([NaN for ta=0:TAB_TRG_NSTEPS[i_trg],
                                        tb=0:TAB_TRG_NSTEPS[i_trg]],
                                          0:TAB_TRG_NSTEPS[i_trg],0:TAB_TRG_NSTEPS[i_trg])
                                    for i_t = 1:3]
                                    for i_eps = 1:TAB_TRG_NB_EPS[i_trg]]
                                    for i_trg = 1:NB_TRG]
#####
# Temporary array used to perform iterations
const TAB_GAMMA_NEW = [[[OffsetArray([NaN for ta=0:TAB_TRG_NSTEPS[i_trg],
                                              tb=0:TAB_TRG_NSTEPS[i_trg]],
                                              0:TAB_TRG_NSTEPS[i_trg],0:TAB_TRG_NSTEPS[i_trg])
                                          for i_t = 1:3]
                                          for i_eps = 1:TAB_TRG_NB_EPS[i_trg]]
                                          for i_trg = 1:NB_TRG]
##################################################
# Exclusion rule for Tensor
# All the Tensor satisfy the exclusion rule that
# Tensor[{-,1,t1},{eps2,l,t2},{eps3,l,t3}] == 0
# whatever l, (eps2,eps3), and (t1,t2,t3)
# We also recall that Tensor is fully antisymmetric
# wrt its three arguments
# Returns true if the triplets of (eps,l) has a vanishing value
##################################################
function exclusion_Tensor(eps1::Int64,l1::Int64,
                          eps2::Int64,l2::Int64,
                          eps3::Int64,l3::Int64) :: Bool
    #####
    if  ((eps1 == 2) && (l1 == 1) && (l2 == l3)) ||
        ((eps2 == 2) && (l2 == 1) && (l1 == l3)) ||
        ((eps3 == 2) && (l3 == 1) && (l1 == l2))
        #####
        return true
        #####
    end
    #####
    return false
    #####
end
##################################################
# Wrapped function that returns Gamma or Theta
# by reading from tab
# We use (i)   the full antisymmetry
#        (ii)  the stationarity in time
#        (iii) only triangular couplings matter
# @IMPROVE -- We could also implement the keyword "trustme"
#             but this part of the code is not the bottleneck
##################################################
function get_Tensor(eps1::Int64,l1::Int64,t1::Int64,
                    eps2::Int64,l2::Int64,t2::Int64,
                    eps3::Int64,l3::Int64,t3::Int64,
                    tab::Array) :: Float64
    #####
    # Checking the range in harmonics
    if  !(0 <= l1 <= LMAX) ||
        !(0 <= l2 <= LMAX) ||
        !(0 <= l3 <= LMAX)
        return 0.0
    end
    #####
    # Checking that the coupling is triangular
    if exclusion_trg(l1,l2,l3)
        return 0.0
    end
    #####
    # We use time stationarity, i.e.
    # Tensor[t1,t2,t3] = Tensor[t1+T,t2+T,t3+T] for any T
    # We shift to new times where the smallest time is 0
    t_shift = min(t1,t2,t3)
    #####
    # Determining the new times
    t1_shift = t1 - t_shift
    t2_shift = t2 - t_shift
    t3_shift = t3 - t_shift
    #####
    # Checking the range of times
    nsteps = get_Nsteps_trg(l1,l2,l3)
    #####
    if  !(0 <= t1_shift <= nsteps) ||
        !(0 <= t2_shift <= nsteps) ||
        !(0 <= t3_shift <= nsteps)
        return 0.0
    end
    #####
    # We now have a non-trivial value to compute
    # that is available in tab
    #####
    # Extracting the array with the needed information
    # to implement the symmetry
    tab_sym_loc = TAB_TRG_EPS_ALL[eps1,l1,eps2,l2,eps3,l3] # Array containing the needed information
    #####
    # Reading the information
    ind_trg   = tab_sym_loc[1]
    ind_eps   = tab_sym_loc[2]
    signature = tab_sym_loc[3]
    perm_1    = tab_sym_loc[4]
    perm_2    = tab_sym_loc[5]
    perm_3    = tab_sym_loc[6]
    #####
    # Determining the times at which to evaluate Tensor
    # ATTENTION, we must perform the permutation
    # over the shifted times
    t_old = t1_shift,
            t2_shift,
            t3_shift
    t_new = t_old[perm_1],
            t_old[perm_2],
            t_old[perm_3]
    #####
    # We now determine which of the times is equal to 0,
    # so that we know which of the arrays to access
    # And the indices of the times for which it should be read
    if     (t_new[1] == 0)
        ind_t = 1
        ta = t_new[2]
        tb = t_new[3]
    elseif (t_new[2] == 0)
        ind_t = 2
        ta = t_new[1]
        tb = t_new[3]
    else # No need to test here -- we know for sure that t_new[3] == 0
        ind_t = 3
        ta = t_new[1]
        tb = t_new[2]
    end
    #####
    # We can now access the correct array
    tab_data = tab[ind_trg][ind_eps][ind_t]
    # Read the stored value
    res_stored = tab_data[ta,tb]
    # And multiply by the signature of the permutation made
    res = signature * res_stored
    #####
    return res # Output
end
##################################################
# Function to read Gamma from TAB_GAMMA
##################################################
function get_Gamma(eps1::Int64,l1::Int64,t1::Int64,
                   eps2::Int64,l2::Int64,t2::Int64,
                   eps3::Int64,l3::Int64,t3::Int64) :: Float64
    #####
    return get_Tensor(eps1,l1,t1,
                      eps2,l2,t2,
                      eps3,l3,t3,
                      TAB_GAMMA)
    #####
end
##################################################
# Function to read Theta from TAB_THETA
##################################################
function get_Theta(eps1::Int64,l1::Int64,t1::Int64,
                   eps2::Int64,l2::Int64,t2::Int64,
                   eps3::Int64,l3::Int64,t3::Int64) :: Float64
    #####
    return get_Tensor(eps1,l1,t1,
                      eps2,l2,t2,
                      eps3,l3,t3,
                      TAB_THETA)
    #####
end
##################################################
# Function to read G_3 from TAB_G_3
##################################################
function get_G_3(eps1::Int64,l1::Int64,t1::Int64,
                 eps2::Int64,l2::Int64,t2::Int64,
                 eps3::Int64,l3::Int64,t3::Int64) :: Float64
    #####
    return get_Tensor(eps1,l1,t1,
                      eps2,l2,t2,
                      eps3,l3,t3,
                      TAB_G_3)
    #####
end
##################################################
# Filling in TAB_TENSOR
# In practice, to fully leverage parallelisation,
# we parallelise over all the 4-uple (i_trg,i_eps,ta,tb)
# @IMPROVE -- Parallelisation is hardcoded
##################################################
function TAB_TENSOR!(tab::Array, compute::F) where {F <: Function}
    # We synchronise over the outside iterations
    # i.e. we should not move on before all iterations are performed
    @sync for i_trg = 1:NB_TRG # Loop over the triangles
        #####
        Nsteps = TAB_TRG_NSTEPS[i_trg] # Number of steps for the current triangle
        nb_eps = TAB_TRG_NB_EPS[i_trg] # Number of eps-triplets to consider
        #####
        # Iterating over all the i_eps, ta, tb
        for i_eps=1:nb_eps,
            ta=0:Nsteps,
            tb=0:Nsteps
            #####
            # Making one calculation
            Threads.@spawn TAB_TENSOR_CALC!(tab,compute,
                                            i_trg,i_eps,
                                            ta,tb)
            #####
        end
    end
end
##################################################
# Function that fills in TAB_TENSOR
# for a given value of (i_trg,i_eps,ta,tb)
# Given that the call to TAB_THETA! is (by far)
# the most expensive part of the computation for ORDER=2,
# we leverage antisymmetry as much as possible
# to reduce the number of costly evaluations.
# We recall that TENSOR are fully antisymmetric,
# i.e. Tensor[a,b,c] = - Tensor[b,a,c] = - Tensor[c,b,a]
##################################################
function TAB_TENSOR_CALC!(tab::Array, compute::F,
                          i_trg::Int64,i_eps::Int64,
                          ta::Int64,tb::Int64)       where {F <: Function}
    #####
    trg = TAB_TRG[i_trg] # Current triangle
    l1, l2, l3 = trg[1], trg[2], trg[3] # Current harmonics
    #####
    tab_eps = TAB_TRG_EPS[i_trg] # Array of the eps-triplets to consider
    #####
    eps = tab_eps[i_eps] # Current eps-triplet
    eps1, eps2, eps3 = eps[1], eps[2], eps[3] # Current epsilons
    #####
    # Giving a name to the arguments p=(eps,l)
    # This is convenient to perform comparisons more easily
    p1 = eps1, l1
    p2 = eps2, l2
    p3 = eps3, l3
    #####
    # Giving a name to the three arrays to be filled
    tab_1 = tab[i_trg][i_eps][1] # For t1=0
    tab_2 = tab[i_trg][i_eps][2] # For t2=0
    tab_3 = tab[i_trg][i_eps][3] # For t3=0
    #####
    # Generically, we must fill in the values
    # v1 = ([p1,0],[p2,ta],[p3,tb])
    # v2 = ([p1,ta],[p2,0],[p3,tb])
    # v3 = ([p1,ta],[p2,tb],[p3,0])
    #####
    # If some of the p=(eps,l) are the same,
    # we can use antisymmetry to limit the number of evaluations
    #####
    # All the pairs are the same, with p1=p2=p3=p
    if (p1 == p2) && (p1 == p3)
        #####
        # Computing ([p,0],[p,ta],[p,tb])
        if (ta <= tb)
            if (ta == 0)
                val_1 = 0.0 # Obvious value
            elseif (ta == tb)
                val_1 = 0.0 # Obvious value
            else
                val_1 = compute_Tensor(eps1,l1,0,  # Corresponds to ([p,0],[p,ta],[p,tb])
                                       eps2,l2,ta,
                                       eps3,l3,tb,
                                       compute)
            end
            #####
            # Filling in the arrays
            # @IMPROVE -- if ta=tb, no need to fill in twice
            tab_1[ta,tb] =   val_1 # ([p,0],[p,ta],[p,tb])
            tab_1[tb,ta] = - val_1 # ([p,0],[p,tb],[p,ta])
            tab_2[ta,tb] = - val_1 # ([p,ta],[p,0],[p,tb])
            tab_2[tb,ta] =   val_1 # ([p,tb],[p,0],[p,ta])
            tab_3[ta,tb] =   val_1 # ([p,ta],[p,tb],[p,0])
            tab_3[tb,ta] = - val_1 # ([p,tb],[p,ta],[p,0])
            #####
        end
    #####
    # The pairs are such that p1=p2=p and p3 != p
    elseif (p1 == p2)
        #####
        # Computing ([p,ta],[p,tb],[p3,0])
        if (ta <= tb)
            if (ta == tb)
                val_3 = 0.0 # Obvious value
            else
                val_3 = compute_Tensor(eps1,l1,ta,  # Corresponds to ([p,ta],[p,tb],[p3,0])
                                       eps2,l2,tb,
                                       eps3,l3,0,
                                       compute)
            end
            #####
            # Filling in the arrays
            # @IMPROVE -- if ta=tb, no need to fill in twice
            tab_3[ta,tb] =   val_3 # ([p,ta],[p,tb],[p3,0])
            tab_3[tb,ta] = - val_3 # ([p,tb],[p,ta],[p3,0])
        end
        #####
        # Computing ([p,0],[p,ta],[p3,tb])
        if (ta == 0)
            val_1 = 0.0 # Obvious value
        else
            val_1 = compute_Tensor(eps1,l1,0,  # Corresponds to ([p,0],[p,ta],[p3,tb])
                                   eps2,l2,ta,
                                   eps3,l3,tb,
                                   compute)
        end
        #####
        # Filling in the arrays
        tab_1[ta,tb] =   val_1 # ([p,0],[p,ta],[p3,tb])
        tab_2[ta,tb] = - val_1 # ([p,ta],[p,0],[p3,tb])
    ######
    # The pairs are such that p1=p3=p and p2 != p
    elseif (p1 == p3)
        #####
        # Computing ([p,ta],[p2,0],[p,tb])
        if (ta <= tb)
            if (ta == tb)
                val_2 = 0.0 # Obvious value
            else
                val_2 = compute_Tensor(eps1,l1,ta,  # Corresponds to ([p,ta],[p2,0],[p,tb])
                                       eps2,l2,0,
                                       eps3,l3,tb,
                                       compute)
            end
            #####
            # Filling in the arrays
            # @IMPROVE -- if ta=tb, no need to fill in twice
            tab_2[ta,tb] =   val_2 # ([p,ta],[p2,0],[p,tb])
            tab_2[tb,ta] = - val_2 # ([p,tb],[p2,0],[p,ta])
        end
        #####
        # Computing ([p,0],[p2,ta],[p,tb])
        if (tb == 0)
            val_1 = 0.0 # Obvious value
        else
            val_1 = compute_Tensor(eps1,l1,0,  # Corresponds to ([p,0],[p2,ta],[p,tb])
                                   eps2,l2,ta,
                                   eps3,l3,tb,
                                   compute)
        end
        #####
        # Filling in the arrays
        tab_1[ta,tb] =   val_1 # ([p,0],[p2,ta],[p,tb])
        tab_3[tb,ta] = - val_1 # ([p,tb],[p2,ta],[p,0]) -- ATTENTION, to the order here
        #####
        # The pairs are such that p2=p3=p and p1 != p
        elseif (p2 == p3)
        #####
        # Computing ([p1,0],[p,ta],[p,tb])
        if (ta <= tb)
            if (ta == tb)
                val_1 = 0.0 # Obvious value
            else
                val_1 = compute_Tensor(eps1,l1,0,  # Corresponds to ([p1,0],[p,ta],[p,tb])
                                       eps2,l2,ta,
                                       eps3,l3,tb,
                                       compute)
            end
            #####
            # Filling in the arrays
            # @IMPROVE -- if ta=tb, no need to fill in twice
            tab_1[ta,tb] =   val_1 # ([p1,0],[p,ta],[p,tb])
            tab_1[tb,ta] = - val_1 # ([p1,0],[p,tb],[p,ta])
        end
        #####
        # Computing ([p1,ta],[p,0],[p,tb])
        if (tb == 0)
            val_2 = 0.0 # Obvious value
        else
            val_2 = compute_Tensor(eps1,l1,ta,  # Corresponds to ([p1,ta],[p,0],[p,tb])
                                   eps2,l2,0,
                                   eps3,l3,tb,
                                   compute)
        end
        #####
        # Filling in the arrays
        tab_2[ta,tb] =   val_2 # ([p1,ta],[p,0],[p,tb])
        tab_3[ta,tb] = - val_2 # ([p1,ta],[p,tb],[p,0])
    #####
    # The pairs are all different,
    # so we must compute everything
    else
        #####
        val_1 = compute_Tensor(eps1,l1,0,  # Corresponds to ([p1,0],[p2,ta],[p3,tb])
                               eps2,l2,ta,
                               eps3,l3,tb,
                               compute)
        #####
        val_2 = compute_Tensor(eps1,l1,ta,  # Corresponds to ([p1,ta],[p2,0],[p3,tb])
                               eps2,l2,0,
                               eps3,l3,tb,
                               compute)
        #####
        val_3 = compute_Tensor(eps1,l1,ta,  # Corresponds to ([p1,ta],[p2,tb],[p3,0])
                               eps2,l2,tb,
                               eps3,l3,0,
                               compute)
        #####
        # Filling in the arrays
        tab_1[ta,tb] = val_1 # ([p1,0],[p2,ta],[p3,tb])
        tab_2[ta,tb] = val_2 # ([p1,ta],[p2,0],[p3,tb])
        tab_3[ta,tb] = val_3 # ([p1,ta],[p2,tb],[p3,tb])
        #####
    end
    #####
    return nothing
end
#################################################
# Wrapper function to evaluate Tensor with the function compute
# It ensures the exact antisymmetry wrt the first two arguments
# @IMPROVE -- There might be a simpler way of symmetrising
#################################################
function compute_Tensor(eps1::Int64,l1::Int64,t1::Int64,
                        eps2::Int64,l2::Int64,t2::Int64,
                        eps3::Int64,l3::Int64,t3::Int64,
                        compute::F) where {F <: Function}
    #####
    arg1 = eps1, l1, t1 # First argument
    arg2 = eps2, l2, t2 # Second argument
    arg3 = eps3, l3, t3 # Third argument
    #####
    # Testing if two arguments are equal
    if  (arg1 == arg2) ||
        (arg1 == arg3) ||
        (arg2 == arg3)
        #####
        return 0.0
        #####
    end
    #####
    # Leveraging the exclusion rule for Tensor
    if exclusion_Tensor(eps1,l1,
                        eps2,l2,
                        eps3,l3)
        #####
        return 0.0
        #####
    end
    #####
    # We now know that all the arguments are different
    # So, we can sort them
    arg1_new, arg2_new, arg3_new, perm_sign = sort3(arg1,arg2,arg3)
    #####
    # Evaluating the function with the appropriate arguments
    return perm_sign * compute(arg1_new[1],arg1_new[2],arg1_new[3],
                               arg2_new[1],arg2_new[2],arg2_new[3],
                               arg3_new[1],arg3_new[2],arg3_new[3])
end
#################################################
# Computing the values of TAB_THETA
##################################################
function TAB_THETA!() :: Nothing
    #####
    TAB_TENSOR!(TAB_THETA, compute_Theta)
    #####
    return nothing
end
#################################################
# Computing the values of TAB_G_3
##################################################
function TAB_G_3!() :: Nothing
    #####
    TAB_TENSOR!(TAB_G_3, compute_G_3)
    #####
    return nothing
end
#################################################
# Computing the values of TAB_GAMMA_NEW
##################################################
function TAB_GAMMA_NEW!() :: Nothing
    #####
    TAB_TENSOR!(TAB_GAMMA_NEW, compute_Gamma)
    #####
    return nothing
end
##################################################
# Copying the content of TAB_GAMMA_NEW into TAB_GAMMA
##################################################
function TAB_GAMMA_COPY!() :: Nothing
    #####
    for i_trg = 1:NB_TRG # Loop over the triangles
        #####
        nb_eps = TAB_TRG_NB_EPS[i_trg] # Number of eps-triplets
        Nsteps = TAB_TRG_NSTEPS[i_trg] # Number of timesteps
        #####
        for i_eps = 1:nb_eps
            #####
            for i_t = 1:3 # Loop over index of the time set to 0
                #####
                tab     = TAB_GAMMA[i_trg][i_eps][i_t]     # Old array 
                tab_new = TAB_GAMMA_NEW[i_trg][i_eps][i_t] # New array
                #####
                # Loop over (ta,tb) in the order appropriate for the array
                for tb=0:Nsteps,
                    ta=0:Nsteps
                    #####
                    tab[ta,tb] = tab_new[ta,tb]
                    #####
                end
            end
        end
    end
    #####
    return nothing
end
##################################################
# Picking the function to compute Gamma
##################################################
if     (ORDER == 0)
    include("Gamma_0.jl")
    const compute_Gamma = compute_Gamma_0
elseif (ORDER == 1)
    include("Gamma_1.jl")
    const compute_Gamma = compute_Gamma_1
elseif (ORDER == 2)
    include("Gamma_2.jl")
    const compute_Gamma = compute_Gamma_2
end