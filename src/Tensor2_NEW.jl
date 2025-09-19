##################################################
# Definition of the array containing the tensor Lambda_bare and Lambda
# Here, we try to be economic on memory consumption
# As a consequence, we only store values for non-trivial
# interaction triangles2, as given by TAB_TRG2
# Then, for a given triangle2, using antisymmetry,
# we store only a limited number of eps-triplets,
# as given by TAB_TRG2_EPS
# Finally, given that Lambda_bare and Lambda are stationary functions of time,
# we may always evaluate the function with
# 0 <= t1, t2, t3 <= get_Nsteps_trg(l1,l2,l3)
# and t1=0 or t2=0 or t3=0
# As a consequece, for a given (l1,l2,l3)
# and a given (eps1,eps2,eps3),
# we will need to store only 3 arrays of size (Nsteps+1)^3,
# depending on whether t1=0, t2=0, t3=0
# This allows for a (very) significant reduction of the memory imprint
##################################################
# Initialising the array with NaN
# This ensures that the code will fail
# if we have not correctly initialised some values
# With the present convention
# TAB[i_trg2][i_eps][1][ta,tb]
#  -->  (t1,t2,t3) = (0,ta,tb)
# TAB[i_trg2][i_eps][2][ta,tb]
#  -->  (t1,t2,t3) = (ta,0,tb)
# TAB[i_trg2][i_eps][3][ta,tb]
#  -->  (t1,t2,t3) = (ta,tb,0)
const TAB_LAMBDA_BARE = [[[OffsetArray([NaN for ta=0:TAB_TRG2_NSTEPS[i_trg2],
                                                tb=0:TAB_TRG2_NSTEPS[i_trg2]],
                                                0:TAB_TRG2_NSTEPS[i_trg2],0:TAB_TRG2_NSTEPS[i_trg2])
                                            for i_t = 1:3]
                                            for i_eps = 1:TAB_TRG2_NB_EPS[i_trg2]]
                                            for i_trg2 = 1:NB_TRG2]
##################################################
const TAB_LAMBDA = [[[OffsetArray([NaN for ta=0:TAB_TRG2_NSTEPS[i_trg2],
                                           tb=0:TAB_TRG2_NSTEPS[i_trg2]],
                                           0:TAB_TRG2_NSTEPS[i_trg2],0:TAB_TRG2_NSTEPS[i_trg2])
                                       for i_t = 1:3]
                                       for i_eps = 1:TAB_TRG2_NB_EPS[i_trg2]]
                                       for i_trg2 = 1:NB_TRG2]
##################################################
# Exclusion rule for Tensor2 
# All the Tensor2 satisfy the exclusion rule that
# Tensor2[{-,1,t1},{eps2,l,t2},{eps3,l,t3}] == 0
# whatever l, (eps2,eps3), and (t1,t2,t3)
# We also recall that Tensor2 is fully antisymmetric
# wrt its two first arguments
# Returns true if the triplets of (eps,l) has a vanishing value
# ATTENTION, this function is (anti)symmetric only wrt its two first indices
##################################################
function exclusion_Tensor2(eps1::Int64,l1::Int64,
                           eps2::Int64,l2::Int64,
                           l3::Int64) :: Bool
    #####
    if  ((eps1 == 2) && (l1 == 1) && (l2 == l3)) ||
        ((eps2 == 2) && (l2 == 1) && (l1 == l3))
        #####
        return true
        #####
    end
    #####
    return false
    #####
end
##################################################
# Wrapped function that returns Lambda_bare or Lambda
# by reading from tab
# We use (i)   the restricted antisymmetry
#        (ii)  the stationarity in time
#        (iii) only triangular couplings matter
# We add a keyword argument, "trustme",
# to allow for some bounds checks to be removed
# The keyword trustme=true is only turned on in compute_Theta
##################################################
function get_Tensor2(eps1::Int64,l1::Int64,t1::Int64,
                     eps2::Int64,l2::Int64,t2::Int64,
                     eps3::Int64,l3::Int64,t3::Int64,
                     tab::Array;
                     trustme=nothing) :: Float64
    #####
    # Checking the range in harmonics
    if trustme == nothing # We perform the bound check on lmax
        if  !(0 <= l1 <= LMAX) ||
            !(0 <= l2 <= LMAX) ||
            !(0 <= l3 <= LMAX)
            return 0.0
        end
    end
    #####
    # Checking that the coupling is triangular
    if trustme == nothing # We perform the bound check on (l1,l2,l3)
        if exclusion_trg(l1,l2,l3)
            return 0.0
        end
    end
    #####
    # We use time stationarity, i.e.
    # Tensor2[t1,t2,t3] = Tensor2[t1+T,t2+T,t3+T] for any T
    # We shift to new times where the smallest time is 0
    t_shift = min(t1,t2,t3)
    #####
    # Determining the new times
    t1_shift = t1 - t_shift
    t2_shift = t2 - t_shift
    t3_shift = t3 - t_shift
    #####
    # Checking the range of times
    if trustme == nothing
        Nsteps = get_Nsteps_trg(l1,l2,l3)
        if  !(0 <= t1_shift <= Nsteps) ||
            !(0 <= t2_shift <= Nsteps) ||
            !(0 <= t3_shift <= Nsteps)
            return 0.0
        end
    end
    #####
    # We now have a non-trivial value to compute
    # that is available in tab
    #####
    # Extracting the array with the needed information
    # to implement the symmetry
    tab_sym_loc = TAB_TRG2_EPS_ALL[eps1,l1,eps2,l2,eps3,l3] # Array containing the needed information
    #####
    # Reading the information
    ind_trg2  = tab_sym_loc[1]
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
    tab_data = tab[ind_trg2][ind_eps][ind_t]
    # Read the stored value
    res_stored = tab_data[ta,tb]
    # And multiply by the signature of the permutation made
    res = signature * res_stored
    #####
    return res # Output
end
####################################################
# Function to read Lambda_bare from TAB_LAMBDA_BARE
####################################################
function get_Lambda_bare(eps1::Int64,l1::Int64,t1::Int64,
                         eps2::Int64,l2::Int64,t2::Int64,
                         eps3::Int64,l3::Int64,t3::Int64) :: Float64
    #####
    return get_Tensor2(eps1,l1,t1,
                       eps2,l2,t2,
                       eps3,l3,t3,
                       TAB_LAMBDA_BARE)
    #####
end
##################################################
# Function to read Lambda from TAB_LAMBDA
# We add a keyword argument, "trustme",
# to turn off some of the sanity checks
# This keyword is only used with compute_Theta
##################################################
function get_Lambda(eps1::Int64,l1::Int64,t1::Int64,
                    eps2::Int64,l2::Int64,t2::Int64,
                    eps3::Int64,l3::Int64,t3::Int64;
                    trustme=nothing) :: Float64
    #####
    return get_Tensor2(eps1,l1,t1,
                       eps2,l2,t2,
                       eps3,l3,t3,
                       TAB_LAMBDA;
                       trustme=trustme)
    #####
end
##################################################
# Filling in TAB_TENSOR2
# In practice, to fully leverage parallelisation,
# we parallelise over the 4-uple (i_trg2,i_eps,ta,tb)
##################################################
function TAB_TENSOR2!(tab::Array, compute::F) where {F <: Function}
    #####
    # We synchronise over the outside iterations
    # i.e. we should not move on before all iterations are performed
    @sync for i_trg2 = 1:NB_TRG2
        #####
        Nsteps = TAB_TRG2_NSTEPS[i_trg2] # Number of steps
        nb_eps = TAB_TRG2_NB_EPS[i_trg2] # Number of eps-triplets to consider
        #####
        # Iterating over all the i_eps, ta, tb
        for i_eps=1:nb_eps,
            ta=0:Nsteps,
            tb=0:Nsteps
            #####
            # Making one calculation
            Threads.@spawn TAB_TENSOR2_CALC!(tab,compute,
                                             i_trg2,i_eps,
                                             ta,tb)
            #####
        end
    end
    #####
    return nothing
end
##################################################
# Function that fills in TAB_TENSOR2
# for a given value of (i_trg2,i_eps,ta,tb)
# To reduce a bit the number of evaluations,
# we leverage antisymmetry
# We recall that TENSOR2 are fully antisymmetric
# wrt their two first indices, i.e.
# Tensor2[a,b,c] = - Tensor2[b,a,c]
##################################################
function TAB_TENSOR2_CALC!(tab::Array, compute::F,
                           i_trg2::Int64,i_eps,
                           ta::Int64,tb::Int64) where {F <: Function}
    #####
    trg2 = TAB_TRG2[i_trg2] # Current triangle2
    l1, l2, l3 = trg2[1], trg2[2], trg2[3] # Current harmonics
    #####
    tab_eps = TAB_TRG2_EPS[i_trg2] # Array of the eps-triplets to consider
    #####
    eps = tab_eps[i_eps] # Current eps-triplet
    eps1, eps2, eps3 = eps[1], eps[2], eps[3] # Current epsilons
    #####
    # Giving a name to the arguments p=(eps,l)
    p1 = eps1, l1
    p2 = eps2, l2
    ####
    # Giving a name to the three arrays to be filled
    tab_1 = tab[i_trg2][i_eps][1] # For t1=0
    tab_2 = tab[i_trg2][i_eps][2] # For t2=0
    tab_3 = tab[i_trg2][i_eps][3] # For t3=0
    #####
    # Generically, we must fill in the values
    # v1 = ([p1,0],[p2,ta],[p3,tb])
    # v2 = ([p1,ta],[p2,0],[p3,tb])
    # v3 = ([p1,ta],[p2,tb],[p3,0])
    #####
    # If the two first p=(eps,l) are the same,
    # we can use antisymmetry to reduce the number of evaluations
    #####
    # The two first pairs are the same
    if (p1 == p2)
        #####
        # Computing ([p,ta],[p,tb],[p3,0])
        if (ta <= tb)
            if (ta == tb)
                val_3 = 0.0 # Obvious value
            else
                val_3 = compute_Tensor2(eps1,l1,ta,  # Corresponds to ([p,ta],[p,tb],[p3,0])
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
            val_1 = compute_Tensor2(eps1,l1,0,  # Corresponds to ([p,0],[p,ta],[p3,tb])
                                    eps2,l2,ta,
                                    eps3,l3,tb,
                                    compute)
        end
        #####
        # Filling in the arrays
        tab_1[ta,tb] =   val_1 # ([p,0],[p,ta],[p3,tb])
        tab_2[ta,tb] = - val_1 # ([p,ta],[p,0],[p3,tb])                
    #####
    # The two first pairs are different
    # so we must compute everything
    else
        #####
        val_1 = compute_Tensor2(eps1,l1,0,  # Corresponds to ([p1,0],[p2,ta],[p3,tb])
                                eps2,l2,ta,
                                eps3,l3,tb,
                                compute)
        #####
        val_2 = compute_Tensor2(eps1,l1,ta,  # Corresponds to ([p1,ta],[p2,0],[p3,tb])
                                eps2,l2,0,
                                eps3,l3,tb,
                                compute)
        #####
        val_3 = compute_Tensor2(eps1,l1,ta,  # Corresponds to ([p1,ta],[p2,tb],[p3,0])
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
# Wrapper function to evaluate Tensor2 with the function compute
# It ensures the exact antisymmetry wrt the first two arguments
#################################################
function compute_Tensor2(eps1::Int64,l1::Int64,t1::Int64,
                         eps2::Int64,l2::Int64,t2::Int64,
                         eps3::Int64,l3::Int64,t3::Int64,
                         compute::F) where {F <: Function}
    #####
    arg1 = eps1, l1, t1 # First argument
    arg2 = eps2, l2, t2 # Second argument
    #####
    if (arg1 == arg2) # The two arguments match
        return 0.0 # The function must vanish
    end
    #####
    # Leveraging the exclusion rule for Tensor2
    if exclusion_Tensor2(eps1,l1,
                         eps2,l2,
                         l3)
        #####
        return 0.0
        #####
    end
    #####
    if (arg1 < arg2) # Sorting the arguments by lexicographic order
        return compute(eps1,l1,t1,
                       eps2,l2,t2,
                       eps3,l3,t3)
    else # We change the order of the two arguments, and multiply by -1
        return -1 * compute(eps2,l2,t2,
                            eps1,l1,t1,
                            eps3,l3,t3)
    end
end
#################################################
# Computing the values of TAB_LAMBDA_BARE
##################################################
function TAB_LAMBDA_BARE!() :: Nothing
    #####
    TAB_TENSOR2!(TAB_LAMBDA_BARE, compute_Lambda_bare)
    #####
    return nothing
end
#################################################
# Computing the values of TAB_LAMBDA
##################################################
function TAB_LAMBDA!() :: Nothing
    #####
    TAB_TENSOR2!(TAB_LAMBDA, compute_Lambda)
    #####
    return nothing
end
