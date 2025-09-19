##################################################
# For a given value of LMAX,
# we determine the list of non-trivial "triangles" (l1,l2,l3)
# to be considered when summing over all the interactions
# Here, (l1,l2,l3) is such that
# + l1 <= l2 <= l3
# + exclusion_trg(l1,l2,l3) == false
# Up to transpositions, this corresponds to all the triplets
# (l1,l2,l3) for which Gamma^L[l1,l2,l3] has a chance of being non-zero,
# i.e. all the triplets over which we will need to sum
# In order to reduce the memory imprint of the program,
# we record as few as possible of these triangles
# We use two ways to reduce the number of triangles to be stored:
# + Gamma^L is fully antisymmetric, so that we can always
# sort its arguments, i.e. impose l1 <= l2 <= l3
# + Only triangular triplets need to be considered
# Given these constraints, we need only to consider
# four types of triangles, (l1,l2,l3), namely:
# + (l1,l1,l1)
# + (l1,l1,l2)  with  l1 < l2
# + (l1,l2,l2)  with  l1 < l2
# + (l1,l2,l3)  with  l1 < l2 < l3
# @IMPROVE -- We have not implemented any anti-symmetry wrt time
# @IMPROVE -- There might be a (much) simpler fashion of dealing
#             with these symmetries
##################################################
# Function that constructs the list of all the non-trivial
# triangles given the value of LMAX
##################################################
function get_tab_trg() :: Vector{Vector{Int64}}
    # First we create a list of all the possible triangles
    tab_all = [zeros(Int64,3) for ind=1:(LMAX+1)^3] # ATTENTION, to the size of the array
    #####
    # Filling in the array
    c = 1 # Counter for the array
    for l1=0:LMAX,
        l2=0:LMAX,
        l3=0:LMAX
        #####
        tab_all[c][1] = l1
        tab_all[c][2] = l2
        tab_all[c][3] = l3
        #####
        c += 1
    end
    #####
    # Filtering the array so that we keep only the triplets
    # that comply with the triangular inequality
    tab_trg = filter(x -> !exclusion_trg(x[1],x[2],x[3]),tab_all)
    #####
    # Removing all the pairs that are the same when sorted
    tab_trg = unique(sort,tab_trg)
    #####
    # Sorting this list of triangles, just because I can
    # We now have at our disposal the list
    # of all the non-trivial triangles
    tab_trg = sort(tab_trg)
    #####
    return tab_trg # Output
end
##################################################
const TAB_TRG = get_tab_trg() # Array containing the list of all the non-trivial triangles
const NB_TRG = size(TAB_TRG)[1] # Total number of non-trivial triangles
##################################################
const TAB_TRG_NSTEPS = zeros(Int64,NB_TRG) # Containing the number of timesteps for each triangle
##################################################
# Filling in the array with the number of timesteps
##################################################
function TAB_TRG_NSTEPS!() :: Nothing
    for i_trg = 1:NB_TRG
        trg = TAB_TRG[i_trg]
        #####
        l1, l2, l3 = trg[1], trg[2], trg[3]
        #####
        TAB_TRG_NSTEPS[i_trg] = get_Nsteps_trg(l1,l2,l3)
    end
    #####
    return nothing
end
##################################################
TAB_TRG_NSTEPS!() # Filling in the array
##################################################
# For a given interaction triangle
# (l1,l2,l3) with l1 <= l2 <= l3
# we now want to determine the list of all
# the non-trivial triplets of (eps1,eps2,eps3)
# that need to be considered so that
# any values of Gamma can be determined using antisymmetry
##################################################
# List of non-trivial (eps1,eps2,eps3)
# to consider for a given triangle
# We recall the convention
# +  <-->  1
# -  <-->  2
##################################################
# Triangle of the form (l1,l1,l1)
# We need four triplets in (eps1,eps2,eps3)
# (+,+,+), (+,+,-), (+,-,-), (-,-,-)
##################################################
const TAB_TRG_EPS_L1L1L1 = [[1,1,1],
                            [1,1,2],
                            [1,2,2],
                            [2,2,2]]
##################################################                       
# Triangle of the form (l1,l1,l2) with l1 < l2
# We need six triplets in (eps1,eps2,eps3)
# (+,+,+), (+,+,-), (+,-,+), (+,-,-), (-,-,+), (-,-,-)
################################################## 
const TAB_TRG_EPS_L1L1L2 = [[1,1,1],
                            [1,1,2],
                            [1,2,1],
                            [1,2,2],
                            [2,2,1],
                            [2,2,2]]
##################################################                       
# Triangle of the form (l1,l2,l2) with l1 < l2
# We need six triplets in (eps1,eps2,eps3)
# (+,+,+), (+,+,-), (+,-,-), (-,+,+), (-,+,-), (-,-,-)
##################################################
const TAB_TRG_EPS_L1L2L2 = [[1,1,1],
                            [1,1,2],
                            [1,2,2],
                            [2,1,1],
                            [2,1,2],
                            [2,2,2]]
##################################################                       
# Triangle of the form (l1,l2,l3) with l1 < l2 < l3
# We need eight triplets in (eps1,eps2,eps3)
# (+,+,+), (+,+,-), (+,-,+), (+,-,-), (-,+,+), (-,+,-), (-,-,+), (-,-,-)
##################################################
const TAB_TRG_EPS_L1L2L3 = [[1,1,1],
                            [1,1,2],
                            [1,2,1],
                            [1,2,2],
                            [2,1,1],
                            [2,1,2],
                            [2,2,1],
                            [2,2,2]]
##################################################
# Array containing all the permutations on [1,2,3]
# This array is used to find the appropriate permutation
# of the eps-triplet, once the triangle (l1,l2,l3) has been fixed
# Here, permutations are defined with the same convention
# as in TAB_TRG_EPS_PERM, i.e.
# Given an initial array, "old", and sorted one, "new",
# the permutation, "p", associated with this sorting is such that
# new[i] = old[p[i]]
##################################################
const TAB_PERM_LIST = [[1,2,3],
                       [1,3,2],
                       [2,1,3],
                       [2,3,1],
                       [3,1,2],
                       [3,2,1]]
const NB_PERM_LIST = size(TAB_PERM_LIST)[1] # Total number of permutations
##################################################
# Function that for a given interaction triangle
# returns the list of the non-trivial
# (eps1,eps2,eps3) to be considered
##################################################
function get_trg_eps(trg::Vector{Int64}) :: Vector{Vector{Int64}}
    #####
    l1, l2, l3 = trg[1], trg[2], trg[3] # Extracting the interaction triangle
    #####
    # For safety, we check that it is indeed an appropriate interaction triangle
    exclusion_trg(l1,l2,l3) && error("get_trg_eps: Not a triangle")
    !(l1 <= l2 <= l3)       && error("get_trg_eps: Not ordered")
    #####
    # Triangle of the form (l1,l1,l1)
    if (l1 == l2 == l3)
        return TAB_TRG_EPS_L1L1L1
    end
    #####
    # Triangle of the form (l1,l1,l2)
    if (l1 == l2)
        return TAB_TRG_EPS_L1L1L2
    end
    #####
    # Triangle of the form (l1,l2,l2)
    if (l2 == l3)
        return TAB_TRG_EPS_L1L2L2
    end
    #####
    # Triangle of the form (l1,l2,l3)
        return TAB_TRG_EPS_L1L2L3
end
##################################################
# Creating an array of length NB_TRG
# containing the list of non-trivial (eps1,eps2,eps3)
# to be considered for the triangles given by TAB_TRG
##################################################
const TAB_TRG_EPS = [get_trg_eps(TAB_TRG[i_trg]) for i_trg = 1:NB_TRG]
const TAB_TRG_NB_EPS = [size(TAB_TRG_EPS[i_trg])[1] for i_trg = 1:NB_TRG] # Number of non-trivial eps-triplets for each triangle
##################################################
# We now prepare the ground to allow for the easy evaluation
# of Gamma using the symmetries at our disposal
# For a given location of evaluation, namely
# ([eps1,l1,t1],[eps2,l2,t2],[eps3,l3,t3]),
# we need to determine
# + the index of the triangle of harmonics associated
# [to be used in TAB_TRG]
# + the index of the associated eps-triplets
# [to be used in TAB_TRG_EPS]
# + the prefactor (+/-1) by which the result must be multiplied
# [this corresponds to the parity of the number of inversions, see Utils.jl.]
# ATTENTION, we set this number to 0 for non-triangular harmonics
# + the permutation used to obtain the rewriting
# [this corresponds to a 3-vector, perm, of Int64, between 1 and 3,
# so that t1_new = [t1,t2,t3][perm[1]], t2_new = [t1,t2,t3][perm[2]], ...] 
##################################################
# In order to accelerate the evaluation of the function,
# we store all these quantities in an array of the form
# [eps1,l1,eps2,l2,eps3,l3]
# where each element reads
# {ind_trg,ind_eps,sign_perm,perm_1,perm_2,perm_3}
# For now, we create a temporary array,
# that we make static afterwards
const TAB_TRG_EPS_ALL_TEMP = OffsetArray([zeros(Int64,6)
                        for eps1=1:2,l1=0:LMAX,
                            eps2=1:2,l2=0:LMAX,
                            eps3=1:2,l3=0:LMAX],
                            1:2,0:LMAX,
                            1:2,0:LMAX,
                            1:2,0:LMAX)
##################################################
# Filling in the array containing the symmetrising mapping
# @IMPROVE -- I did not care at all about performance and allocations
##################################################
function TAB_TRG_EPS_ALL!()
   #####
   # Loop over all the harmonics and spins
    for l1=0:LMAX,
        l2=0:LMAX,
        l3=0:LMAX,
        eps1=1:2,
        eps2=1:2,
        eps3=1:2
        #####
        # Checking that the index is triangular
        # Then, we need to determine (ind_trg,ind_eps,sign,perm)
        if !exclusion_trg(l1,l2,l3)
            #####
            # First, we need to find the triangle
            # to which the current (l1,l2,l3) is associated
            # To do so, we simply scan over TAB_TRG
            # and find the first occurrence where the sets of harmonics are the same
            # @IMPROVE -- This is a lazy implementation
            #####
            # Initialising the index
            # We start with zero, so that we would spot bugs
            ind_trg = 0
            my_trg = [l1,l2,l3] # Current harmonics
            my_trg_sorted = sort(my_trg) # Sorting the harmonics
            #####
            # Loop over the triangles to try and find the matching one
            for ind_trg_current = 1:NB_TRG
                current_trg = TAB_TRG[ind_trg_current] # Current triangle
                if (my_trg_sorted == current_trg) # The two triangles match
                    ind_trg = ind_trg_current # We have found the position of the triangle
                    break # Exiting the for loop
                end
            end
            #####
            # We have found the position of the associated triangle
            # We now determine the associated permutation
            # The array perm_trg is such that
            # the element "i" of the sorted array
            # is in position perm_trg[i] of the initial array
            # ATTENTION, given this convention,
            # we must use the julia function 'sortperm' on the initial array
            perm_trg = sortperm(my_trg)
            #####
            # We now act with this permutation on the input
            # (eps1,eps2,eps3)
            my_eps = [eps1,eps2,eps3]
            my_eps_permuted = my_eps[perm_trg]
            #####
            # From the permuted my_eps_permuted,
            # we now need to determine the associated
            # non-trivial eps-triplets
            # This is delicate.
            # Indeed, we need to find a permutation, perm_eps, such that
            #  + my_trg_sorted is left unchanged by perm_eps [ATTENTION, this is for 'my_trg_sorted']
            #  + my_eps_permuted becomes an element of TAB_TRG_EPS[ind_trg] by perm_eps
            # To find this permutation, we simply try every possible permutations,
            # as stored in TAB_PERM_LIST, and stop as soon as once we find a permutation
            # that conserves my_trg and makes my_eps_permuted an elements of TAB_TRG_EPS[ind_trg]
            # @IMPROVE -- It must be possible to write a smarter/simpler code.
            #####
            nb_eps = TAB_TRG_NB_EPS[ind_trg] # Number of eps-triplets to compare with
            tab_eps = TAB_TRG_EPS[ind_trg] # List of the non-trivial eps-triplets to compare with
            #####
            # Initialising the indices
            ind_perm_eps = 0 # Index of the permutation, in TAB_PERM_LIST. We are searching for this
            ind_eps = 0 # Index of the non-trivial eps-triplet, in tab_eps. We are searching for this
            #####
            # Loop over the possible permutations
            for ind_perm_eps_current = 1:NB_PERM_LIST # Loop over all possible permutations
                #####
                perm_loc = TAB_PERM_LIST[ind_perm_eps_current] # Current permutation that we are trying
                #####
                # We now act with this permutation over the triangle
                # as well as over epsilon
                trg_loc = my_trg_sorted[perm_loc] # Permuting the triangle of harmonics. ATTENTION, this is for 'my_trg_sorted'
                eps_loc = my_eps_permuted[perm_loc] # Attempting the eps-triplet. ATTENTION, this is for 'my_eps_permuted'
                #####
                # We now need to check that
                # (i)  trg_loc == my_trg_sorted, i.e the triangle has not changed [ATTENTION, this is for 'my_trg_sorted']
                # (ii) eps_loc appears within tab_eps
                if (trg_loc == my_trg_sorted) # The triangle has not changed. ATTENTION, this is for 'my_trg_sorted'
                    #####
                    # We now try to see if we can find eps_loc within tab_eps
                    for ind_eps_current = 1:nb_eps # Loop over all the non-trivial eps
                        current_eps = tab_eps[ind_eps_current] # Current eps
                        if (eps_loc == current_eps) # The two eps-triplets match
                            #####
                            # Fantastic, we have found a match
                            ind_perm_eps = ind_perm_eps_current # Index of the appropriate permutation in TAB_PERM_LIST
                            ind_eps = ind_eps_current # Index of the appropriate eps-triplet in tab_eps
                            #####
                            break # Leaving the eps-loop
                        end
                    end
                end
                #####
                # If we have found a solution,
                # we may leave the loop over the permutations
                if (ind_perm_eps != 0) && (ind_eps != 0) # The indices have been updated
                    break # Leaving the loop over permutations
                end
                #####
            end
            #####
            perm_eps = TAB_PERM_LIST[ind_perm_eps] # Appropriate permutation of the eps-triplet
            #####
            # We now construct the total permutation
            # It is composed of the composition
            # of the permutation perm_trg then perm_eps
            # The final perm_tot is defined such that
            # the element "i" of the sorted array
            # is in position perm_tot[i] of the initial array
            # Let us consider a sequence of sorting of array
            # A  -->_{a}  B  -->_{b}  C with the associated 'perm', a and b.
            # With our convention, we have
            #   + B[i] = A[a[i]]
            #   + C[i] = B[b[i]]
            # In order to find the total perm, we want to find
            # a permutation, c, such that
            #   + C[i] = A[c[i]]
            # Injecting the two previous relations into one another, we find
            #   + C[i] = B[b[i]] = A[a[b[i]]]
            # So, we ultimately find
            #   + c[i] = a[b[i]]
            # ATTENTION, this is a very delicate step.
            #####
            perm_tot = perm_trg[perm_eps] # Constructing the total permutation
            #####
            sign_perm = get_signature(perm_tot) # Signature of the permutation
            #####
            # Filling in the temporary array
            tab_sym_loc = TAB_TRG_EPS_ALL_TEMP[eps1,l1,eps2,l2,eps3,l3] # Pointer to the array to be filled
            tab_sym_loc[1] = ind_trg
            tab_sym_loc[2] = ind_eps
            tab_sym_loc[3] = sign_perm
            tab_sym_loc[4] = perm_tot[1]
            tab_sym_loc[5] = perm_tot[2]
            tab_sym_loc[6] = perm_tot[3]
        else
            #####
            # The harmonics considered are not triangular.
            # So, we fill in the arrays with "bad" values,
            # so that they would be caught by the code, if ever accessed
            ind_trg = 0 # Bad index outside of the range
            ind_eps = 0 # Bad index outside of the range
            sign_perm = 0 # Bad signature (Should only be +/- 1)
            perm_tot = 0, 0, 0 # Bad permutation (Should be a permutation of [1,2,3])
            #####
            # Filling in the temporary array
            tab_sym_loc = TAB_TRG_EPS_ALL_TEMP[eps1,l1,eps2,l2,eps3,l3] # Pointer to the array to be filled
            tab_sym_loc[1] = ind_trg
            tab_sym_loc[2] = ind_eps
            tab_sym_loc[3] = sign_perm
            tab_sym_loc[4] = perm_tot[1]
            tab_sym_loc[5] = perm_tot[2]
            tab_sym_loc[6] = perm_tot[3]
        end
    end
    #####
end
##################################################
TAB_TRG_EPS_ALL!() # Preparing all the symmetrising arrays
###################################################
# Now, we can create a static version of the array
# to allow for even faster readings
# The elements are Tuple, i.e. statically sized
const TAB_TRG_EPS_ALL = OffsetArray([ Tuple(TAB_TRG_EPS_ALL_TEMP[eps1,l1,eps2,l2,eps3,l3])
                        for eps1=1:2,l1=0:LMAX,
                            eps2=1:2,l2=0:LMAX,
                            eps3=1:2,l3=0:LMAX],
                            1:2,0:LMAX,
                            1:2,0:LMAX,
                            1:2,0:LMAX)