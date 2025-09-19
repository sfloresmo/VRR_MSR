##################################################
# For a given value of LMAX,
# we determine the list of non-trivial "triangles2" (l1,l2,l3)
# to be considered when constructing the tensor Lambda and Omega
# ATTENTION, these tensors are antisymmetric only wrt their two first indices
# Here, (l1,l2,l3) is such that
# + l1 <= l2 <= l3
# + exclusion_trg(l1,l2,l3) == false
# Up to transpositions, this corresponds to all the triplets
# (l1,l2,l3) for which Lambda^L[l1,l2,l3]/Omega^L[l1,l2,l3] has a chance of being non-zero,
# i.e. all the triplets over which we will need to sum
# In order to reduce the memory imprint of the program,
# we record as few as possible of these triangles2
# We use two ways to reduce the number of triangles2 to be stored:
# + Lambda^L/Omega^L is antisymmetric wrt its two first indices.
# As a consequence, we can always sort the two first arguments,
# i.e. impose l1 <= l2
# + Only triangular triplets need to be considered
# Given these constraints, we need only to consider
# two types of triangles2
# + (l1,l1,l3)  [no constraints on l3]
# + (l1,l2,l3)  [l1 < l2, no constraints on l3]
##################################################
# Function that constructs the list of all the non-trivial
# triangles2 given the value of LMAX
##################################################
function get_tab_trg2() :: Vector{Vector{Int64}}
    # First we create a list of all the possible triangles2
    tab_all2 = [zeros(Int64,3) for ind=1:(LMAX+1)^3] # ATTENTION, to the size of the array
    #####
    # Filling in the array
    c = 1 # Counter for the array
    for l1=0:LMAX,
        l2=0:LMAX,
        l3=0:LMAX
        #####
        tab_all2[c][1] = l1
        tab_all2[c][2] = l2
        tab_all2[c][3] = l3
        #####
        c += 1
    end
    #####
    # Filtering the array so that we keep only the triplets
    # that comply with the triangular inequality
    tab_trg2 = filter(x -> !exclusion_trg(x[1],x[2],x[3]),tab_all2)
    #####
    # Removing all the pairs that are the same
    # when the two first arguments are sorted
    tab_trg2 = unique(x -> [min(x[1],x[2]),max(x[1],x[2]),x[3]],tab_trg2)
    #####
    # Sorting this list of triangles2, just because I can.
    # We now have at our disposal the list
    # of all the non-trivial triangles2
    tab_trg2 = sort(tab_trg2)
    #####
    return tab_trg2 # Output
end
##################################################
const TAB_TRG2 = get_tab_trg2() # Array containing the list of all the non-trivial triangles2
const NB_TRG2 = size(TAB_TRG2)[1] # Total number of non-trivial triangles2
##################################################
const TAB_TRG2_NSTEPS = zeros(Int64,NB_TRG2) # Containing the number of timesteps for each triangle
##################################################
# Filling in the array with the number of timesteps
##################################################
function TAB_TRG2_NSTEPS!() :: Nothing
    for i_trg2 = 1:NB_TRG2
        trg2 = TAB_TRG2[i_trg2]
        #####
        l1, l2, l3 = trg2[1], trg2[2], trg2[3]
        #####
        TAB_TRG2_NSTEPS[i_trg2] = get_Nsteps_trg(l1,l2,l3)
    end
    #####
    return nothing
end
##################################################
TAB_TRG2_NSTEPS!() # Filling in the array
##################################################
# For a given interaction triangle2
# (l1,l2,l3) with l1 <= l2 [no constraints on l3]
# we now want to determine the list of all the non-trivial
# triplets of (eps1,eps2,eps3) that need to be considered
# We recall the convention
# +  <-->  1
# -  <-->  2
##################################################
# Triangle2 of the form (l1,l1,l3) [no constraints on l3]
# We need six triplets in (eps1,eps2)
# (+,+,+), (+,+,-), (+,-,+), (+,-,-), (-,-,+), (-,-,-)
##################################################
const TAB_TRG2_EPS_L1L1L3 = [[1,1,1],
                             [1,1,2],
                             [1,2,1],
                             [1,2,2],
                             [2,2,1],
                             [2,2,2]]
##################################################
# Triangle2 of the form (l1,l2,l3) [l1 < l2, no constraints on l3]
# We need eight triplets in (eps1,eps2,eps3)
# (+,+,+), (+,+,-), (+,-,+), (+,-,-), (-,+,+), (-,+,-), (-,-,+), (-,-,-)
##################################################
const TAB_TRG2_EPS_L1L2L3 = [[1,1,1],
                             [1,1,2],
                             [1,2,1],
                             [1,2,2],
                             [2,1,1],
                             [2,1,2],
                             [2,2,1],
                             [2,2,2]]
##################################################
# Function that for a given triangle2
# returns the list of the non-trivial (eps1,eps2,eps3) to be considered
##################################################
function get_trg2_eps(trg2::Vector{Int64}) :: Vector{Vector{Int64}}
    #####
    l1, l2, l3 = trg2[1], trg2[2], trg2[3] # Extracting the triangle2
    #####
    # For safety, we check that it is indeed an appropriate interaction triangle2
    exclusion_trg(l1,l2,l3) && error("get_trg2_eps: Not a triangle2")
    !(l1 <= l2)             && error("get_trg2_eps: Not ordered in two first arguments")
    #####
    # Triangle2 of the form (l1,l1,l3) [no constraints on l3]
    if (l1 == l2)
        return TAB_TRG2_EPS_L1L1L3
    end
    #####
    # Triangle2 of the form (l1,l2,l3) [with l1 < l2, no constraints on l3]
    return TAB_TRG2_EPS_L1L2L3
end
##################################################
# Creating an array of length NB_TRG2
# containing the list of the non-trivial (eps1,eps2,eps3)
# to be considered for the triangles2 given by TAB_TRG2
##################################################
const TAB_TRG2_EPS = [get_trg2_eps(TAB_TRG2[i_trg2]) for i_trg2 = 1:NB_TRG2]
const TAB_TRG2_NB_EPS = [size(TAB_TRG2_EPS[i_trg2])[1] for i_trg2 = 1:NB_TRG2] # Number of non-trivial eps-triplets for each triangle2
##################################################
# We now prepare the ground for the easy evaluation of Lambda and Omega
# using the symmetries at our disposal
# For a given location of evaluation, namely
# ([eps1,l1,t1],[eps2,l2,t2],[eps3,l3,t3])
# we need to determine:
# + the index of the triangle2 of harmonics associated
# [to be used in TAB_TRG2]
# + the index of the associated eps-triplets
# [to be used in TAB_TRG2_EPS]
# + the prefactor (+/-1) by which the result must be multiplied
# [this corresponds to the parity of the number of inversions]
# ATTENTION, we set this number to 0 for non-triangular harmonics
# the permutation used to obtain the rewriting
# [this corresponds to a 3-vector, perm, of Int64, between 1 and 3,
# so that t1_new = [t1,t2,t3][perm[1]], t2_new = [t1,t2,t3][perm[2]], ...] 
##################################################
# In order to accelerate the evaluation of the function,
# we store all these quantities in an array of the form
# [eps1,l1,eps2,l2,eps3,l3]
# where each element reads
# {ind_trg2,ind_eps,sign_perm,perm_1,perm_2,perm_3}
# For now, we create a temporary array,
# that we make static afterwards
const TAB_TRG2_EPS_ALL_TEMP = OffsetArray([zeros(Int64,6)
                        for eps1=1:2,l1=0:LMAX,
                            eps2=1:2,l2=0:LMAX,
                            eps3=1:2,l3=0:LMAX],
                            1:2,0:LMAX,
                            1:2,0:LMAX,
                            1:2,0:LMAX)
##################################################
# Filling the array containing the symmetrising mappings
# @IMPROVE -- I did not care at all about performance and allocations
##################################################
function TAB_TRG2_EPS_ALL!()
    #####
    # Loops over all the harmonics and spins
    for l1=0:LMAX,
        l2=0:LMAX,
        l3=0:LMAX,
        eps1=1:2,
        eps2=1:2,
        eps3=1:2
        #####
        # Checking that the index is triangular
        # Then, we need to determine (ind_trg2,ind_eps,sign,perm)
        #####
        if !exclusion_trg(l1,l2,l3)
            #####
            # First, we need to find the triangle2 to which
            # the current (l1,l2,l3) is associated
            # To do so, we simply scan over TAB_TRG2
            # and find the first occurrence where the sets of harmonics are the same
            # @IMPROVE -- This is a lazy implementation
            #####
            # Initialising the index
            # We start with zero, so that we would spot bugs
            ind_trg2 = 0
            my_trg2 = [l1,l2,l3] # Current harmonics
            #####
            # Sorting the array wrt to the two first indices
            # We also record whether we had to perform a permutation
            if (l1 <= l2)
                my_trg2_sorted = [l1,l2,l3]
            else
                my_trg2_sorted = [l2,l1,l3]
            end
            #####
            # Loop over the triangles2 to find the matching one
            for ind_trg2_current = 1:NB_TRG2
                current_trg2 = TAB_TRG2[ind_trg2_current] # Current triangle2
                if (my_trg2_sorted == current_trg2) # The two triangles2 match
                    ind_trg2 = ind_trg2_current # We have found the position of the triangle2
                    break # Exiting the for loop
                end
            end
            #####
            # We have found the position of the associated triangle2
            # We now determine the associated permutation
            # The array perm_trg2 is such that
            # the element "i" of the sorted array
            # is in position perm_trg2[i] of the initial array
            if (l1 <= l2) # No permutation performed
                perm_trg2 = [1,2,3]
            else # We performed a permtuation
                perm_trg2 = [2,1,3]
            end
            #####
            # We now act with this permutation on the input
            # (eps1,eps2,eps3)
            my_eps = [eps1,eps2,eps3]
            my_eps_permuted = my_eps[perm_trg2]
            #####
            # From the permuted my_eps_permuted,
            # we now need to determine the associated non-trivial eps-triplet
            # We need to find a permutation, perm_eps, such that
            # + my_trg2_sorted is left unchanged by perm_eps [ATTENTION, this is for 'my_trg2_sorted']
            # + my_eps_permuted becomes an element of TAB_TRG2_EPS[ind_trg2] by perm_eps
            # To find this permutation, we simply try every possible permutations,
            # as stored in TAB_PERM_LIST, and stop as soon as once we find a permutation
            # that conserves my_trg2 and makes my_eps_permuted an elements of TAB_TRG2_EPS[ind_trg2]
            # @IMPROVE -- It must be possible to write a smarter/simpler code.
            #####
            nb_eps = TAB_TRG2_NB_EPS[ind_trg2] # Number of eps-triplets to compare with
            tab_eps = TAB_TRG2_EPS[ind_trg2] # List of the non-trivial eps-triplets to compare with
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
                # We now act with this permutation over the triangle2
                # as well as over epsilon
                trg2_loc = my_trg2_sorted[perm_loc] # Permuting the triangle2 of harmonics. ATTENTION, this is for 'my_trg2_sorted'
                eps_loc = my_eps_permuted[perm_loc] # Attempting the eps-triplet. ATTENTION, this is for 'my_eps_permuted'
                #####
                # We now need to check that
                # (i)  trg2_loc == my_trg2_sorted, i.e the triangle2 has not changed [ATTENTION, this is for 'my_trg2_sorted']
                # (ii) eps_loc appears within tab_eps
                if (trg2_loc == my_trg2_sorted) # The triangle2 has not changed. ATTENTION, this is for 'my_trg2_sorted'
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
            # of the permutation perm_trg2 then perm_eps
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
            perm_tot = perm_trg2[perm_eps] # Constructing the total permutation
            #####
            sign_perm = get_signature(perm_tot) # Signature of the permutation
            #####
            # Filling in the temporary array
            tab_sym_loc = TAB_TRG2_EPS_ALL_TEMP[eps1,l1,eps2,l2,eps3,l3] # Pointer to the array to be filled
            tab_sym_loc[1] = ind_trg2
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
            ind_trg2 = 0 # Bad index outside of the range
            ind_eps = 0 # Bad index outside of the range
            sign_perm = 0 # Bad signature (Should only be +/- 1)
            perm_tot = 0, 0, 0 # Bad permutation (Should be a permutation of [1,2,3])
            #####
            # Filling in the temporary array
            tab_sym_loc = TAB_TRG2_EPS_ALL_TEMP[eps1,l1,eps2,l2,eps3,l3] # Pointer to the array to be filled
            tab_sym_loc[1] = ind_trg2
            tab_sym_loc[2] = ind_eps
            tab_sym_loc[3] = sign_perm
            tab_sym_loc[4] = perm_tot[1]
            tab_sym_loc[5] = perm_tot[2]
            tab_sym_loc[6] = perm_tot[3]
        end
    end
end
##################################################
TAB_TRG2_EPS_ALL!() # Preparing all the symmetrising arrays
##################################################
# Now, we can create a static version of the array
# to allow for even faster readings
# The elements are Tuple, i.e. statically sized
const TAB_TRG2_EPS_ALL = OffsetArray([ Tuple(TAB_TRG2_EPS_ALL_TEMP[eps1,l1,eps2,l2,eps3,l3])
                        for eps1=1:2,l1=0:LMAX,
                            eps2=1:2,l2=0:LMAX,
                            eps3=1:2,l3=0:LMAX],
                            1:2,0:LMAX,
                            1:2,0:LMAX,
                            1:2,0:LMAX)

