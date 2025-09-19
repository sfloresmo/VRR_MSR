##################################################
# Returns the signature of a permutation over 3 elements
# We do this by computing the number of inversions, see
# https://en.wikipedia.org/wiki/Parity_of_a_permutation
##################################################
function get_signature(perm::Vector{Int64}) :: Int64
    n = length(perm) # Length of the permutation -- In practice, will be 3
    c = 0 # Counter for the number of inversions
    for i=1:n
        for j=(i+1):n
            if (perm[i] > perm[j]) # We have found an inversion
                c += 1 # Updating the counter
            end
        end
    end
    #####
    if iseven(c)
        return 1
    else
        return -1
    end
end
#################################################
# Simple function used to swap two arguments
#################################################
function swap(arg1,arg2)
    return arg2, arg1
end
#################################################
# Function to sort efficiently a sequence of three t-uples
# We use a sorting network
# cf https://bertdobbelaere.github.io/sorting_networks.html#N3L3D3
# We also return the signature of the permutation
# ATTENTION, it (may) assume that all the arguments are different
#################################################
function sort3(arg1,arg2,arg3)
    #####
    # Copying the inputs
    arg1_loc = arg1
    arg2_loc = arg2
    arg3_loc = arg3
    #####
    perm_sign = 1 # Initialising the signature of the permutation
    #####
    if (arg1_loc > arg3_loc)
        perm_sign *= -1
        arg1_loc, arg3_loc = swap(arg1_loc,arg3_loc)
    end
    #####
    if (arg1_loc > arg2_loc)
        perm_sign *= -1
        arg1_loc, arg2_loc = swap(arg1_loc,arg2_loc)
    end
    #####
    if (arg2_loc > arg3_loc)
        perm_sign *= -1
        arg2_loc, arg3_loc = swap(arg2_loc,arg3_loc)
    end
    #####
    return arg1_loc, arg2_loc, arg3_loc, perm_sign
    #####
end