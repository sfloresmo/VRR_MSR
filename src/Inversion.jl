##################################################
# We perform the inverstion
# of a lower triangular Toeplitz matrix
# Its inverse is a lower triangular Toeplitz matrix
# Both matrices can be specified by a single vector,
# namely the first column of the matrix
##################################################
# Function to invert TAB_INVR
# and dumping the result in TAB_R_NEW
# @IMPROVE -- Code should be factored
# @IMPROVE -- Implement a faster algorithm [most likely un-needed]
##################################################
function INV_TAB_INVR!() :: Nothing
    for l=0:LMAX
        inv_0 = 1 / TAB_INVR[l][0] # Inverting the first element
        TAB_R_NEW[l][0] = inv_0 / DT^(2) # Filling in the array. ATTENTION to the unit
        #####
        for k=1:get_Nsteps_G(l) # Iteration over the other elements
            s = 0.0 # Initialising the sum
            for i=1:k
                s += TAB_R_NEW[l][k-i] * TAB_INVR[l][i] # ATTENTION, no need for a change of unit here
            end
            TAB_R_NEW[l][k] = - inv_0 * s # ATTENTION, no need for a change of unit here
        end
    end
    #####
    return nothing
end
#Sofia: check units when inverting
