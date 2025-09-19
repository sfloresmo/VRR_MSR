##################################################
# Computing new values for TAB_R_NEW and TAB_GAMMA_NEW
# using the values stored in TAB_R and TAB_GAMMA
##################################################
function TAB_NEW!() :: Nothing
    #####
    TAB_PREPARE!() # Computing Tensors as needed
    #####
    TAB_R_NEW!() # Computing the new value of R
    TAB_GAMMA_NEW!() # Computing the new value of Gamma
    #####
    return nothing
end
##################################################
# Function to copy the new arrays
# [TAB_R_NEW,TAB_GAMMA_NEW]
# into the old ones
# [TAB_R,TAB_GAMMA]
##################################################
function TAB_COPY!() :: Nothing
    #####
    TAB_R_COPY!()
    TAB_GAMMA_COPY!()
    #####
    return nothing
end
##################################################
# Performing one iteration of the fix-point search
##################################################
function iterate!() :: Nothing
    #####
    println("Performing iteration | ",get_ITER()+1) #
    #####
    TAB_NEW!() # Computing the new arrays
    TAB_COPY!() # Copying the new arrays into the old ones
    #####
    increase_ITER!() # We have performed an iteration
    #####
    write_generic!(get_ITER()) # Dumping the result
    #####
    return nothing
end
##################################################
# Picking the Prepare function 
##################################################
if     (ORDER == 0)
    include("Prepare_0.jl")
    const TAB_PREPARE! = TAB_PREPARE_0!
elseif (ORDER == 1)
    include("Prepare_1.jl")
    const TAB_PREPARE! = TAB_PREPARE_1!
elseif (ORDER == 2)
    include("Prepare_2.jl")
    const TAB_PREPARE! = TAB_PREPARE_2!
end