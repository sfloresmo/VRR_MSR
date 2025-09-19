##################################################
# Computing the Tensors needed in the iteration of Order 0
##################################################
function TAB_PREPARE_0!() :: Nothing
    #####
    TAB_LAMBDA_BARE!() # Computing Lambda_bare
    TAB_LAMBDA!() # Computing Lambda
    TAB_SIGMA!() # Computing Sigma
    #####
    return nothing
end