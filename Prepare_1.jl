##################################################
# Computing the Tensors needed in the iteration of Order 1
##################################################
function TAB_PREPARE_1!() :: Nothing
    #####
    # print("TAB_LAMBDA_BARE |")
    TAB_LAMBDA_BARE!() # Computing Lambda_bare
    # print("TAB_LAMBDA |")
    TAB_LAMBDA!() # Computing Lambda
    print("TAB_G_3 |")
    @time TAB_G_3!() # Computing G_3
    # print("TAB_SIGMA |")
    TAB_SIGMA!() # Computing Sigma
    #####
    return nothing
end