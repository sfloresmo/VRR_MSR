##################################################
# Computing the Tensors needed in the iteration of Order 2
##################################################
function TAB_PREPARE_2!() :: Nothing
    #####
    # print("TAB_LAMBDA_BARE |")
    TAB_LAMBDA_BARE!() # Computing Lambda_bare
    # print("TAB_LAMBDA |")
    TAB_LAMBDA!() # Computing Lambda
    print("TAB_THETA |")
    @time TAB_THETA!() # Computing Theta
    print("TAB_G_3 |")
    @time TAB_G_3!() # Computing G_3
    # print("TAB_SIGMA |")
    TAB_SIGMA!() # Computing Sigma
    #####
    return nothing
end
