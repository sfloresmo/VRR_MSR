##################################################
# Value of the diagram W2
##################################################
function get_W2(l1::Int64, l2::Int64, l3::Int64) :: Float64
    #####
    # Checking triangle
    #####
    if exclusion_trg(l1,l2,l3)
        return 0.0
    end
    #####
    return 1 / (2 * l1 + 1)
end