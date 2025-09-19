##################################################
# Value of the diagram W3
# The value of this diagram depends on l1, ..., l6
##################################################
function get_W3(l1::Int64,l2::Int64,l3::Int64,
                l4::Int64,l5::Int64,l6::Int64) :: Float64
    #####
    # Checking triangles
    #####
    if exclusion_trg(l1,l4,l5) || exclusion_trg(l2,l4,l6) || exclusion_trg(l3,l5,l6)
        return 0.0
    end
    #####
    res = (-1)^(1 + l4 + l5 + l6) * my_wig6j(l1,l2,l3,l6,l5,l4)
    #####
    return res # Output
end

