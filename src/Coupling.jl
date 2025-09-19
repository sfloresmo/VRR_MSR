##################################################
# Definition of the coupling coefficients
# ATTENTION, we assume that we are in the single-population case
##################################################
# Coupling coefficients for COUPLING=Quad
##################################################
function get_Jl_Quad(l::Int64) :: Float64
    (l == 2) ? J2_COUPLING : 0.0
end
##################################################
# Picking the appropriate coupling coefficient
##################################################
if (COUPLING == "Quad")
    const get_Jl = get_Jl_Quad
end
