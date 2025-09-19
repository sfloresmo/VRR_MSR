##################################################
# Implementing the Elsasser coefficients
# We follow Appendix in Fouvry+2019
##################################################
# Coefficient Lambda appearing in E^L
# See Eq. (56) in Fouvry+2019
##################################################
function _Lambda(l1::Int64,
                 l2::Int64,
                 l3::Int64) :: Float64
    #####
    return sqrt(((2*l1+1)*(2*l2+1)*(2*l3+1))/(4 * pi))
end
##################################################
# Coefficient Delta appearing in E^L
# See Eq. (56) in Fouvry+2019
##################################################
function _Delta(l1::Int64,
                l2::Int64,
                l3::Int64) :: Float64
    #####
    return sqrt(((l1+l2+l3+2)*(l1+l2+l3+4))/(4*(l1+l2+l3+3))) *
           sqrt((l1+l2-l3+1)*(l3+l1-l2+1)*(l2+l3-l1+1))
end
##################################################
# Exclusion rule for the E^L coefficients "the triangular inequality"
# Returns true if the triplet of harmonics has a vanishing value
##################################################
function exclusion_trg(l1::Int64,
                       l2::Int64,
                       l3::Int64) :: Bool
    #####
    if iseven(l1+l2+l3)
        return true
    end
    #####
    if (!(abs(l1-l2) < l3 < (l1+l2)))
        return true
    end
    #####
    return false
end
##################################################
# Returns the E^L coefficients, E^L[l1,l2,l3]
# We load E^L in Gamma_bare only once, initially
##################################################
function get_EL(l1::Int64,
                l2::Int64,
                l3::Int64) :: Float64
    #####
    # We do not compute E^L that do not satisfy the triangular inequality
    # Otherwise, Delta might be ill-defined
    if exclusion_trg(l1,l2,l3)
        return 0.0
    end
    #####
    # Following Eq. (55) of Fouvry+2019,
    # the Elsasser coefficients read
    # E^L = Lambda * Delta * Wigner
    val_Lambda = _Lambda(l1,l2,l3)
    val_Delta  =  _Delta(l1,l2,l3)
    val_Wigner = my_wig3j(l1+1,l2+1,l3+1,0,0,0)
    #####
    val_EL = val_Lambda * val_Delta * val_Wigner # Elsasser coefficients
    #####
    return val_EL # Output
end
#######################################################
# Computing Wigner3j symbols with package CGcoefficient
#######################################################
function my_wig3j(l1::Int64,l2::Int64,l3::Int64,
                  m1::Int64,m2::Int64,m3::Int64) :: Float64
    #####
    return f3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3)
    #####
end
#######################################################
# Computing Wigner6j symbols with package CGcoefficient
#######################################################
function my_wig6j(l1::Int64,l2::Int64,l3::Int64,
                  l4::Int64,l5::Int64,l6::Int64) :: Float64
    #####
    return f6j(2*l1,2*l2,2*l3,2*l4,2*l5,2*l6)
    #####
end

