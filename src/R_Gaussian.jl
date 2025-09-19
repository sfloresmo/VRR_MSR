##################################################
# Definition of the "Gaussian" propagator
# This is used to initialise the fix-point search
# We follow the "Gaussian" expression from Eq. (18) in Fouvry+2019
##################################################
const TAB_R_GAUSSIAN = OffsetArray([OffsetArray([0.0 for t=0:get_Nsteps_G(l)],0:get_Nsteps_G(l)) for l=0:LMAX],0:LMAX) # Gaussian response function
##################################################
# Filling in the Gaussian response function
# R_gaussian[l,t] = exp(- (1/2) * (t/Tl)^2)
##################################################
function TAB_R_GAUSSIAN!() :: Nothing
    #####
    # For l=0,1, we keep a constant R_l
    # Because these harmonics do not decay
    for t=0:get_Nsteps_G(0)
        TAB_R_GAUSSIAN[0][t] = 1.0
    end
    #####
    for t=0:get_Nsteps_G(1)
        TAB_R_GAUSSIAN[1][t] = 1.0
    end
    #####
    for l=2:LMAX # ATTENTION, the loop starts at l=2
        Tl = TAB_TL[l]
        for t=0:get_Nsteps_G(l)
            #####
            t_phys = t * DT # Time in physical units
            #####
            TAB_R_GAUSSIAN[l][t] = exp( - (1/2) * (t_phys / Tl)^2) # Gaussian correlation function
            #####
        end
    end
    #####
    return nothing
end
##################################################
TAB_R_GAUSSIAN!() # Pre-computing the bare response function
