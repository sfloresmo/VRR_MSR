# ATTENTION, this file must be loaded after Gamma_bare.jl
# so that TAB_GAMMA_BARE is already initialised
##################################################
# Initialising TAB_R with TAB_R_BARE
##################################################
function init_TAB_R_BARE!() :: Nothing
    #####
    for l=0:LMAX
        for t=0:get_Nsteps_G(l)
            #####
            TAB_R[l][t] = TAB_R_BARE[l][t]
            #####
        end
    end
    #####
    return nothing
end
##################################################
# Initialising TAB_R with TAB_R_GAUSSIAN
##################################################
function init_TAB_R_GAUSSIAN!() :: Nothing
    #####
    for l=0:LMAX
        for t=0:get_Nsteps_G(l)
            #####
            TAB_R[l][t] = TAB_R_GAUSSIAN[l][t]
            #####
        end
    end
    #####
    return nothing
end
##################################################
# Initialising TAB_GAMMA with TAB_GAMMA_BARE
# Side-note: We could also initialise with 0.0, I believe?
##################################################
function init_TAB_GAMMA_BARE!() :: Nothing
    #####
    for i_trg = 1:NB_TRG # Loop over the triangles
        #####
        tab_eps = TAB_TRG_EPS[i_trg] # Array of the eps-triplets to consider
        nb_eps = TAB_TRG_NB_EPS[i_trg] # Number of eps-triplets to consider
        #####
        for i_eps = 1:nb_eps # Loop over the eps-triplets
            #####
            eps = tab_eps[i_eps] # Current eps-triplet
            eps1, eps2, eps3 = eps[1], eps[2], eps[3] # Current epsilons
            #####
            val_bare = TAB_GAMMA_BARE[i_trg][i_eps] # Bare value used to initialise
            #####
            for i_t = 1:3 # Loop over the times
                #####
                tab = TAB_GAMMA[i_trg][i_eps][i_t] # Array to be initialised
                #####
                fill!(tab,0.0) # Filling with zero values everywhere
                #####
                tab[0,0] = val_bare # For t1=t2=t3=0, we set to the bare value
                #####
            end
        end
    end
    #####
    return nothing
end
##################################################
# Picking the correct initialisation functions
# For now, we do not have a "Gaussian" initialisation for Gamma.
# Hence, we always use the bare one
##################################################
if (INIT == "Bare")
    const init_TAB_R! = init_TAB_R_BARE!
    const init_TAB_GAMMA! = init_TAB_GAMMA_BARE!
elseif (INIT == "Gaussian")
    const init_TAB_R! = init_TAB_R_GAUSSIAN!
    const init_TAB_GAMMA! = init_TAB_GAMMA_BARE!
end
##################################################
# Initialising TAB_R and TAB_GAMMA
# for ITER = 0
##################################################
function init_from_calc!() :: Nothing
    #####
    init_TAB_R!()
    init_TAB_GAMMA!()
    #####
    write_generic!(0) # Writing the initial file
    #####
    return nothing
end
##################################################
# Initialising TAB_R and TAB_GAMMA
# for ITER > 0
##################################################
function init_from_file!() :: Nothing
    #####
    load_generic!(get_ITER()) # Reading the dumped file
    #####
    return nothing
end
##################################################
# Current iteration number
# We put it in an array so that we can update it
const ITER = [ITER_INIT]
##################################################
# Returns the current iteration number
##################################################
function get_ITER() :: Int64
    return ITER[1]
end
##################################################
# Increases the iteration number
##################################################
function increase_ITER!() :: Nothing
    ITER[1] += 1
    return nothing
end
##################################################
# If we are using ITER = 0,
# then we initialise TAB_R and TAB_GAMMA
# with the initialisation functions
##################################################
if   (ITER_INIT == 0)
    init_from_calc!()
else
    init_from_file!()
end