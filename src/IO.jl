const PREFIX = pwd() * "/data/" # Prefix of the files
##################################################
# Returns the namefile for a given iteration
##################################################
function get_namefile(iter::Int64) :: String
    return PREFIX * "data_VRR_MSR_ORDER_"*string(ORDER)*"_LMAX_"*string(LMAX)*
           "_TC_OVER_DT_"*string(TC_OVER_DT)*"_TMAX_OVER_TC_"*string(TMAX_OVER_TC)*
           "_INIT_"*string(INIT)*"_CUT_"*string(CUT)*
           "_ITER_"*string(iter)*".hf5"
end
##################################################
# Writing the generic elements of the data
# in particular the command-line arguments
# and the constants
# ATTENTION, here we create the file, so we use the mode "w" 
##################################################
function write_generic!(iter) :: Nothing
    #####
    namefile = get_namefile(iter) # Name of the dump
    #####
    file = h5open(namefile,"w") # Opening the file and erasing everything if it already existed
    #####
    # Command-line arguments
    write(file,"PARALLEL",PARALLEL)
    write(file,"LMAX",LMAX)
    write(file,"COUPLING",COUPLING)
    write(file,"TC_OVER_DT",TC_OVER_DT)
    write(file,"TMAX_OVER_TC",TMAX_OVER_TC)
    write(file,"ORDER",ORDER)
    write(file,"INIT",INIT)
    write(file,"ITER_INIT",ITER_INIT)
    write(file,"NB_ITER",NB_ITER)
    #####
    # In Constants.jl
    write(file,"G",G)
    write(file,"MBH",MBH)
    write(file,"MTOT",MTOT)
    write(file,"NPART",NPART)
    write(file,"M_STAR",M_STAR)
    write(file,"A_STAR",A_STAR)
    write(file,"ECC_STAR",ECC_STAR)
    write(file,"L_STAR",L_STAR)
    write(file,"H2_COUPLING",H2_COUPLING)
    write(file,"J2_COUPLING",J2_COUPLING)
    write(file,"B2_COUPLING",B2_COUPLING)
    write(file,"C0",C0)
    write(file,"TC",TC)
    #####
    # From Nsteps.jl
    write(file,"TAB_TL",OffsetArrays.no_offset_view(TAB_TL)) # ATTENTION, this is an OffsetArray
    write(file,"DT",DT)
    write(file,"NSTEPS_MAX",NSTEPS_MAX)
    write(file,"TAB_NSTEPS",OffsetArrays.no_offset_view(TAB_NSTEPS)) # ATTENTION, this is an OffsetArray
    write(file,"TAB_T_OVER_TC",OffsetArrays.no_offset_view(TAB_T_OVER_TC)) # ATTENTION, this is an OffsetArray
    #####
    # From Triangles.jl
    write(file,"TAB_TRG",hcat(TAB_TRG...)) # Converting from Vector{Vector{Int64}} to Matrix{Int64}
    write(file,"NB_TRG",NB_TRG)
    write(file,"TAB_TRG_NSTEPS",TAB_TRG_NSTEPS)
    write(file,"TAB_TRG_NB_EPS",TAB_TRG_NB_EPS)
    #####
    # From Triangles.jl
    write(file,"TAB_TRG2",hcat(TAB_TRG2...)) # Converting from Vector{Vector{Int64}} to Matrix{Int64}
    write(file,"NB_TRG2",NB_TRG2)
    write(file,"TAB_TRG2_NSTEPS",TAB_TRG2_NSTEPS)
    #####
    # Dumping TAB_R
    for l=0:LMAX
        write(file,"TAB_R_L_"*string(l),OffsetArrays.no_offset_view(TAB_R[l])) # ATTENTION, this is an OffsetArray
    end
    #####
    # Dumping TAB_GAMMA and TAB_G_3
    for i_trg = 1:NB_TRG
        nb_eps = TAB_TRG_NB_EPS[i_trg] # Number of eps-triplets
        #####
        for i_eps = 1:nb_eps
            for i_t = 1:3
                write(file,"TAB_GAMMA_TRG_"*string(i_trg)*
                            "_EPS_"*string(i_eps)*
                            "_TIME_"*string(i_t),
                            OffsetArrays.no_offset_view(TAB_GAMMA[i_trg][i_eps][i_t])) # ATTENTION, this is an OffsetArray
                write(file,"TAB_G_3_TRG_"*string(i_trg)*
                            "_EPS_"*string(i_eps)*
                            "_TIME_"*string(i_t),
                            OffsetArrays.no_offset_view(TAB_G_3[i_trg][i_eps][i_t])) # ATTENTION, this is an OffsetArray
            end
        end
    end
    #####
    close(file)
    #####
    return nothing
    #####
end
##################################################
# Load from a given iteration
# More precisely, this reads the array TAB_R and TAB_GAMMA
# from the file
# @IMPROVE -- We should add bound checks
##################################################
function load_generic!(iter::Int64) :: Nothing
    #####
    namefile = get_namefile(iter)
    #####
    file = h5open(namefile,"r") # Reading from the file
    #####
    # Loading TAB_R
    for l=0:LMAX
        #####
        Nsteps = get_Nsteps_G(l) # Number of steps for the current l
        #####
        name_R_read = "TAB_R_L_"*string(l) # Name of the dataset to read
        # Reading the dumped tab_R
        # We offset the array so as to not get confused
        tab_R_read = OffsetArray(read(file,name_R_read),0:Nsteps) # Reading the dumped data
        #####
        tab_R_tofill = TAB_R[l] # Name of the array to fill
        #####
        # Filling in the array
        for t=0:Nsteps
            tab_R_tofill[t] = tab_R_read[t]
        end
        #####
    end
    #####
    # Loading TAB_GAMMA
    #####
    for i_trg = 1:NB_TRG # Loop over the triangles
        #####
        nb_eps = TAB_TRG_NB_EPS[i_trg] # Number of eps-triplets
        Nsteps = TAB_TRG_NSTEPS[i_trg] # Number of timesteps
        #####
        for i_eps = 1:nb_eps
            for i_t = 1:3
                #####
                name_Gamma_read = "TAB_GAMMA_TRG_"*string(i_trg)*
                                           "_EPS_"*string(i_eps)*
                                            "_TIME_"*string(i_t)
                #####
                # Reading the dumped tab_Gamma
                # We offset the array so as not to get confused
                tab_Gamma_read = OffsetArray(read(file,name_Gamma_read),0:Nsteps,0:Nsteps)
                #####
                tab_Gamma_tofill = TAB_GAMMA[i_trg][i_eps][i_t] # Name of the array to be filled
                #####
                # Filling in the array
                for tb=0:Nsteps,
                    ta=0:Nsteps
                    #####
                    tab_Gamma_tofill[ta,tb] = tab_Gamma_read[ta,tb]
                    #####
                end
            end
        end
    end
    #####
    close(file)
    #####
    return nothing
end


