# Reading 2-point correlation function
function read_prediction2(l::Int64, iter::Int64)
    # Name of the file
    namefile = "data/data_VRR_MSR_ORDER_"*string(ORDER)*"_LMAX_"*string(LMAX)*
           "_TC_OVER_DT_"*string(TC_OVER_DT)*"_TMAX_OVER_TC_"*string(TMAX_OVER_TC)*
           "_INIT_"*string(INIT)*"_CUT_"*string(CUT)*
           "_ITER_"*string(iter)*".hf5"
    
    #####
    file = h5open(namefile,"r")

    # Rescaled time with respect to Tc
    tabt = read(file["TAB_T_OVER_TC"])
    
    # Reading the two-point correlation function
    tabCl2 = read(file["TAB_R_L_$(l)"])
    
    # Output
    return tabt, tabCl2

    #####
    close(file)
end

# Reading 3-point correlation function for triangle (l,l,l) = (1,2,2)
function read_prediction3(iter::Int64)
    # Name of the file
    namefile = "data/data_VRR_MSR_ORDER_"*string(ORDER)*"_LMAX_"*string(LMAX)*
           "_TC_OVER_DT_"*string(TC_OVER_DT)*"_TMAX_OVER_TC_"*string(TMAX_OVER_TC)*
           "_INIT_"*string(INIT)*"_CUT_"*string(CUT)*
           "_ITER_"*string(iter)*".hf5"
    
    #####
    file = h5open(namefile,"r")
    
    # Reading the two-point correlation function
    i_trg = 2 #trianlge (l,l,l) = (1,2,2)
    i_eps = 1 #spins (+ + +)
    i_t = 1 #times (0,ta,tb)
    tabCl3 = read(file["TAB_G_3_TRG_"*string(i_trg)*
                            "_EPS_"*string(i_eps)*
                            "_TIME_"*string(i_t)])
    
    # Output
    return ((4*pi)^1.5) * tabCl3 #added skewness normalisation

    #####
    close(file)
end

# Reading the 2-point prediction for the harmonics l=2, 3, 4
time2, values2 = read_prediction2(2,NB_ITER)
time3, values3 = read_prediction2(3,NB_ITER)
time4, values4 = read_prediction2(4,NB_ITER)

# Reading 3-point correlation function for triangle (l,l,l) = (1,2,2)
values3point = read_prediction3(NB_ITER)


