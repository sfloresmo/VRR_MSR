# julia Run_Create_Figs.jl --lmax 7 --coupling Quad --Tc_over_dt 10 --Tmax_over_Tc 1 --order 1 --init Gaussian --iter_init 0 --nb_iter 10
##################################################
include("../src/Main.jl")
include("../src/Create_Figs.jl")
##################################################

# Creating the plot 2-point correlation
plot(time2, values2, label="l=2", xlabel="t/Tc", ylabel="C_l/C0", title="Normalised Two-point Correlation")
plot!(time3, values3, label="l=3")
plot!(time4, values4, label="l=4")

# Save the plot 
figs_directory = string("figs/")
savefig(figs_directory * "figure_two_point_ORDER_"*string(ORDER)*"_LMAX_"*string(LMAX)*"_ITER_"*string(NB_ITER)*".hf5.png")

# Creating the plot 3-point correlation
plot() 
contour!(values3point, fill=true, title="Normalised Three-point Correlation (skewness)", xticks=false, yticks=false)

# Save the plot 
figs_directory = string("figs/")
savefig(figs_directory * "figure_three_point_ORDER_"*string(ORDER)*"_LMAX_"*string(LMAX)*"_ITER_"*string(NB_ITER)*".hf5.png")