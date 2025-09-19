# julia -t 8 --check-bounds=no run/Run.jl --parallel true --lmax 10 --coupling Quad --Tc_over_dt 10 --Tmax_over_Tc 1 --order 2 --init Bare --iter_init 0 --nb_iter 10 --cut false
##################################################
include("../src/Main.jl")
##################################################
for i=1:NB_ITER
    @time iterate!()
    flush(stdout) # Forcing the printing
end
# println(TAB_R[2])