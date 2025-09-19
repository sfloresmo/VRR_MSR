##################################################
# Printing the code's parameters
##################################################
println("PARALLEL | ",PARALLEL," | ORDER | ",ORDER," | COUPLING | ",COUPLING," | INIT | ",INIT," | ITER | ",ITER," | CUT | ",CUT)
println("LMAX | ",LMAX," | TC_OVER_DT | ",TC_OVER_DT," | TMAX_OVER_TC | ",TMAX_OVER_TC," | NB_TRG | ",NB_TRG," | NB_TRG2 | ",NB_TRG2)
println("TAB_GAMMA | ",Base.format_bytes(Base.summarysize(TAB_GAMMA))," | TAB_LAMBDA | ",Base.format_bytes(Base.summarysize(TAB_LAMBDA)))
println("TAB_NSTEPS | ",TAB_NSTEPS)
##################################################
flush(stdout) # forcing the print