##################################################
# Macro to add parallelisation
# This macro is useful so that for PARALLEL = false,
# no threads are used. We can then peacefully check
# for un-wanted allocations
# @IMPROVE -- Here, we use the stdlib parallelisation
#             In some regimes, @batch from Polyester.jl
#             could be faster [most likely un-needed given the cost of the evaluations]
# @IMPROVE -- Not sure that the code for the macro is correct.
#             https://discourse.julialang.org/t/117808
##################################################
if PARALLEL # We asked for parallelisation
    macro parallel(ex::Expr)
        return esc(:( Threads.@threads $ex ))
    end
else        # No parallelisation
    macro parallel(ex::Expr)
        return esc(:( $(ex) ))
    end
end