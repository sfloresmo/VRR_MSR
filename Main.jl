##################################################
include("Packages.jl") # Loading the packages
include("Args.jl") # Parsing the command-line
include("Parallel.jl") # Macro to deal with parallelisation
include("Constants.jl") # Physical constants
include("Utils.jl") # Utility functions
include("Elsasser.jl") # Elsasser coefficients
include("W2.jl") # To compute the diagram W2
include("Coupling.jl") # Coupling coefficients
include("Nsteps.jl") # Depth in time
include("R_bare.jl") # Bare response function
include("R_Gaussian.jl") # Gaussian response function
include("R.jl") # Dressed response function
include("G.jl") # Dressed propagator
include("Triangles.jl") # Triangular triplets of harmonics
include("Gamma_bare.jl") # Bare vertex
include("Theta.jl") # Second order contribution to Gamma
include("G_3.jl") # 3-point correlation function
include("Tensor.jl") # Dressed vertex and Theta
include("Triangles2.jl") # Triangular triplets of type 2
include("Tensor2.jl") # Wrapper for the computation of Lambda
include("Lambda.jl") # Lambda
include("Sigma.jl") # Self-energy
include("Ranges_time.jl") # Ranges in time
include("IO.jl") # Dumping/Reading files
include("Init.jl") # Initialising R and Gamma
include("Inversion.jl") # To inverse Toeplitz matrices
include("Iteration.jl") # To iterate and search for fix-points
include("Printing.jl") # Printing the parameters of the code
##################################################