##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--parallel"
    help = "Parallelisation: true/false"
    arg_type = Bool
    default = true
    "--lmax"
    help = "Maximum harmonic number"
    arg_type = Int64
    default = 4
    "--coupling"
    help = "Coupling considered: Quad"
    arg_type = String
    default = "Quad"
    "--Tc_over_dt"
    help = "Discretisation in time wrt Tc"
    arg_type = Int64
    default = 1
    "--Tmax_over_Tc"
    help = "Overall duration wrt Tc"
    arg_type = Int64
    default = 2
    "--order"
    help = "Order of the approximation: 0/1/2"
    arg_type = Int64
    default = 2
    "--init"
    help = "Initialisation used: Bare/Gaussian"
    arg_type = String
    default = "Bare"
    "--iter_init"
    help = "Iterations to read from. If 0, we initialise with INIT"
    arg_type = Int64
    default = 0
    "--nb_iter"
    help = "Number of iterations performed"
    arg_type = Int64
    default = 10
    "--cut"
    help = "Whether to cut in time: true/false"
    arg_type = Bool
    default = false
end
##################################################
parsed_args = parse_args(tabargs)
##################################################
# General parameters
##################################################
const PARALLEL = parsed_args["parallel"] # Boolean to determine if we parallelise or not
const LMAX = parsed_args["lmax"] # Maximum harmonic number
const COUPLING = parsed_args["coupling"] # Coupling considered
const TC_OVER_DT = parsed_args["Tc_over_dt"] # Discretisation step
const TMAX_OVER_TC = parsed_args["Tmax_over_Tc"] # Overall duration of the time discretisation
const ORDER = parsed_args["order"] # Order of the approximation
const INIT = parsed_args["init"] # Initialisation used
const ITER_INIT = parsed_args["iter_init"] # Initial iterations. If > 0, we read from the files
const NB_ITER = parsed_args["nb_iter"] # Number of iterations performed
const CUT = parsed_args["cut"] # Whether to cut in time
##################################################
# Checking the sanity of some of the parameters
##################################################
if (COUPLING != "Quad")
    error("Supported COUPLING: Quad")
end
#####
if  (ORDER != 0) &&
    (ORDER != 1) &&
    (ORDER != 2)
    error("Supported ORDER: 0/1/2")
end
#####
if  (INIT != "Bare") &&
    (INIT != "Gaussian")
    error("Supported INIT: Bare/Gaussian")
end
#####
if (ITER_INIT < 0)
    error("Supported ITER_INIT: must be >= 0")
end