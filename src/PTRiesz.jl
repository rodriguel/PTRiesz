include("RWM.jl")
include("potential.jl")
include("utility.jl")
include("display.jl")
include("monotone.jl")

using JSON
using Dates
using ArgParse
#using SMTPClient

PATH_CLUSTER = "tests_cluster"
PATH_CLUSTER_OUTPUT = "tests_cluster/output"
PATH_LOCAL = "tests"
PATH_LOCAL_OUTPUT = "tests/output"

ARGS_SETTING = ArgParseSettings()

@add_arg_table! ARGS_SETTING begin
    "-s"
        arg_type = Float64
        help = "Homogeneity parameter"
    "-N"
        arg_type = Int
        help = "Number of particles"
    "-b"
        arg_type = Float64
        help = "Inverse temperature"
    "-n"
        arg_type = Int
        help = "Number of total samples"
        required = true
    "-c"
        arg_type = Int
        help = "Number of chains (default = number of available threads)"
        default = Threads.nthreads()
    "--brange"
        arg_type = Float64
        help = "Range of inverse temperatures (format 'begin step end')"
        nargs = 3
    "--blist"
        arg_type = Float64
        help = "List of (inverse) temperatures"
        nargs = '*'
    "--Nrange"
        arg_type = Int
        help = "Range of number of particles (format 'begin step end')"
        nargs = 3
    "--Nlist"
        arg_type = Int
        help = "List of number of particles"
        nargs = '*'
    "--slist"
        arg_type = Float64
        help = "List of homogenenity parameters"
        nargs = '*'
    "--srange"
        arg_type = Float64
        help = "Range of homogenenity parameters (format 'begin step end')"
        nargs = 3
    "--Trange"
        arg_type = Float64
        nargs = 3
    "--Tlist"
        arg_type = Float64
        nargs = '*'
        
end

function my_parse(s)
    args = parse_args(ARGS, s)
    # Parsing the homogeneity parameter(s)
    s = args["s"];  
    sr = args["srange"]; sl = args["slist"]
    if isempty(sr) && isempty(sl)
        ss = [s]
    else 
        ss = (isempty(sr)) ? sl : sr[1]:sr[2]:sr[3] 
    end 

    # Parsing the temperature(s)
    β = args["b"]
    br = args["brange"]; bl = args["blist"]
    βs = (isempty(br)) ? bl : br[1]:br[2]:br[3]

    Tr = args["Trange"]; Tl = args["Tlist"]
    Tr = (isempty(Tr)) ? Tl : Tr[1]:Tr[2]:Tr[3]
    if isempty(βs)
        βs = 1 ./ Tr
    end

    # Parsing number(s) of particles
    N = args["N"]
    Nr = args["Nrange"]; Nl = args["Nlist"]
    Ns = (isempty(Nr)) ? Nl : Nr[1]:Nr[2]:Nr[3] 
    # Number of chains and total sample
    nchains = args["c"]; nsamples = args["n"]
    return N, Array(Ns), s, Array(ss), β, Array(βs), nsamples, nchains
end 
