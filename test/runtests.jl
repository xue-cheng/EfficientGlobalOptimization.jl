using Test, Statistics
using EfficientGlobalOptimization
using GaussianProcesses

include("test_testfunction.jl")
include("test_scaler.jl")
include("test_optim.jl")
include("test_kriging.jl")
include("test_ego.jl")
include("test_cons.jl")
include("test_pareto.jl")