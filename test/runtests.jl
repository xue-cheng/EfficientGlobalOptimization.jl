using Test, Statistics
import EfficientGlobalOptimization
EGO = EfficientGlobalOptimization
include("test_testfunction.jl")
include("test_optimizer.jl")
include("test_ego.jl")
include("test_pareto.jl")