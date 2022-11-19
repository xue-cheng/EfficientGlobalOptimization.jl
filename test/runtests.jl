# File: EfficientGlobalOptimization.jl/test/runtests.jl
# Copyright (c) 2019-2022 XUE Cheng
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Test, Statistics
using EfficientGlobalOptimization
using GaussianProcesses
using MAT

include("test_testfunction.jl")
include("test_scaler.jl")
include("test_optim.jl")
include("test_kriging.jl")
include("test_ego.jl")
include("test_cons.jl")
include("test_pareto.jl")