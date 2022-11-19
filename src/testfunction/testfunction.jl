# File: EfficientGlobalOptimization.jl/src/testfunction/testfunction.jl
# Copyright (c) 2019-2022 XUE Cheng
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

abstract type TestFunction{N,M} end

ndims(::TestFunction{N,M}) where{N,M} = (N,M)
constrains(::TestFunction) = ()

include("single_obj_functions.jl")
include("multi_obj_functions.jl")