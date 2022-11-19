# File: EfficientGlobalOptimization.jl/src/aquisition/aquisition.jl
# Copyright (c) 2019-2022 XUE Cheng
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

abstract type EGOAquisition end

update_parameters!(::EGOAquisition, args...) = error("NOT IMPLEMENTED")
objective(::EGOAquisition, m::Kriging, x) = error("NOT IMPLEMENTED")

const φ = StatsFuns.normpdf
const Φ = StatsFuns.normcdf

include("single_obj_aquisitions.jl")
include("constrained_aquisition.jl")
