# File: EfficientGlobalOptimization.jl/src/aquisition/single_obj_aquisitions.jl
# Copyright (c) 2019-2022 XUE Cheng
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

abstract type SingleObjAquisition{S}<:EGOAquisition end

include("single_obj_improvement.jl")
