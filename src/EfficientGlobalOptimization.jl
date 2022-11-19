# File: EfficientGlobalOptimization.jl/src/EfficientGlobalOptimization.jl
# Copyright (c) 2019-2022 XUE Cheng
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module EfficientGlobalOptimization

using Statistics: mean!, mean, std
using ElasticArrays: ElasticArray
using Sobol
using StatsFuns
using LinearAlgebra
using NLopt, ForwardDiff, DiffResults

import GaussianProcesses
import Base: append!, ndims, length, isempty, getindex, setproperty!, getproperty

const GP = GaussianProcesses

export StandardScaler, MinMaxScaler
export GlobalOptimizer, LocalOptimizer, StochasticSearchOptimizer
export CRSOptimizer, ISRESOptimizer, LBFGSOptimizer, BOBYQAOptimizer
export minimize, maximize
export Kriging, tune!, predict_full, getx, gety, get_samples
export KrigingTuner
export sampling, LHS, GRID
export NoFeasibleInfill
export EGO, optimum, acquire, history, lowerbounds, upperbounds
export Acquire, EGOAquisition, ProbabilityOfImprovement, ExpectedImprovement
export constrains, ConstrainedAquisition
#export ego_save
export ParetoSet

@enum EGOSense Min = 1 Max = -1
const _sym2sense = Dict(Symbol(i) => i for i in instances(EGOSense))
function sym2sense(s::Symbol)
    try
        return _sym2sense[s]
    catch
        error("invalid sense, must be :Min or :Max")
    end
end

include("scaler/scaler.jl")
include("optimizer/optimizer.jl")
include("kriging/kriging.jl")
include("sampling/sampling.jl")
include("aquisition/aquisition.jl")
include("acquire/acquire.jl")
include("krgtuner.jl")
include("EGO.jl")
include("testfunction/testfunction.jl")
include("paretoset.jl")
end # module
