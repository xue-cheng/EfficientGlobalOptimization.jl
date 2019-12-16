module EfficientGlobalOptimization


using KrigingModel

export KrigingTuner
export sampling, LHS, GRID
export EGO, optimum, acquire, history, lowerbounds, upperbounds
export Acquire, EGOAquisition, ProbabilityOfImprovement, ExpectedImprovement
#export ego_save
export ParetoSet

@enum EGOSense Min=1 Max=-1
const _sym2sense = Dict(Symbol(i)=>i for i in instances(EGOSense))
function sym2sense(s::Symbol)
    try
        return _sym2sense[s]
    catch
        error("invalid sense, must be :Min or :Max")
    end
end

include("sampling/sampling.jl")
include("aquisition/aquisition.jl")
include("acquire/acquire.jl")
include("krgtuner.jl")
include("EGO.jl")
#include("fileio.jl")
include("testfunction/testfunction.jl")
include("paretoset.jl")
end # module
