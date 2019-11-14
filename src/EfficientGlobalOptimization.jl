module EfficientGlobalOptimization


export EGO, optimum, acquire, init_sampling, history, lowerbounds, upperbounds
export Acquire, EGOAquisition, ProbabilityOfImprovement, ExpectedImprovement
export MAPGPOptimizer
export wrap_function
export LBFGSOptimizer, CRSOptimizer, ISRESOptimizer
export ego_save
export ParetoSet

@enum EGOSense Min Max
const _sym2sense = Dict(Symbol(i)=>i for i in instances(EGOSense))
function sym2sense(s::Symbol)
    try
        return _sym2sense[s]
    catch
        error("invalid sense, must be :Min or :Max")
    end
end

include("optimizer/optimizer.jl")
include("gpmodel/gpmodel.jl")
include("aquisition/aquisition.jl")
include("acquire/acquire.jl")
include("EGO.jl")
include("fileio.jl")
include("testfunction/testfunction.jl")
include("paretoset.jl")
end # module
