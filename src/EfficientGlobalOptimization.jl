module EfficientGlobalOptimization


export EGO,optimum,acquire,init_sampling
export Acquire, EGOAquisition, ProbabilityOfImprovement, ExpectedImprovement
export MAPGPOptimizer
export wrap_function
export LBFGSOptimizer, CRSOptimizer, ISRESOptimizer

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
include("acquire/acquire.jl")
include("EGO.jl")

include("testfunction/testfunction.jl")

end # module
