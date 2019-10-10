import GaussianProcesses

const GP = GaussianProcesses

abstract type GPOptimizer end

optimizemodel!(o::GPOptimizer, model::GP.GPBase, n::Int=1) = nothing

include("MAPGPOptimizer.jl")

y_max(model::GP.GPBase) = isempty(model.y) ? -Inf : maximum(model.y)
y_min(model::GP.GPBase) = isempty(model.y) ?  Inf : minimum(model.y)


function predict_gp(gp::GP.GPBase, x::AbstractVector)
    m, v = GP.predict_f(gp, reshape(x,:,1))
    m[1], v[1]
end

predict_gp(gp::GP.GPBase, x::AbstractMatrix) = GP.predict_f(gp, x)