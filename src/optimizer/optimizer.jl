import NLopt, ForwardDiff, DiffResults
import Base: setproperty!, getproperty

include("utils.jl")

abstract type NLOptimizer end

minimize(o::NLOptimizer, f, g, lb, ub, x) = minimize(o, wrap_function(f,g), lb, ub, x)
maximize(o::NLOptimizer, f, g, lb, ub, x) = maximize(o, wrap_function(f,g), lb, ub, x)
isglobal(o::NLOptimizer) = error("not implemented")
need_gradient(o::NLOptimizer) = error("not implemented")


include("stochastic_opt.jl")
include("gradient_opt.jl")
include("gradfree_opt.jl")