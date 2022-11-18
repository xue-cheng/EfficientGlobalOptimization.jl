
struct NoFeasibleInfill <: Exception
    msg::String
end

mutable struct Acquire
    global_opt::Optimizer
    local_opt::Optimizer
end

function Acquire(;
    global_optimizer=ISRESOptimizer(; maxeval=10000, ftol_rel=1e-2),
    local_optimizer=LBFGSOptimizer(; maxeval=1000),
)
    return Acquire(global_optimizer, local_optimizer)
end

function doacquire(opt::Acquire, obj, lb::AbstractVector, ub::AbstractVector)
    fx, x, ret = maximize(opt.global_opt, obj, :Auto, lb, ub, (lb + ub) / 2)
    fx, x, ret = maximize(opt.local_opt, obj, :Auto, lb, ub, x)
    return x, fx
end

function acquire(
    krg::Kriging,
    a::EGOAquisition,
    opt::Acquire,
    lb::AbstractVector,
    ub::AbstractVector,
    x_cache::AbstractVector...,
)
    f = x -> begin
        y = objective(a, krg, x)
        for xc in x_cache
            y *= 1.0 - exp(-sum((x - xc) .^ 2))
        end
        y
    end
    x, fx = doacquire(opt, f, lb, ub)
    if fx == 0
        throw(NoFeasibleInfill("Cannot find feasible infill. Kriging model is converged."))
    end
    return x, fx
end
