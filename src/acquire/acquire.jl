
mutable struct Acquire
    global_opt::NLOptimizer
    local_opt::NLOptimizer
end

function Acquire(; global_optimizer=ISRESOptimizer(maxeval=10000, ftol_rel=1e-2),
                            local_optimizer=LBFGSOptimizer(maxeval=100))
    Acquire(global_optimizer, local_optimizer)
end

function doacquire(opt::Acquire, obj, lb::AbstractVector, ub::AbstractVector)
    fx, x, ret = maximize(opt.global_opt, obj, :Auto, lb, ub, lb)
    fx, x, ret = maximize(opt.local_opt,  obj, :Auto, lb, ub, x)
    return x
end

function acquire(gp::GP.GPBase, 
    a::EGOAquisition,
    opt::Acquire,
    lb::AbstractVector,
    ub::AbstractVector,
    x_cache::AbstractVector{V}) where {V<:AbstractVector}

    f = x -> begin 
        y = a(predict_gp(gp, x)...)
        for xc in x_cache
            y *= 1.0 - exp(-sum((x-xc).^2))
        end
        y
    end
    doacquire(opt, f, lb,ub)
end

function acquire(gp::GP.GPBase, 
    a::EGOAquisition,
    opt::Acquire,
    lb::AbstractVector,
    ub::AbstractVector,
    x_cache::AbstractMatrix)
    f = x -> begin 
        y = a(predict_gp(gp, x)...)
        for xc in eachcol(x_cache)
            y *= 1.0 - exp(-sum((x-xc).^2))
        end
        y
    end
    doacquire(opt, f, lb,ub)
end


function acquire(gp::GP.GPBase, 
    a::EGOAquisition,
    opt::Acquire,
    lb::AbstractVector,
    ub::AbstractVector)
    doacquire(opt, (x)->a(predict_gp(gp, x)...), lb,ub)
end