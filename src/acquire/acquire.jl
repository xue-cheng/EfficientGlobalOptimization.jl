
mutable struct Acquire
    global_opt::KrigingModel.Optimizer
    local_opt::KrigingModel.Optimizer
end

function Acquire(; global_optimizer = ISRESOptimizer(maxeval = 10000, ftol_rel = 1e-2),
                   local_optimizer = LBFGSOptimizer(maxeval = 1000))
    Acquire(global_optimizer, local_optimizer)
end

function doacquire(opt::Acquire, obj, lb::AbstractVector, ub::AbstractVector)
    fx, x, ret = maximize(opt.global_opt, obj, :Auto, lb, ub, (lb + ub) / 2)
    fx, x, ret = maximize(opt.local_opt,  obj, :Auto, lb, ub, x)
    return x, fx
end

function acquire(krg::Kriging{N,1}, 
    a::SingleObjAquisition,
    opt::Acquire,
    lb::AbstractVector,
    ub::AbstractVector,
    x_cache::AbstractVector...) where {N}
    f = x->begin 
        m, v = predict_full(krg, x)
        y = a(m[1], v[1])
        for xc in x_cache
            y *= 1.0 - exp(-sum((x - xc).^2))
        end
        y
    end
    doacquire(opt, f, lb, ub)
end
