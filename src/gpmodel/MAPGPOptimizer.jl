import Base:isempty, empty!


mutable struct MAPGPOptimizer{O} <: GPOptimizer
    count::Int
    options::Dict{Symbol, Any}
    nlopt::O
end

function MAPGPOptimizer(;
    global_interval::Int=10,
    global_optimizer::NLOptimizer=ISRESOptimizer(;maxeval=10000),
    local_interval::Int=1,
    local_optimzier::NLOptimizer=LBFGSOptimizer(;maxeval=1000),
    domean = true, kern = true, noise = false, lik = false, 
    meanbounds = nothing, kernbounds = nothing,
    noisebounds = nothing, likbounds = nothing)

    domean && isnothing(meanbounds) && error("`meanbounds` must be given")
    kern && isnothing(kernbounds) && error("`kernbounds` must be given")
    noise && isnothing(noisebounds) && error("`noisebounds` must be given")
    lik && isnothing(likbounds) && error("`likbounds` must be given")

    model_options = Dict{Symbol,Any}(
        :domean=>domean, :kern=>kern, :noise=>noise, :lik=>lik,
        :meanbounds=>meanbounds, :kernbounds=>kernbounds, 
        :noisebounds=>noisebounds, :likbounds=>likbounds)
    nlopt = (
        (global_interval,global_optimizer), 
        (local_interval,local_optimzier)
    )
    MAPGPOptimizer(0, model_options, nlopt)
end

function optimizemodel!(o::MAPGPOptimizer, model::GP.GPBase, n::Int=1)
    for (itv, optimizer) in o.nlopt
        if ceil((o.count+n)/itv) > ceil(o.count/itv) 
            optimizemodel!(model, optimizer;o.options...)
        end
    end
    o.count += n
end

function check_bounds(name, x, lb, ub)
    d = ub-lb
    if any(abs.((x .- ub) ./ d) .< 1e-3)
        @warn "Found optimum `$name`(=$x) is close to the upper boundary (=$ub), You may try enlarging the `$name` bounds"
    elseif any(abs.((x .- lb) ./ d) .< 1e-3)
        @warn "Found optimum `$name`(=$x) is close to the lower boundary (=$lb), You may try enlarging the `$name` bounds"
    end
end

function optimizemodel!(gp::GP.GPBase, opt::NLOptimizer; 
    domean, kern, noise, lik, meanbounds, kernbounds, noisebounds, likbounds)
    params_kwargs = GP.get_params_kwargs(gp; domean=domean, kern=kern, noise=noise, lik=lik)
    f = if need_gradient(opt)
        (x, g) -> begin
            GP.set_params!(gp, x; params_kwargs...)
            GP.update_target_and_dtarget!(gp; params_kwargs...)
            @. g = gp.dtarget
            gp.target
        end
    else
        (x, g) -> begin
            GP.set_params!(gp, x; params_kwargs...)
            GP.update_target!(gp; params_kwargs...)
            gp.target
        end
    end
    lb, ub = GP.bounds(
        gp, noisebounds, meanbounds, kernbounds, likbounds; 
        domean=domean, kern=kern, noise=noise, lik=lik)
    fx, x, ret = maximize(opt, f, lb, ub, GP.get_params(gp; params_kwargs...))
    GP.set_params!(gp, x; params_kwargs...)
    GP.update_target_and_dtarget!(gp; params_kwargs...)
    domean && check_bounds("mean", GP.get_params(gp.mean), meanbounds...)
    kern && check_bounds("kernel", GP.get_params(gp.kernel), kernbounds...)
    noise && check_bounds("logNoise", GP.get_params(gp.logNoise), noisebounds...)
    lik && checkbounds("lik", GP.get_params(gp.lik), likbounds...)
    return fx, x, ret
end