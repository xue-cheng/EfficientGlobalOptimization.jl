import Base: append!, ndims, length, isempty, getindex
import Sobol

export EGO,optimum,acquire,acquire!,init_sampling

struct EGO{S}
    gp
    gp_optimizer
    aquisition
    acquire
    lowerbounds
    upperbounds
end

function EGO(aquisition::EGOAquisition{S},
    lowerbounds::AbstractVector,
    upperbounds::AbstractVector;
    mean::GP.Mean=GP.MeanConst(0.), 
    kernel::GP.Kernel=GP.SEArd(zeros(length(lowerbounds)), 5.),
    logNoise::Float64=-5.0,
    model_optimizer::GPOptimizer=MAPGPOptimizer(meanbounds = [-10., 10.], kernbounds = [fill(-10., length(lowerbounds)+1), fill(10., length(lowerbounds)+1)]),
    acquire::Acquire = Acquire()) where{S}

    d = length(lowerbounds)
    gp = GP.ElasticGPE(d; mean=mean, kernel=kernel, logNoise=logNoise)
    EGO{S}(gp, model_optimizer, aquisition, acquire, lowerbounds, upperbounds)
end

ndims(ego::EGO) = size(ego.gp.x, 1)
length(ego::EGO) = length(ego.gp.y)
isempty(ego::EGO) = isempty(ego.gp.y)

function init_sampling(ego::EGO, n::Int)
    d = length(ego.lowerbounds)
    x = Matrix{Float64}(undef, d, n)
    s = Sobol.SobolSeq(ego.lowerbounds,ego.upperbounds)
    for i = 1:n
        Sobol.next!(s, @view x[:,i])
    end
    x
end

function append!(ego::EGO, x::AbstractMatrix, y::AbstractArray)
    append!(ego.gp, x, y)
    update_parameters!(ego.aquisition, ego.gp)
    optimizemodel!(ego.gp_optimizer, ego.gp, length(y))
end

function append!(ego::EGO, x::AbstractVector, y::Float64)
    append!(ego.gp, reshape(x,:,1), [y])
    update_parameters!(ego.aquisition, ego.gp)
    optimizemodel!(ego.gp_optimizer, ego.gp, 1)
end

getindex(ego::EGO, i) = ego.gp.y[i], ego.gp.x[:,i]

function getindex(ego::EGO, symbol::Symbol, i)
    if symbol===:x
        ego.gp.x[:,i]
    elseif symbol===:y
        ego.gp.y[i]
    else
        throw(ArgumentError("invalid index: $symbol, must be `:x` or `:y`"))
    end
end

optimum(ego::EGO{Min}) = ego[argmin(ego.gp.y)]

optimum(ego::EGO{Max}) = ego[argmax(ego.gp.y)]

function acquire(ego::EGO, x_cache)
    acquire(ego.gp, ego.aquisition, ego.acquire_optimizer, ego.lowerbounds, ego.upperbounds, x_cache)
end
function acquire(ego::EGO)
    acquire(ego.gp, ego.aquisition, ego.acquire_optimizer, ego.lowerbounds, ego.upperbounds)
end