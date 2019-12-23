import Base: append!, ndims, length, isempty, getindex
struct NoFeasibleInfill <: Exception
    msg::String
end

struct EGO{S,KRG <: Kriging,AQ <: SingleObjAquisition,AC <: Acquire,LB,UB}
    krg::KRG
    tuner::KrigingTuner
    aquisition::AQ
    acquire::AC
    lowerbounds::LB
    upperbounds::UB
    function EGO(kriging::Kriging{N,1},
        tuner::KrigingTuner,
        lowerbounds::AbstractVector,
        upperbounds::AbstractVector;
        acquire::Acquire = Acquire(),
        aquisition::SingleObjAquisition{S} = ExpectedImprovement(:Min)) where {N,S}
        length(lowerbounds) != N && throw(ArgumentError("length of `lowerbounds` must be $N"))
        length(upperbounds) != N && throw(ArgumentError("length of `upperbounds` must be $N"))
        lob = copy(lowerbounds)
        upb = copy(upperbounds)
        new{S,typeof(kriging),typeof(aquisition),typeof(acquire),typeof(lob),typeof(upb)}(kriging,
            tuner,
            aquisition,
            acquire,
            lob,
            upb)
    end
end

Base.length(ego::EGO) = length(ego.krg.gps[1].y)
Base.isempty(ego::EGO) = isempty(ego.krg.gps[1].y)
lowerbounds(ego::EGO) = ego.lowerbounds
upperbounds(ego::EGO) = ego.upperbounds

function Base.append!(ego::EGO, x::AbstractMatrix, y::AbstractVector)
    append!(ego.krg, x, transpose(y))
    update_parameters!(ego.aquisition, ego.krg)
    tunekrg!(ego.krg, ego.tuner)
end

Base.append!(ego::EGO, x::AbstractVector, y::Float64) = append!(ego, reshape(x, :, 1), [y])

optimum(ego::EGO{Min}) = ego.krg[argmin(ego.krg.gps[1].y)]

optimum(ego::EGO{Max}) = ego.krg[argmax(ego.krg.gps[1].y)]

function acquire(ego::EGO, x_cache::AbstractVector...; retry::Integer = 5, verbose::Bool=false)
    x, fx = acquire(ego.krg, ego.aquisition, ego.acquire, ego.lowerbounds, ego.upperbounds, x_cache...)
    if fx==0 
        if retry > 1
            if verbose
                @info "Cannot find feasible infill. retry($(retry-1))"
            end
            return acquire(ego, x_cache...; retry=retry-1, verbose=verbose)
        else
            throw(NoFeasibleInfill("Cannot find feasible infill. Kriging model is converged."))
        end
    end
    if verbose
        @info "Infill acquired, x=$x, fx=$fx"
    end
    x, fx
end

history(ego::EGO) = get_samples(ego.krg)