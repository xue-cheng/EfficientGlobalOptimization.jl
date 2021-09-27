abstract type ImprovementAquisition{S} <: SingleObjAquisition{S} end

function update_parameters!(a::ImprovementAquisition{Min}, krg::Kriging{N,1}) where {N}
    yo = krg.yscaler \ [minimum(krg.gps[1].y)]
    return a.τ = yo[1]
end

function update_parameters!(a::ImprovementAquisition{Max}, krg::Kriging{N,1}) where {N}
    yo = krg.yscaler \ [maximum(krg.gps[1].y)]
    return a.τ = yo[1]
end

improvement(a::ImprovementAquisition{Min}, y) = a.τ - y

improvement(a::ImprovementAquisition{Max}, y) = y - a.τ

mutable struct ProbabilityOfImprovement{S} <: ImprovementAquisition{S}
    τ::Float64
end

ProbabilityOfImprovement(s::Symbol) = ProbabilityOfImprovement(sym2sense(s))
ProbabilityOfImprovement(s::EGOSense) = ProbabilityOfImprovement{s}(0.0)

function objective(a::ProbabilityOfImprovement{Min}, model::Kriging{N,1}, x) where {N}
    m, v = predict_full(model, x)
    μ = m[1]
    σ² = v[1]
    if σ² == 0
        float(μ < a.τ)
    else
        Φ((a.τ - μ) / σ²)
    end
end

function objective(a::ProbabilityOfImprovement{Max}, model::Kriging{N,1}, x) where {N}
    m, v = predict_full(model, x)
    μ = m[1]
    σ² = v[1]
    if σ² == 0
        float(μ > a.τ)
    else
        Φ((μ - a.τ) / σ²)
    end
end
mutable struct ExpectedImprovement{S} <: ImprovementAquisition{S}
    τ::Float64
end

ExpectedImprovement(s::Symbol) = ExpectedImprovement(sym2sense(s))
ExpectedImprovement(s::EGOSense) = ExpectedImprovement{s}(0.0)

function objective(a::ExpectedImprovement{Min}, model::Kriging{N,1}, x) where {N}
    m, v = predict_full(model, x)
    μ = m[1]
    σ² = v[1]
    δ = a.τ - μ
    if σ² == 0
        max(δ, 0.0)
    else
        σ = √σ²
        n = δ / σ
        δ * Φ(n) + σ * φ(n)
    end
end

function objective(a::ExpectedImprovement{Max}, model::Kriging{N,1}, x) where {N}
    m, v = predict_full(model, x)
    μ = m[1]
    σ² = v[1]
    δ = μ - a.τ
    if σ² == 0
        max(δ, 0.0)
    else
        σ = √σ²
        n = δ / σ
        δ * Φ(n) + σ * φ(n)
    end
end
