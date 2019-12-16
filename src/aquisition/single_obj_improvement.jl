abstract type ImprovementAquisition{S}<:SingleObjAquisition{S} end

function update_parameters!(a::ImprovementAquisition{Min}, krg::Kriging{N,1})where{N} 
    yo = krg.yscaler \ [minimum(krg.gps[1].y)]
    a.τ = yo[1]
end

function update_parameters!(a::ImprovementAquisition{Max}, krg::Kriging{N,1})where{N} 
    yo = krg.yscaler \ [maximum(krg.gps[1].y)]
    a.τ = yo[1]
end

improvement(a::ImprovementAquisition{Min}, y) = a.τ - y

improvement(a::ImprovementAquisition{Max}, y) = y - a.τ

mutable struct ProbabilityOfImprovement{S} <: ImprovementAquisition{S} 
    τ::Float64
end

ProbabilityOfImprovement(s::Symbol) = ProbabilityOfImprovement(sym2sense(s))
ProbabilityOfImprovement(s::EGOSense) = ProbabilityOfImprovement{s}(0.0)

@inline function (a::ProbabilityOfImprovement{Min})(μ, σ²)
    σ² == 0 && return float(μ < a.τ)
    Φ((a.τ - μ)/σ²)
end

@inline function (a::ProbabilityOfImprovement{Max})(μ, σ²)
    σ² == 0 && return float(μ > a.τ)
    Φ((μ - a.τ)/σ²)
end

mutable struct ExpectedImprovement{S} <: ImprovementAquisition{S} 
    τ::Float64
end

ExpectedImprovement(s::Symbol) = ExpectedImprovement(sym2sense(s))
ExpectedImprovement(s::EGOSense) = ExpectedImprovement{s}(0.0)

@inline function (a::ExpectedImprovement{Min})(μ, σ²)
    δ = a.τ - μ
    σ² == 0 && return max(δ, 0.)
    σ = √σ²
    n = δ/σ
    δ*Φ(n) + σ*φ(n)
end

@inline function (a::ExpectedImprovement{Max})(μ, σ²)
    δ = μ - a.τ
    σ² == 0 && return max(δ, 0.)
    σ = √σ²
    n = δ/σ
    δ*Φ(n) + σ*φ(n)
end
