import GaussianProcesses: φ, Φ
export ProbabilityOfImprovement, ExpectedImprovement
abstract type ImprovementAquisition{S}<:EGOAquisition{S} end

update_parameters!(a::ImprovementAquisition{Min}, gp::GP.GPBase) = a.τ = y_min(gp)
update_parameters!(a::ImprovementAquisition{Max}, gp::GP.GPBase) = a.τ = y_max(gp)

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
