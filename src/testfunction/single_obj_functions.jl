
abstract type SingleObjectiveTestFunction{N} <: TestFunction{N,1} end

struct Branin_Forrester <: SingleObjectiveTestFunction{2} end
function (f::Branin_Forrester)(x::AbstractVector)
    term1 = (x[2] - 5.1 / (4 * π^2) * x[1]^2 + 5 / π * x[1] - 6)^2
    term2 = 10 * (1 - 1 / (8π)) * cos(x[1])
    return term1 + term2 + 10 + 5x[1]
end

lowerbounds(::Branin_Forrester) = Float64[-5.0, 0.0]
upperbounds(::Branin_Forrester) = Float64[10.0, 15.0]
optimum(::Branin_Forrester) = -16.64402157084319, [-3.689285272296118, 13.629987729088747]

struct Ackley{N} <: SingleObjectiveTestFunction{N}
    a::Float64
    b::Float64
    c::Float64
    Ackley(n; a=20.0, b=0.2, c=2π) = new{n}(a, b, c)
end

function (f::Ackley{N})(x::AbstractVector) where {N}
    sum1 = 0.0
    sum2 = 0.0
    for i in 1:N
        sum1 += x[i]^2
        sum2 += cos(f.c * x[i])
    end
    return -f.a * exp(-f.b * sqrt(sum1 / N)) - exp(sum2 / N) + f.a + ℯ
end

lowerbounds(::Ackley{N}) where {N} = fill(-32.768, N)
upperbounds(::Ackley{N}) where {N} = fill(32.768, N)
optimum(::Ackley{N}) where {N} = 0.0, zeros(N)

struct Rosenbrock{N} <: SingleObjectiveTestFunction{N} end
Rosenbrock(n::Int) = Rosenbrock{n}()
function (f::Rosenbrock{N})(x::AbstractVector) where {N}
    sum = 0.0
    for i in 1:(N - 1)
        xᵢ = x[i]
        xⱼ = x[i + 1]
        sum += 100.0 * (xⱼ - xᵢ^2)^2 + (xᵢ - 1)^2
    end
    return sum
end

lowerbounds(::Rosenbrock{N}) where {N} = fill(-5.0, N)
upperbounds(::Rosenbrock{N}) where {N} = fill(15.0, N)
optimum(::Rosenbrock{N}) where {N} = 0.0, ones(N)

struct SumSphere{N} <: SingleObjectiveTestFunction{N} end

SumSphere(n::Integer) = SumSphere{n}()

(f::SumSphere{N})(x::AbstractVector) where {N} = sum(i -> i * x[i]^2, 1:N)
lowerbounds(::SumSphere{N}) where {N} = fill(-10.0, N)
upperbounds(::SumSphere{N}) where {N} = fill(10.0, N)
optimum(::SumSphere{N}) where {N} = 0.0, zeros(N)

abstract type ConstrainedSingleObjectTestFunction{N,C} <: SingleObjectiveTestFunction{N} end

constrains(f::ConstrainedSingleObjectTestFunction) = f.g

struct PressureVesselDesign <: ConstrainedSingleObjectTestFunction{4,4}
    g::NTuple{4,Function}
    function PressureVesselDesign()
        return new((
            x -> -x[1] + 0.0193 * x[3],
            x -> -x[2] + 0.00954 * x[3],
            x -> -pi * x[3]^2 * x[4] - (4 / 3) * pi * x[3]^3 + 1296000,
            x -> x[4] - 240,
        ))
    end
end

function (f::PressureVesselDesign)(x::AbstractVector)
    d1, d2, r, L = x
    return 0.6224 * d1 * r * L + 1.7781 * d2 * r^2 + 3.1661 * d1^2 * L + 19.84 * d1^2 * r
end

lowerbounds(::PressureVesselDesign) = [0.0625, 0.0625, 10.0, 10.0]
upperbounds(::PressureVesselDesign) = [99 * 0.0625, 99 * 0.0625, 200.0, 200.0]
function optimum(::PressureVesselDesign)
    f = 6059.714335048436
    x = [0.8125, 0.4375, 42.0984455958549, 176.6365958424394]
    return f, x
end

struct CantileverBeamDesign <: ConstrainedSingleObjectTestFunction{5,1}
    g::NTuple{1,Function}
    function PressureVesselDesign()
        return new((
            x -> 61 / x[1]^3 + 37 / x[2]^3 + 19 / x[3]^3 + 7 / x[4]^3 + 1 / x[5]^3 - 1,
        ))
    end
end

(f::CantileverBeamDesign)(x::AbstractVector) = 0.0624 * sum(x)

lowerbounds(::CantileverBeamDesign) = fill(0.01, 5)
upperbounds(::CantileverBeamDesign) = fill(100.0, 5)
function optimum(::CantileverBeamDesign)
    f = 1.339956367
    x = [6.0160159, 5.3091739, 4.4943296, 3.5014750, 2.15266533]
    return f, x
end
