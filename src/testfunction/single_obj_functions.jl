
abstract type SingleObjectiveTestFunction{N}<:TestFunction{N,1} end

struct Branin_Forrester <: SingleObjectiveTestFunction{2} end
function (f::Branin_Forrester)(x::AbstractVector) 
    term1 = (x[2]-5.1/(4*π^2)*x[1]^2+5/π*x[1]-6)^2
    term2 =  10*(1-1/(8π))*cos(x[1])
    term1 + term2 + 10 + 5x[1]
end

lowerbounds(::Branin_Forrester) = Float64[-5.0,  0.0]
upperbounds(::Branin_Forrester) = Float64[10.0, 15.0]
optimum(::Branin_Forrester)= -16.64402157084319,  [-3.689285272296118,13.629987729088747]


struct Ackley{N} <: SingleObjectiveTestFunction{N}
    a::Float64
    b::Float64
    c::Float64
    Ackley(n;a=20.0, b=0.2, c=2π) = new{n}(a,b,c)
end

function (f::Ackley{N})(x::AbstractVector) where {N}
    sum1 = 0.0
    sum2 = 0.0
    for i = 1:N
        sum1 += x[i]^2
        sum2 += cos(f.c*x[i])
    end
    -f.a*exp(-f.b*sqrt(sum1/N)) - exp(sum2/N) + f.a + ℯ
end

lowerbounds(::Ackley{N}) where{N} = fill(-32.768, N)
upperbounds(::Ackley{N}) where{N} = fill( 32.768, N)
optimum(::Ackley{N}) where{N} = 0.0, zeros(N)

struct Rosenbrock{N} <: SingleObjectiveTestFunction{N} end
Rosenbrock(n::Int) = Rosenbrock{n}()
function (f::Rosenbrock{N})(x::AbstractVector) where {N}
    sum = 0.0
    for i = 1:N-1
        xᵢ = x[i]
        xⱼ = x[i+1]
        sum += 100.0*(xⱼ-xᵢ^2)^2+(xᵢ-1)^2
    end
    sum
end

lowerbounds(::Rosenbrock{N}) where{N} = fill(-5.0, N)
upperbounds(::Rosenbrock{N}) where{N} = fill(15.0, N)
optimum(::Rosenbrock{N}) where{N} = 0.0, ones(N)


struct SumSphere{N} <: SingleObjectiveTestFunction{N} end

SumSphere(n::Integer) = SumSphere{n}()

(f::SumSphere{N})(x::AbstractVector) where {N} = sum(i->i*x[i]^2,1:N)
lowerbounds(::SumSphere{N}) where{N} = fill(-10.0, N)
upperbounds(::SumSphere{N}) where{N} = fill( 10.0, N)
optimum(::SumSphere{N}) where{N} = 0.0, zeros(N)