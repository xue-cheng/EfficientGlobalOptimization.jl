import Statistics:mean

abstract type MultiObjectiveTestFunction{N, M} <: TestFunction{N,M} end

"""
    -A≤x≤A. 
Values of `A` from 10 to 10⁵ have been used successfully. 
Higher values of `A` increase the difficulty of the problem.

# References
[1] Schaffer, J. David (1984). Multiple Objective Optimization with Vector Evaluated Genetic Algorithms. Proceedings of the First Int. Conference on Genetic Algortihms, Ed. G.J.E Grefensette, J.J. Lawrence Erlbraum (PhD). Vanderbilt University. OCLC 20004572
[2] https://en.wikipedia.org/wiki/Test_functions_for_optimization#Test_functions_for_multi-objective_optimization
"""
struct Schaffer <: MultiObjectiveTestFunction{1,2} end

(f::Schaffer)(x::AbstractVector) = [x[1]^2, (x[1]-2)^2]

lowerbounds(::Schaffer) = Float64[-10.]
upperbounds(::Schaffer) = Float64[10.]

"""
# References
[1] Zitzler, E., Deb, K., & Thiele, L. (2000). Comparison of multiobjective evolutionary algorithms: Empirical results. Evolutionary Computation, 8(2), 173–195.
"""
struct ZDT1{N} <: MultiObjectiveTestFunction{N,2} end
function (f::ZDT1{N})(x::AbstractVector) where {N}
    f = Vector{Float64}(undef,2)
    f[1] = x[1]
    g = 1.0 + 9.0*mean(x[2:N])
    h = 1.0 - sqrt(f[1]/g)
    f[2] = g*h
    f
end

lowerbounds(::ZDT1{N}) where {N} = zeros(N)
upperbounds(::ZDT1{N}) where {N} = ones(N)

"""
# References
[1] Zitzler, E., Deb, K., & Thiele, L. (2000). Comparison of multiobjective evolutionary algorithms: Empirical results. Evolutionary Computation, 8(2), 173–195.
"""
struct ZDT2{N} <: MultiObjectiveTestFunction{N,2} end

function (f::ZDT2{N})(x::AbstractVector) where {N}
    f = Vector{Float64}(undef,2)
    f[1] = x[1]
    g = 1.0 + 9.0*mean(x[2:N])
    h = 1.0 - (f[1]/g)^2
    f[2] = g*h
    f
end

lowerbounds(::ZDT2{N}) where {N} = zeros(N)
upperbounds(::ZDT2{N}) where {N} = ones(N)

"""
# References
[1] Zitzler, E., Deb, K., & Thiele, L. (2000). Comparison of multiobjective evolutionary algorithms: Empirical results. Evolutionary Computation, 8(2), 173–195.
"""
struct ZDT3{N} <: MultiObjectiveTestFunction{N,2} end

function (f::ZDT3{N})(x::AbstractVector) where {N}
    f = Vector{Float64}(undef,2)
    f[1] = x[1]
    g = 1.0 + 9.0*mean(x[2:N])
    h = 1.0 - sqrt(f[1]/g)-(f[1]/g)*sin(10π*f[1])
    f[2] = g*h
    f
end
lowerbounds(::ZDT3{N}) where {N} = zeros(N)
upperbounds(::ZDT3{N}) where {N} = ones(N)

function ZDT(prob_id::Int, param::Int)
    if param < 2
        error("invalid `param`=$param, must be greater than 1")
    end
    if prob_id == 1
        ZDT1{param}()
    elseif prob_id == 2
        ZDT2{param}()
    elseif prob_id == 3
        ZDT3{param}()
    else
        error("invalid `prob_id`=$prob_id, must be one of (1, 2, 3)")
    end
end

