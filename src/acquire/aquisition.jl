import ElasticArrays:ElasticArray


abstract type EGOAquisition{S} end

update_parameters!(a::EGOAquisition, gp::GP.GPBase) = nothing

include("improvement.jl")
