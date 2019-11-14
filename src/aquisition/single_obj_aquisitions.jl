

abstract type EGOAquisition{S} end

update_parameters!(a::EGOAquisition, gp::GP.GPBase) = nothing

include("single_obj_improvement.jl")
