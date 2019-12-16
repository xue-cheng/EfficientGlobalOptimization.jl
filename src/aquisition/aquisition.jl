abstract type EGOAquisition end

update_parameters!(a::EGOAquisition, args...) = error("NOT IMPLEMENTED")

import GaussianProcesses: φ, Φ
include("single_obj_aquisitions.jl")
#include("multi_obj_aquisitions.jl")
