abstract type EGOAquisition end

update_parameters!(::EGOAquisition, args...) = error("NOT IMPLEMENTED")
objective(::EGOAquisition, m::Kriging, x) = error("NOT IMPLEMENTED")

import StatsFuns

const φ = StatsFuns.normpdf
const Φ = StatsFuns.normcdf

include("single_obj_aquisitions.jl")
include("constrained_aquisition.jl")
