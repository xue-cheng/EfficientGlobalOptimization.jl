import Base: ndims

abstract type TestFunction{N,M} end
ndims(::TestFunction{N,M}) where{N,M} = (N,M)
constrains(::TestFunction) = ()
include("single_obj_functions.jl")
include("multi_obj_functions.jl")