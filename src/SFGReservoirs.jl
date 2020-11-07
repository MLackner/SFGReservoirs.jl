module SFGReservoirs

using DifferentialEquations
using SFGAnalysis
using StaticArrays

include("sfspectra.jl")
include("reservoirs.jl")

export  model,
        State

end # module
