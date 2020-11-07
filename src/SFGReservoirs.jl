module SFGReservoirs

using DifferentialEquations
using SFGAnalysis
using StaticArrays
using Optim

include("sfspectra.jl")
include("reservoirs.jl")

export  model,
        State

end # module
