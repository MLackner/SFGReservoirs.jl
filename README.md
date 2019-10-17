# SFGReservoirs

[![Build Status](https://travis-ci.com/MLackner/SFGReservoirs.jl.svg?branch=master)](https://travis-ci.com/MLackner/SFGReservoirs.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/MLackner/SFGReservoirs.jl?svg=true)](https://ci.appveyor.com/project/MLackner/SFGReservoirs-jl)
[![Codecov](https://codecov.io/gh/MLackner/SFGReservoirs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MLackner/SFGReservoirs.jl)
[![Coveralls](https://coveralls.io/repos/github/MLackner/SFGReservoirs.jl/badge.svg?branch=master)](https://coveralls.io/github/MLackner/SFGReservoirs.jl?branch=master)

## Installation

In the REPL package mode type:
```julia
pkg> add https://github.com/MLackner/SFGReservoirs.jl#master
```

You will also need the `DifferentialEquations` package.
```julia
pkg> add DifferentialEquations
```

## Usage

Load packages:
```julia
using DifferentialEquations
using SparseArrays
```

Define the number of states in the model. We initialize the model with everything
in the ground state.
```julia
N = 4      # total number of states
u0 = zeros(N)
u0[2] = 1.0
# u0[1]:   heat bath
# u0[2]:   ground state
# u0[3]:   pumped state
# u0[4:N]: additional states
```

Set the rate constants `k[i,j]` between the states `i` and `j`. This is done so
that `k[i,j]` denotes the rate constant in 1/ps from state `i` to state `j`.
```julia
k = spzeros(N,N)
k[2,3] = 0.01     # ground state --> pumped state
k[3,4] = 0.05     # pumped state --> state 4
k[4,3] = 0.01     # state 4      --> pumped state
k[3,1] = 0.001    # pumped state --> heat bath
k[4,1] = 0.001    # state 4      --> heat bath
k[1,2] = 0.0001   # heat bath    --> ground state
```

Set time span and solve the model:
```julia
tstart = -30.0
tend = 270.0
tspan = (tstart, tend)

prob = ODEProblem(model, u0, tspan, k, maxiter=1)
sol = solve(prob)
```
