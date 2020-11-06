
function model(du, u, p::Tuple, t)
    FWHM = 15.0   # defined it as a global
    k, x, t0 = p
    t0 *= 1000 # input has to be smaller for optimization purposes
    f(t,t0) = gauss(t, t0, FWHM)

    # number of reservoirs
    N = length(du)

    # reset deltas to 0
    du .= 0.0

    # main loop
    for i in 1:N, j in 1:N
        i == j && continue
        du[i] -= k[i,j] * u[i]
        du[i] += k[j,i] * u[j]
        x[i,j] ≠ 0 && (du[i] -= u[i] * x[i,j] * f(t, t0)) # pumping
        x[j,i] ≠ 0 && (du[i] += u[j] * x[j,i] * f(t, t0)) # pumping don't flip j and i
    end

    nothing
end

function parameters_to_rate_matrix!(M, A, p)
    for i in eachindex(A)
        A[i] ≠ 0 && (M[i] = p[A[i]])
    end
    nothing
end

function parameters_to_rate_matrix!(M, a::Dict{Tuple{T,T}, N}, p) where {T, N}
    for key in keys(a)
        i = a[key] # parameter index
        M[key...] = p[i] # set the value
    end
    nothing
end

time_zero_from_parameters(p) = @view p[end]

num_states(a::Dict) = maximum(maximum.(keys(a)))
num_states(a::Array{T,2}) where T = size(a,1)

function gauss(t, t0, σ)
    exp(-(t - t0)^2 / σ^2)
end

struct State
    target::Int
    source::Array{Int,1}
end
State() = State(0, Int[])

function population_to_attenuation!(a, u, states)
    for i in eachindex(states)
        if states[i].target == 0
            a[i] = 1.0
            continue
        end
        a[i] = 0.0
        for source_state in states[i].source
            a[i] += u[source_state]
        end
        u_target = u[states[i].target]
        a[i] -= u_target
    end
end

function attenuated_sfspectrum!(y, x, a, A, ω, Γ)
    a .*= A
    sfspectrum!(y, x, a, ω, Γ)
end

function solution_to_attenuated_spectra!(Y, a, u, wn, dl, spec_params, sol, states)
    Threads.@threads for i in eachindex(dl)
        u .= sol[:,i]
        population_to_attenuation!(a, u, states)
        y = @view Y[:,i]
        attenuated_sfspectrum!(y, wn, a, spec_params[1], spec_params[2], spec_params[3])
    end
end
    