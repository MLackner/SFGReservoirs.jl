
function model(du, u, p::Tuple, t)
    k, x, t0, σ = p
    t0 *= 1000 # input has to be smaller for optimization purposes
    pump_strength = gauss(t, t0, σ)

    # number of reservoirs
    N = length(du)

    # reset deltas to 0
    du .= 0.0

    # main loop
    for i in 1:N, j in 1:N
        i == j && continue
        du[i] -= k[i,j] * u[i]
        du[i] += k[j,i] * u[j]
        x[i,j] ≠ 0 && (du[i] -= u[i] * x[i,j] * pump_strength) # pumping
        x[j,i] ≠ 0 && (du[i] += u[j] * x[j,i] * pump_strength) # pumping don't flip j and i
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

function solution_to_attenuated_spectra!(Y, a, wn, dl, spec_params, sol, states)
    Threads.@threads for i in eachindex(dl)
        _a = @view a[:,i]
        _u = @view sol[:,i]
        _y = @view Y[:,i]
        population_to_attenuation!(_a, _u, states)
        attenuated_sfspectrum!(_y, wn, _a, spec_params[1], spec_params[2], spec_params[3])
    end
end

function fit(data, spec_params, u0, p_initial, assignments, state_array)
    function loss(p)
            SFGReservoirs.parameters_to_rate_matrix!(k, ak, p)
            SFGReservoirs.parameters_to_rate_matrix!(x, ax, p)

            # these should be at fixed position
            t0 = p[end-1]
            σ  = p[end]

            tspan = extrema(delaytimes)
            params = (k,x,t0,σ)
            prob = ODEProblem(model, u0, tspan, params)
            sol = solve(prob; saveat=delaytimes)

            SFGReservoirs.solution_to_attenuated_spectra!(
                Y, a, wavenumbers, delaytimes, spec_params, sol, state_array
            )

            @. R = Y - Y′
            sum(abs2.(R))
    end

    # unpack data
    delaytimes, wavenumbers, Y′ = data
    ak, ax = assignments

    n_states = length(u0)
    # preallocations
    k  = @MMatrix zeros(n_states, n_states)
    x  = @MMatrix zeros(n_states, n_states)
    # caclulated values
    Y = zeros(ComplexF64, (length(wavenumbers), length(delaytimes)))
    # residuals
    R = zeros(ComplexF64, (length(wavenumbers), length(delaytimes)))
    # attenuation factors
    a = zeros(length(state_array), length(delaytimes))

    optimize(loss, p_initial, NelderMead(),
            Optim.Options(
                time_limit=120,
            );
    )
end