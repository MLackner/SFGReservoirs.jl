using SFGReservoirs
using SFGAnalysis
using Test
using StaticArrays
using DifferentialEquations
using BenchmarkTools


begin
      N = 4 # number of states
      p = [0.1, 0.2, 30.0, -1.0] # [k1, k2, x, t0]
      k = @MMatrix zeros(N,N)
      assignments_k = Dict(
            (1,2) => 1, # the entry k[1,2] should map to p[1]
            (4,1) => 2,
            (4,2) => 2,
      )
      assignments_k_matrix = @MMatrix zeros(Int, 4,4)
      assignments_k_matrix[1,2] = 1
      assignments_k_matrix[4,1] = 2
      assignments_k_matrix[4,2] = 2

      @btime SFGReservoirs.parameters_to_rate_matrix!(k, assignments_k, p)
      @btime SFGReservoirs.parameters_to_rate_matrix!(k, assignments_k_matrix, p)
end

@testset "Model Setup" begin
      N = 4 # number of states
      p = [0.1, 0.2, 30.0, -1.0] # [k1, k2, x, t0]
      k = @MMatrix zeros(N,N)
      assignments_k = Dict(
            (1,2) => 1, # the entry k[1,2] should map to p[1]
            (4,1) => 2,
            (4,2) => 2,
      )
      SFGReservoirs.parameters_to_rate_matrix!(k, assignments_k, p)
      @test k[1,2] == p[1]
      @test k[4,1] == p[2]
      @test k[4,2] == p[2]
      @test size(k) == (4,4)

      k2 = @MMatrix zeros(N,N)
      assignments_k_matrix = @MMatrix zeros(Int, 4,4)
      assignments_k_matrix[1,2] = 1
      assignments_k_matrix[4,1] = 2
      assignments_k_matrix[4,2] = 2
      SFGReservoirs.parameters_to_rate_matrix!(k2, assignments_k_matrix, p)
      @test k2 == k


      x = @MMatrix zeros(N,N)
      assignments_x = Dict(
            (1,2) => 3,
            (2,1) => 3,
      )
      SFGReservoirs.parameters_to_rate_matrix!(x, assignments_x, p)
end

@test SFGReservoirs.gauss(0, 0, 3) == 1
@test SFGReservoirs.gauss(100, 0, 3) ≈ 0
@test SFGReservoirs.gauss(4, 1, 3) ≈ 1/ℯ

begin
      du = @MArray zeros(3)
      u  = @MArray zeros(3)
      k  = @MMatrix rand(3,3)
      x  = @MMatrix zeros(3,3)
      t0 = 0.0
      x[1,2] = 0.1
      x[2,1] = 0.1
      p  = (k,x,t0)
      t  = 0.0
      @btime model(du, u, p, t)
end

@testset "General Model" begin
      tstart = -30.0
      tend = 270.0
      tspan = (tstart, tend)

      # All in ground state and pumping hard. We shoud end up with
      # 1/2 in both ground and excited state.
      u0 = [1.0, 0, 0]
      N = length(u0)
      k = @MMatrix zeros(N,N)
      x = @MMatrix zeros(N,N)
      x[1,2] = 1.0
      x[2,1] = 1.0
      t0 = 0
      p = (k,x,t0)
      prob = ODEProblem(model, u0, tspan, p)
      sol = solve(prob)
      @test isapprox(sol[end][1], 0.5, atol=1e-4) &&
            isapprox(sol[end][2], 0.5, atol=1e-4)
      @test sum(sol[end]) ≈ 1.0

      # All in 2nd state. let flow to 1st and 3nd state. Let the flow
      # to state 3 be twice as fast.
      # We should end up with 1/3 and 2/3 in states 1 and 3 respectively.
      u0 = [0, 1.0, 0]
      N = length(u0)
      k = @MMatrix zeros(N,N)
      x = @MMatrix zeros(N,N)
      k[2,1] = 1.0
      k[2,3] = 2.0
      t0 = 0
      p = (k,x,t0)
      prob = ODEProblem(model, u0, tspan, p)
      sol = solve(prob)
      @test isapprox(sol[end][1], 1/3, atol=1e-4) &&
            isapprox(sol[end][3], 2/3, atol=1e-4)
      @test sum(sol[end]) ≈ 1.0

      # BenchmarkTools
      @btime solve($prob)
end


@testset "Population to Attenuation" begin
      # let's have a 4-state model
      #     1: ground state
      #     2: excited state (shows up in spectrum)
      #     3: excited state (shows up in spectrum)
      #     4: excited state (does not show up in spectrum and this state
      #                       does not contribute to SF intensity of states 2,3)
      # populations
      u = [0.7, 0.1, 0.05, 0.15]
      s1 = SFGReservoirs.State()
      s2 = SFGReservoirs.State(2, [1,3])
      s3 = SFGReservoirs.State(3, [1,2])
      states1 = (s2, s3)
      states2 = (s1, s2, s3)
      a1 = @MArray [0.0, 0.0]
      a2 = @MArray [0.0, 0.0, 0.0]
      SFGReservoirs.population_to_attenuation!(a1, u, states1)
      @test all(a1 .≈ [0.65, 0.75])
      SFGReservoirs.population_to_attenuation!(a2, u, states2)
      @test all(a2 .≈ [1.0, 0.65, 0.75])
end

begin
      u = [0.7, 0.1, 0.05, 0.15]
      s1 = SFGReservoirs.State()
      s2 = SFGReservoirs.State(2, [1,3])
      s3 = SFGReservoirs.State(3, [1,2])
      states1 = (s2, s3)
      states2 = (s1, s2, s3)
      a1 = @MArray [0.0, 0.0]
      a2 = @MArray [0.0, 0.0, 0.0]
      @btime SFGReservoirs.population_to_attenuation!(a1, u, states1)
      @btime SFGReservoirs.population_to_attenuation!(a2, u, states2)
end

# attenuated spectrum
begin
      x = collect(LinRange(2830, 3000, 512))
      a = [0.65, 1.1]
      y = zeros(ComplexF64, 512)
      A = [5, -3]
      ω = [2870, 2940]
      Γ = [8, 8]
      @btime SFGReservoirs.attenuated_sfspectrum!(y, x, a, A, ω, Γ)
      # plot(x, real(y))
      @btime sfspectrum!(y, x, A, ω, Γ)
      # plot!(x, real(y))
      nothing
end


begin
      n_states = 4
      u0 = [1.0, 0.0, 0.0, 0.0]
      tspan = (-30.0, 300.0)
      delaytimes = range(-30.0, 300.0, length=121)  |> collect
      wavenumbers = range(2830, 3000.0, length=512) |> collect

      A = [0.2, 5, -3]
      ω = [2845.0, 2870, 2940]
      Γ = [5, 8, 8.0]
      spec_params = (A, ω, Γ)

      p = [
            0.01,       # k21, k31
            0.03,       # k23
            0.02,       # k24, k34
            0.04,       # k32
            0.005,      # k41
            0.03,       # x12, x21
            0.01,       # t0
      ]

      k  = @MMatrix zeros(n_states, n_states)
      ak = @SMatrix [
            0 0 0 0
            1 0 2 3
            1 4 0 3
            5 0 0 0
      ]

      x  = @MMatrix zeros(n_states, n_states)
      ax = @SMatrix [
            0 6 0 0
            6 0 0 0
            0 0 0 0
            0 0 0 0
      ]

      t0 = p[7]

      SFGReservoirs.parameters_to_rate_matrix!(k, ak, p)
      SFGReservoirs.parameters_to_rate_matrix!(x, ax, p)

      params = (k,x,t0)
      prob = ODEProblem(model, u0, tspan, params)
      # saving at the delaytimes is performance critical
      @btime solve(prob; saveat=delaytimes)
      sol = solve(prob; saveat=delaytimes)

      display(plot(sol))

      # states represented in the spectrum
      s1 = SFGReservoirs.State()
      s2 = SFGReservoirs.State(2, [1,3])
      s3 = SFGReservoirs.State(3, [1,2])
      statez = [s1, s2, s3]

      # preallocations
      Y = zeros(ComplexF64, (length(wavenumbers), length(delaytimes)))
      a = @MArray zeros(length(statez))
      u = @MArray zeros(length(u0))

      @btime SFGReservoirs.solution_to_attenuated_spectra!(
            Y, a, u, wavenumbers, delaytimes, spec_params, sol, statez
      )

      display(heatmap(real.(Y)))
end