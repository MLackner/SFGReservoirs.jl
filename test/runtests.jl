using SFGReservoirs
using Test

using DifferentialEquations
using SparseArrays


@testset "sfspectra.jl" begin
    x = 2850:3020

    # Two oscillators with opposite signs should
    # yield no resonance
    A = [1, -1]
    ω = [2900, 2900]
    Γ = [5, 5]
    nr = 5
    y = sfspectrum(x, A, ω, Γ, nr=nr)
    @test length(unique(y)) == 1
end


@testset "reservoirs.jl" begin

    tstart = -30.0
    tend = 270.0
    tspan = (tstart, tend)

    # Case: Pumping from 1 --> 3 and nothing else shoud yield .5 in 1 and 3
    u0 = [1.0, 0, 0, 0, 0]

    p = spzeros(5,5)
    p[1,3] = 1.0

    prob = ODEProblem(model_A, u0, tspan, p)
    sol = solve(prob);
    @test isapprox(sol[end][1], 0.5, atol=1e-5) &&
          isapprox(sol[end][3], 0.5, atol=1e-5)

    # Case: Start with all in 2 and let there be an equilibrium in states 2,3,4
    #       All three states should hold 1/3
    u0 = [0, 1.0, 0, 0, 0]

    p = spzeros(5,5)
    p[2,3] = 1.0
    p[2,4] = 1.0
    p[3,2] = 1.0
    p[3,4] = 1.0
    p[4,2] = 1.0
    p[4,3] = 1.0

    prob = ODEProblem(model_A, u0, tspan, p)
    sol = solve(prob);
    @test isapprox(sol[end][2], 1/3, atol=1e-5) &&
          isapprox(sol[end][3], 1/3, atol=1e-5) &&
          isapprox(sol[end][4], 1/3, atol=1e-5)

    # Case: Start with all in 3 and pump. We should end up with .5 in 1 and 3.
    u0 = [0, 0, 1.0, 0, 0]

    p = spzeros(5,5)
    p[1,3] = 1.0

    prob = ODEProblem(model_A, u0, tspan, p)
    sol = solve(prob);
    @test isapprox(sol[end][1], 0.5, atol=1e-5) &&
          isapprox(sol[end][3], 0.5, atol=1e-5)

    # Case: Pump all 1-->3 dissipate equally to 2,3,4 and drain to h.
    #       We should end up with all in h.
    u0 = [1.0, 0, 0, 0, 0]

    p = spzeros(5,5)
    p[1,3] = 1.0
    p[2,3] = 1.0
    p[2,4] = 1.0
    p[3,2] = 1.0
    p[3,4] = 1.0
    p[4,2] = 1.0
    p[4,3] = 1.0
    p[2,5] = 1.0
    p[3,5] = 1.0
    p[4,5] = 1.0

    prob = ODEProblem(model_A, u0, tspan, p)
    sol = solve(prob);
    @test isapprox(sol[end][5], 1.0, atol=1e-4)

    # Case: Pump all 1-->3 dissipate equally to 2,3,4 and drain to h-->1
    #       We should end up with all in 1.
    u0 = [1.0, 0, 0, 0, 0]
    p = spzeros(5,5)
    p[1,3] = 1.0
    p[2,3] = 1.0
    p[2,4] = 1.0
    p[3,2] = 1.0
    p[3,4] = 1.0
    p[4,2] = 1.0
    p[4,3] = 1.0
    p[5,1] = 1.0
    p[2,5] = 1.0
    p[3,5] = 1.0
    p[4,5] = 1.0

    prob = ODEProblem(model_A, u0, tspan, p)
    sol = solve(prob);
    @test isapprox(sol[end][1], 1.0, atol=1e-4)
end


@testset "General Model" begin
      tstart = -30.0
      tend = 270.0
      tspan = (tstart, tend)

      # All in ground state and pumping hard. We shoud end up with
      # 1/2 in both ground and excited state.
      u0 = [0, 1.0, 0]
      p = spzeros(3,3)
      p[2,3] = 1.0
      prob = ODEProblem(model, u0, tspan, p)
      sol = solve(prob)
      @test isapprox(sol[end][2], 0.5, atol=1e-4) &&
            isapprox(sol[end][3], 0.5, atol=1e-4)
      @test sum(sol[end]) ≈ 1.0

      # All in excited state and pumping hard. We shoud end up with
      # 1/2 in both ground and excited state.
      u0 = [0, 0.0, 1.0]
      p = spzeros(3,3)
      p[2,3] = 1.0
      prob = ODEProblem(model, u0, tspan, p)
      sol = solve(prob)
      @test isapprox(sol[end][2], 0.5, atol=1e-4) &&
            isapprox(sol[end][3], 0.5, atol=1e-4)
      @test sum(sol[end]) ≈ 1.0

      # All in ground state. Pump and let equilibrate in all additional
      # modes. Including the ground state we should have 1/10 everywhere.
      N = 11 # total number of states
      u0 = zeros(N)
      u0[2] = 1.0
      p = spzeros(N,N)
      p[2,3] = 1.0
      for i = 3:N, j=3:N
            i == j && continue
            p[i,j] = 1.0
      end
      prob = ODEProblem(model, u0, tspan, p, maxiter=1)
      sol = solve(prob)
      for i = 2:11
            @test isapprox(sol[end][i], 1/10, atol=1e-4)
      end
      @test sum(sol[end]) ≈ 1.0

      # All in excited state and pumping hard. We shoud end up with
      # 1/2 in both ground and excited state.
      u0 = [0, 0.0, 1.0]
      p = spzeros(3,3)
      p[2,3] = 1.0
      prob = ODEProblem(model, u0, tspan, p)
      sol = solve(prob)
      @test isapprox(sol[end][2], 0.5, atol=1e-4) &&
            isapprox(sol[end][3], 0.5, atol=1e-4)
      @test sum(sol[end]) ≈ 1.0

      # All in ground state. Pump and let equilibrate in all additional
      # modes + heat bath. Heat bath routes to ground state again. We should
      # end up with all in ground state.
      N = 6 # total number of states
      u0 = zeros(N)
      u0[2] = 1.0
      p = spzeros(N,N)
      p[2,3] = 1.0
      for i = 3:N
            p[i,1] = 1.0 # all to heat bath
      end
      p[1,2] = 1.0 # heat bath to ground state
      for i = 3:N, j=3:N
            i == j && continue
            p[i,j] = 1.0
      end
      global p
      prob = ODEProblem(model, u0, tspan, p, maxiter=1)
      sol = solve(prob)
      @show sol[end]
      @test sol[end][2] ≈ 1.0
      @test sum(sol[end]) ≈ 1.0
end
