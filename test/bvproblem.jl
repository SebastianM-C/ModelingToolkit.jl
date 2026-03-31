using ModelingToolkit
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using BoundaryValueDiffEqMIRK
using Test

solvers = [MIRK4]

@testset "Simple DAE" begin
    @variables u(t) v(t)
    eqs = [D(u) ~ v,
           0 ~ u^2 + v^2 - 1]
    @mtkcompile sys = System(eqs, t)

    tspan = (0.0, 1.0)
    u0map = [u => 0.5]
    v0 = sqrt(1 - 0.5^2)

    prob = ODEProblem(sys, u0map, tspan)
    osol = solve(prob, Rodas5P())

    bvp = BVProblem(sys, u0map, tspan; guesses = [v => v0])
    for solver in solvers
        sol = solve(bvp, solver(), dt = 0.01)
        # algebraic constraint u² + v² = 1 should be approximately satisfied
        @test all(isapprox.(sol[u .^ 2 + v .^ 2], 1.0; atol = 1e-2))
        # endpoint should match ODE reference
        @test isapprox(sol(1.0, idxs=[u, v]), osol(1.0, idxs=[u, v]); atol = 0.01)
    end
end

@testset "Cartesian Pendulum" begin
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t)
    eqs = [D(D(x)) ~ λ * x
            D(D(y)) ~ λ * y - g
            x^2 + y^2 ~ 1]
    @mtkcompile pend = System(eqs, t)

    tspan = (0.0, 1.5)
    u0map = [x => 1, y => 0]
    pmap = [g => 1]
    guess = [λ => 1]

    prob = ODEProblem(pend, [u0map; pmap], tspan; guesses = guess)
    osol = solve(prob, Rodas5P())

    bvp = BVProblem(pend, [u0map; pmap], tspan; guesses = guess)

    for solver in solvers
        sol = solve(bvp, solver(), dt = 0.001)
        @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
        conditions = getfield.(equations(pend)[3:end], :rhs)
        @test isapprox([sol[conditions][1]; sol[x][1] - 1; sol[y][1]], zeros(5), atol = 0.001)
    end
end
