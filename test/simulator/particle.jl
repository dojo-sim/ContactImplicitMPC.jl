@testset "Simulation: Particle (3D)" begin
    res_path = joinpath(@__DIR__, "../../src/dynamics")
    include(joinpath(res_path, "particle/model.jl"))
    model = deepcopy(particle)
    ContactControl.instantiate_base!(model, joinpath(res_path, "particle/flat/base.jld2"))
    ContactControl.instantiate_dynamics!(model, joinpath(res_path, "particle/flat/dynamics.jld2"))
    ContactControl.instantiate_residual!(model, joinpath(res_path, "particle/flat/residual.jld2"))
    @load joinpath(res_path, "particle/flat/sparse_jacobians.jld2") rz_sp rθ_sp

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    q0 = @SVector [0.0, 0.0, 1.0]
    q1 = @SVector [0.0, 0.0, 1.0]

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.q[end], 0.0, atol = 1.0e-6))

    ## DROP (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(
            warmstart = true,
            z_warmstart = 0.001,
            κ_warmstart = 0.001))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.q[end], 0.0, atol = 1.0e-6))

    ## SLIDE
    # initial conditions
    v1 = @SVector [1.0, 2.0, 0.0]
    q1 = @SVector [0.0, 0.0, 1.0]
    q0 = q1 - h * v1

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test isapprox.(sim.q[end][3], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.q[end] - sim.q[end-1]) ./ h, 0.0, atol = 1.0e-6))

    ## DROP (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(
            warmstart = true,
            z_warmstart = 0.001,
            κ_warmstart = 0.001))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test isapprox.(sim.q[end][3], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.q[end] - sim.q[end-1]) ./ h, 0.0, atol = 1.0e-6))
end

@testset "Simulation: Particle (3D) in Quadratic Bowl" begin
    res_path = joinpath(@__DIR__, "../../src/dynamics")
    include(joinpath(res_path, "particle/model.jl"))
    model = deepcopy(particle_quadratic)
    ContactControl.instantiate_base!(model, joinpath(res_path, "particle/quadratic/base.jld2"))
    ContactControl.instantiate_dynamics!(model, joinpath(res_path, "particle/quadratic/dynamics.jld2"))
    ContactControl.instantiate_residual!(model, joinpath(res_path, "particle/quadratic/residual.jld2"))
    @load joinpath(res_path, "particle/quadratic/sparse_jacobians.jld2") rz_sp rθ_sp

    # time
    h = 0.01
    T = 1000

    ## Forwward Velocity
    # initial conditions
    q1 = @SVector [1.0, 0.5, 2.0]
    q0 = @SVector [1.1, 0.5, 2.0]

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.q[end], 0.0, atol = 1.0e-2))

    ## Forward Velocity (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(
            warmstart = true,
            z_warmstart = 0.001,
            κ_warmstart = 0.001))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.q[end], 0.0, atol = 1.0e-2))

    ## Drop
    # initial conditions
    q1 = @SVector [1.0, 0.5, 2.0]
    q0 = @SVector [1.0, 0.5, 2.0]

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.q[end][3], 0.0, atol = 1.0e-3))
    @test all(isapprox.((sim.q[end] - sim.q[end-1]) ./ h, 0.0, atol = 1.0e-3))

    ## DROP (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(
            warmstart = true,
            z_warmstart = 0.001,
            κ_warmstart = 0.001))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.q[end][3], 0.0, atol = 1.0e-3))
    @test all(isapprox.((sim.q[end] - sim.q[end-1]) ./ h, 0.0, atol = 1.0e-3))
end
