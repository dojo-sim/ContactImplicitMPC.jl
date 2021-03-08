@testset "Simulation: Particle (3D)" begin
    include("particle.jl")

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    q0 = @SVector [0.0, 0.0, 1.0]
    q1 = @SVector [0.0, 0.0, 1.0]

    # simulator
    sim = simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = SimulatorOptions(warmstart = false))

    # simulate
    status = simulate!(sim)
    @show sim.q[end]
    @test status
    @test all(isapprox.(sim.q[end], 0.0, atol = 1.0e-6))

    ## DROP (warmstart)
    # simulator
    sim = simulator(model, q0, q1, h, T,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_tol = 1.0e-8),
        sim_opts = SimulatorOptions(
            warmstart = false,
            z_warmstart = 0.001,
            κ_warmstart = 0.001))

    # simulate
    status = simulate!(sim)
    @show sim.q[end]
    @test status
    @test all(isapprox.(sim.q[end], 0.0, atol = 1.0e-6))

    ## SLIDE
    # initial conditions
    v1 = @SVector [1.0, 2.0, 0.0]
    q1 = @SVector [0.0, 0.0, 1.0]
    q0 = q1 - h * v1

    # simulator
    sim = simulator(model, q0, q1, h, T,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = SimulatorOptions(warmstart = false))

    # simulate
    status = simulate!(sim)
    @show sim.q[end]
    @test status
    @test isapprox.(sim.q[end][3], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.q[end] - sim.q[end-1]) ./ h, 0.0, atol = 1.0e-6))

    ## DROP (warmstart)
    # simulator
    sim = simulator(model, q0, q1, h, T,
        rz = rz_sp,
        rθ = rθ_sp,
        ip_opts = InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_tol = 1.0e-8),
        sim_opts = SimulatorOptions(
            warmstart = false,
            z_warmstart = 0.001,
            κ_warmstart = 0.001))

    # simulate
    status = simulate!(sim)
    @show sim.q[end]
    @test status
    @test isapprox.(sim.q[end][3], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.q[end] - sim.q[end-1]) ./ h, 0.0, atol = 1.0e-6))
end
