@testset "Simulation: Particle (3D)" begin
    model = ContactControl.get_model("particle")

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
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
    model = get_model("particle", surf = "quadratic")

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
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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

model = ContactControl.get_model("particle_2D")

@testset "Simulation: Particle (2D)" begin
    model = ContactControl.get_model("particle_2D")

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    q0 = @SVector [0.0, 1.0]
    q1 = @SVector [0.0, 1.0]

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
    v1 = @SVector [1.0, 0.0]
    q1 = @SVector [0.0, 1.0]
    q0 = q1 - h * v1

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test isapprox.(sim.q[end][2], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.q[end] - sim.q[end-1]) ./ h, 0.0, atol = 1.0e-6))

    ## DROP (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
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
    @test isapprox.(sim.q[end][2], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.q[end] - sim.q[end-1]) ./ h, 0.0, atol = 1.0e-6))
end

@testset "Simulation: Particle (2D) slope" begin
    model = get_model("particle_2D", surf = "slope")

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    q0 = @SVector [0.0, 1.0]
    q1 = @SVector [0.0, 1.0]

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
        ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test sim.q[end][1] < 0.0
    @test sim.q[end][2] < 0.0
end
