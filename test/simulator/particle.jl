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
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))

    ## DROP (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
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
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))

    ## SLIDE
    # initial conditions
    v1 = @SVector [1.0, 2.0, 0.0]
    q1 = @SVector [0.0, 0.0, 1.0]
    q0 = q1 - h * v1

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test isapprox.(sim.traj.q[end][3], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.traj.q[end] - sim.traj.q[end-1]) ./ h, 0.0, atol = 1.0e-6))

    ## DROP (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
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
    @test isapprox.(sim.traj.q[end][3], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.traj.q[end] - sim.traj.q[end-1]) ./ h, 0.0, atol = 1.0e-6))
end

@testset "Simulation: Particle (3D), Code-Gen Solver" begin
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
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8, solver = :mgs_solver),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))

    ## SLIDE
    # initial conditions
    v1 = @SVector [1.0, 2.0, 0.0]
    q1 = @SVector [0.0, 0.0, 1.0]
    q0 = q1 - h * v1

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8, solver = :mgs_solver),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test isapprox.(sim.traj.q[end][3], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.traj.q[end] - sim.traj.q[end-1]) ./ h, 0.0, atol = 1.0e-6))
end

@testset "Simulation: Particle (3D), Approximate Residual Jacobian" begin
    model = ContactControl.get_model("particle")

    dir = joinpath(pwd(), "src/dynamics/particle")
    model2 = deepcopy(particle)

    path_base = joinpath(dir, "dynamics/base.jld2")
    path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
    path_res = joinpath(dir, "quadratic_approx/residual.jld2")
    path_jac = joinpath(dir, "quadratic_approx/sparse_jacobians.jld2")
    path_linearized = joinpath(dir, "quadratic_approx/linearized.jld2")

    instantiate_base!(model2, path_base)

    expr_dyn = generate_dynamics_expressions(model2, derivs = true)
    save_expressions(expr_dyn, path_dyn, overwrite=true)
    instantiate_dynamics!(model2, path_dyn, derivs = true)

    expr_res, rz_sp, rθ_sp = generate_residual_expressions(model2, jacobians = :approx)
    save_expressions(expr_res, path_res, overwrite=true)
    @save path_jac rz_sp rθ_sp
    instantiate_residual!(model2, path_res, jacobians = :approx)
    model2.spa.rz_sp = copy(rz_sp)
    model2.spa.rθ_sp = copy(rθ_sp)

    r0 = zeros(num_var(model))
    model.res.r!(r0, 5.0 * ones(num_var(model)), 0.1 * ones(num_data(model)), 1.0)
    model.res.rz!(model.spa.rz_sp, 5.0 * ones(num_var(model)), 0.1 * ones(num_data(model)))
    model.res.rθ!(model.spa.rθ_sp, ones(num_var(model)), ones(num_data(model)))

    r1 = zeros(num_var(model2))
    model2.res.r!(r1, 5.0 * ones(num_var(model2)), 0.1 * ones(num_data(model2)), 1.0)
    model2.res.rz!(model2.spa.rz_sp, 5.0 * ones(num_var(model2)), 0.1 * ones(num_data(model2)))
    model2.res.rθ!(model2.spa.rθ_sp, ones(num_var(model2)), ones(num_data(model2)))

    norm(r0 - r1) < 1.0e-8
    norm(model.spa.rz_sp - model2.spa.rz_sp) < 1.0e-8
    norm(model.spa.rθ_sp - model2.spa.rθ_sp) < 1.0e-8

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    q0 = @SVector [0.0, 0.0, 1.0]
    q1 = @SVector [0.0, 0.0, 1.0]

    # simulator
    sim = ContactControl.simulator(model2, q0, q1, h, T,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))
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
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test abs(sim.traj.q[end][1]) < 0.05
    @test abs(sim.traj.q[end][2]) < 0.05
    @test abs(sim.traj.q[end][3]) < 0.001

    ## Forward Velocity (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
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
    @test abs(sim.traj.q[end][1]) < 0.05
    @test abs(sim.traj.q[end][2]) < 0.05
    @test abs(sim.traj.q[end][3]) < 0.001

    ## Drop
    # initial conditions
    q1 = @SVector [1.0, 0.5, 2.0]
    q0 = @SVector [1.0, 0.5, 2.0]

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test abs(sim.traj.q[end][1]) < 0.05
    @test abs(sim.traj.q[end][2]) < 0.05
    @test abs(sim.traj.q[end][3]) < 0.001

    ## DROP (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
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
    @test abs(sim.traj.q[end][1]) < 0.05
    @test abs(sim.traj.q[end][2]) < 0.05
    @test abs(sim.traj.q[end][3]) < 0.001
end

@testset "Simulation: Particle (3D) Quadratic Bowl, Approximate Residual Jacobian" begin
    dir = joinpath(pwd(), "src/dynamics/particle")
    model = deepcopy(particle_quadratic)

    path_base = joinpath(dir, "dynamics/base.jld2")
    path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
    path_res = joinpath(dir, "quadratic_approx/residual.jld2")
    path_jac = joinpath(dir, "quadratic_approx/sparse_jacobians.jld2")
    path_linearized = joinpath(dir, "quadratic_approx/linearized.jld2")

    instantiate_base!(model, path_base)

    expr_dyn = generate_dynamics_expressions(model, derivs = true)
    save_expressions(expr_dyn, path_dyn, overwrite=true)
    instantiate_dynamics!(model, path_dyn, derivs = true)

    expr_res, rz_sp, rθ_sp = generate_residual_expressions(model, jacobians = :approx)
    save_expressions(expr_res, path_res, overwrite=true)
    @save path_jac rz_sp rθ_sp
    instantiate_residual!(model, path_res, jacobians = :approx)
    model.spa.rz_sp = copy(rz_sp)
    model.spa.rθ_sp = copy(rθ_sp)

    # time
    h = 0.01
    T = 1000

    ## DROP
    # initial conditions
    q1 = @SVector [1.0, 0.5, 2.0]
    q0 = @SVector [1.1, 0.5, 2.0]

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test abs(sim.traj.q[end][1]) < 0.05
    @test abs(sim.traj.q[end][2]) < 0.05
    @test abs(sim.traj.q[end][3]) < 0.001

end

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
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))

    ## DROP (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
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
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))

    ## SLIDE
    # initial conditions
    v1 = @SVector [1.0, 0.0]
    q1 = @SVector [0.0, 1.0]
    q0 = q1 - h * v1

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test isapprox.(sim.traj.q[end][2], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.traj.q[end] - sim.traj.q[end-1]) ./ h, 0.0, atol = 1.0e-6))

    ## DROP (warmstart)
    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
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
    @test isapprox.(sim.traj.q[end][2], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.traj.q[end] - sim.traj.q[end-1]) ./ h, 0.0, atol = 1.0e-6))
end

@testset "Simulation: Particle (2D) On Slope" begin
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
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8, κ_tol = 1.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    status = ContactControl.simulate!(sim)
    @test status
    @test sim.traj.q[end][1] < 0.0
    @test sim.traj.q[end][2] < 0.0
end
