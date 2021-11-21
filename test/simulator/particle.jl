@testset "Simulation: Particle (3D)" begin
    s = ContactImplicitMPC.get_simulation("particle", "flat_3D_lc", "flat_lc")
    model = s.model

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    v1 = [0.0, 0.0, 0.0]
    q1 = [0.0, 0.0, 1.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    status = simulate!(sim, q1, v1)
    @test status
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))

    ## SLIDE
    # initial conditions
    v1 = [1.0, 2.0, 0.0]
    q1 = [0.0, 0.0, 1.0]

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test isapprox.(sim.traj.q[end][3], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.traj.q[end] - sim.traj.q[end-1]) ./ h, 0.0, atol = 1.0e-6))
end

@testset "Simulation: Particle (3D), Code-Gen Solver" begin
    s = ContactImplicitMPC.get_simulation("particle", "flat_3D_lc", "flat_lc")
    model = s.model

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    v1 = [0.0, 0.0, 0.0]
    q1 = [0.0, 0.0, 1.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))

    ## SLIDE
    # initial conditions
    v1 = [1.0, 2.0, 0.0]
    q1 = [0.0, 0.0, 1.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test isapprox.(sim.traj.q[end][3], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.traj.q[end] - sim.traj.q[end-1]) ./ h, 0.0, atol = 1.0e-6))
end

@testset "Simulation: Particle (3D), Approximate Residual Jacobian" begin

    s = ContactImplicitMPC.get_simulation("particle", "flat_3D_lc", "flat_lc")
    model = s.model

    dir_dyn = joinpath(module_dir(), "src/dynamics/particle")
    dir_sim = joinpath(module_dir(), "src/simulation/particle")
    model2 = deepcopy(ContactImplicitMPC.particle)
    env2 = environment_3D_flat()
    s2 = Simulation(model2, env2)

    path_base = joinpath(dir_dyn, "dynamics/base.jld2")
    path_dyn = joinpath(dir_dyn, "dynamics/dynamics.jld2")
    path_res = joinpath(dir_sim, "quadratic_approx/residual.jld2")
    path_jac = joinpath(dir_sim, "quadratic_approx/sparse_jacobians.jld2")

    instantiate_base!(s2.model, path_base)

    expr_dyn = generate_dynamics_expressions(s2.model, derivs = true)
    save_expressions(expr_dyn, path_dyn, overwrite=true)
    instantiate_dynamics!(s2.model, path_dyn, derivs = true)

    expr_res, rz_sp, rθ_sp = generate_residual_expressions(s2.model, s2.env, jacobians = :approx)
    save_expressions(expr_res, path_res, overwrite=true)
    @save path_jac rz_sp rθ_sp
    instantiate_residual!(s2, path_res, path_jac, jacobians = :approx)

    r0 = zeros(num_var(s.model, s.env))
    s.res.r!(r0, 5.0 * ones(num_var(s.model, s.env)), 0.1 * ones(num_data(model)), 1.0)
    s.res.rz!(s.rz, 5.0 * ones(num_var(s.model, s.env)), 0.1 * ones(num_data(model)))
    s.res.rθ!(s.rθ, ones(num_var(s.model, s.env)), ones(num_data(model)))

    r1 = zeros(num_var(s2.model, s2.env))
    s2.res.r!(r1, 5.0 * ones(num_var(s2.model, s2.env)), 0.1 * ones(num_data(s2.model)), 1.0)
    s2.res.rz!(s2.rz, 5.0 * ones(num_var(s2.model, s2.env)), 0.1 * ones(num_data(s2.model)))
    s2.res.rθ!(s2.rθ, ones(num_var(s2.model, s2.env)), ones(num_data(s2.model)))

    @test norm(r0 - r1) < 1.0e-8
    @test norm(s.rz - s2.rz) < 1.0e-8
    @test norm(s.rθ - s2.rθ) < 1.0e-8

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    v1 = [0.0, 0.0, 0.0]
    q1 = [0.0, 0.0, 1.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))
end

@testset "Simulation: Particle (3D) in Quadratic Bowl" begin
    s = get_simulation("particle", "quadratic_bowl_3D_lc", "quadratic")
    s.model.μ_world = 0.1

    # time
    h = 0.01
    T = 1000

    # initial conditions
    q1 = [1.0, 0.5, 2.0]
    v1 = [0.1; 0.0; 0.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    @time status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test abs(sim.traj.q[end][1]) < 0.05
    @test abs(sim.traj.q[end][2]) < 0.05
    @test abs(sim.traj.q[end][3]) < 0.001

    ## Drop
    # initial conditions
    q1 = [1.0; 0.5; 2.0]
    v1 = [0.0, 0.0, 0.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test abs(sim.traj.q[end][1]) < 0.05
    @test abs(sim.traj.q[end][2]) < 0.05
    @test abs(sim.traj.q[end][3]) < 0.001
end

@testset "Simulation: Particle (3D) Quadratic Bowl, Approximate Residual Jacobian" begin
    dir_dyn = joinpath(module_dir(), "src/dynamics/particle")
    dir_sim = joinpath(module_dir(), "src/simulation/particle")
    model2 = deepcopy(ContactImplicitMPC.particle)
    model2.μ_world = 0.1
    env2 = deepcopy(ContactImplicitMPC.quadratic_bowl_3D_lc)
    s2 = Simulation(model2, env2)

    path_base = joinpath(dir_dyn, "dynamics/base.jld2")
    path_dyn = joinpath(dir_dyn, "dynamics/dynamics.jld2")
    path_res = joinpath(dir_sim, "quadratic_approx/residual.jld2")
    path_jac = joinpath(dir_sim, "quadratic_approx/sparse_jacobians.jld2")

    instantiate_base!(s2.model, path_base)

    expr_dyn = generate_dynamics_expressions(s2.model, derivs = true)
    save_expressions(expr_dyn, path_dyn, overwrite=true)
    instantiate_dynamics!(s2.model, path_dyn, derivs = true)

    expr_res, rz_sp, rθ_sp = generate_residual_expressions(s2.model, s2.env, jacobians = :approx)
    save_expressions(expr_res, path_res, overwrite=true)
    @save path_jac rz_sp rθ_sp
    instantiate_residual!(s2, path_res, path_jac, jacobians = :approx)

    # time
    h = 0.01
    T = 1000

    ## DROP
    # initial conditions
    q1 = [1.0, 0.5, 2.0]
    v1 = [0.1, 0.0, 0.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s2, T, h=h)

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test abs(sim.traj.q[end][1]) < 0.05
    @test abs(sim.traj.q[end][2]) < 0.05
    @test abs(sim.traj.q[end][3]) < 0.001
end

@testset "Simulation: Particle (2D)" begin
    s = get_simulation("particle_2D", "flat_2D_lc", "flat_lc")

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    v1 = [0.0, 0.0]
    q1 = [0.0, 1.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))

    ## SLIDE
    # initial conditions
    v1 = [1.0, 0.0]
    q1 = [0.0, 1.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test isapprox.(sim.traj.q[end][2], 0.0, atol = 1.0e-6)
    @test all(isapprox.((sim.traj.q[end] - sim.traj.q[end-1]) ./ h, 0.0, atol = 1.0e-6))
end

@testset "Simulation: Particle (2D) On Slope" begin
    s = get_simulation("particle_2D", "slope1_2D_lc", "slope")
    s.model.μ_world = 0.1

    # time
    h = 0.01
    T = 100

    ## DROP
    # initial conditions
    v1 = [0.0, 0.0]
    q1 = [0.0, 1.0]

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
    @test sim.traj.q[end][1] < 0.0
    @test sim.traj.q[end][2] < 0.0
end

