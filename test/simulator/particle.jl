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
#
# dir = joinpath(pwd(), "src/dynamics/particle")
# model = deepcopy(particle_quadratic)
#
# path_base = joinpath(dir, "dynamics/base.jld2")
# path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
# path_res = joinpath(dir, "quadratic_approx/residual.jld2")
# path_jac = joinpath(dir, "quadratic_approx/sparse_jacobians.jld2")
# path_linearized = joinpath(dir, "quadratic_approx/linearized.jld2")
#
# instantiate_base!(model, path_base)
#
# expr_dyn = generate_dynamics_expressions(model, derivs = true)
# save_expressions(expr_dyn, path_dyn, overwrite=true)
# instantiate_dynamics!(model, path_dyn, derivs = true)
#
# expr_res, rz_sp, rθ_sp = generate_residual_expressions(model, jacobians = :approx)
# save_expressions(expr_res, path_res, overwrite=true)
# @save path_jac rz_sp rθ_sp
# instantiate_residual!(model, path_res, jacobians = :approx)
#
# model.res.rz!(model.spa.rz_sp, rand(num_var(model)), rand(num_data(model)))
# model.res.rθ!(model.spa.rθ_sp, rand(num_var(model)), rand(num_data(model)))
#
# rθzz = zeros(num_var(model), num_data(model))
# num_data(model)
# idx = collect([(1:2model.dim.q + model.dim.u + model.dim.w)..., num_data(model)])
# rθzz[1:model.dim.q, idx] =
# model.dyn.dθ(1.0, q0s, q1s, u1s, w1s, λ1s, q2s)
#
#
#
#
#
#
#
#
#
#
# nq = model.dim.q
# nu = model.dim.u
# nw = model.dim.w
# nc = model.dim.c
# nb = model.dim.b
# ncf = nc * dim(model.env)
#
# # Declare variables
# @variables q0[1:nq]
# @variables q1[1:nq]
# @variables u1[1:nu]
# @variables w1[1:nw]
# @variables λ1[1:ncf]
# @variables q2[1:nq]
# @variables h
#
# # Expressions
# expr = Dict{Symbol, Expr}()
#
# # Dynamics
# d = dynamics(model, h, q0, q1, u1, w1, λ1, q2)
# d = Symbolics.simplify.(d)
#
# # Functions
# dd = eval(build_function(d, h, q0, q1, u1, w1, λ1, q2)[1])
#
# dθ = Symbolics.jacobian(d, [q0; q1; u1; w1])#, simplify = true)
#
# # expr[:dq2]  = build_function(dq2, h, q0, q1, u1, w1, λ1, q2)[1]
# # expr[:dλ1]  = build_function(dλ1, h, q0, q1, u1, w1, λ1, q2)[1]
# ddθ = eval(build_function(dθ, h, q0, q1, u1, w1, λ1, q2)[1])
# ddθ(1.0, q0s, q1s, u1s, w1s, λ1s, q2s)
#
# if derivs
#     dq2 = Symbolics.jacobian(d, q2, simplify = true)
#     dλ1 = Symbolics.jacobian(d, λ1, simplify = true)
#     dθ = Symbolics.jacobian(d, [q0; q1; u1; w1; h])#, simplify = true)
#
#     expr[:dq2]  = build_function(dq2, h, q0, q1, u1, w1, λ1, q2)[1]
#     expr[:dλ1]  = build_function(dλ1, h, q0, q1, u1, w1, λ1, q2)[1]
#     expr[:dθ]   = build_function(dθ, h, q0, q1, u1, w1, λ1, q2)[1]
# end
#



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
