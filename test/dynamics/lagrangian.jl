@testset "Lagrangian Dynamics" begin
	# Code gen
	dir = joinpath(module_dir(), "src", "dynamics", "rigidbody")
	# include(joinpath(dir, "model.jl"))

	# Analytical Model
	model = deepcopy(ContactImplicitMPC.rigidbody)

	path_base = joinpath(dir, "test_dynamics/base.jld2")
	path_dyn = joinpath(dir, "test_dynamics/dynamics.jld2")

	expr_base = generate_base_expressions(model,
		M_analytical = true,
		mapping = ContactImplicitMPC.G_func,
		nv = model.dim.q - 1)

	save_expressions(expr_base, path_base, overwrite=true)
	instantiate_base!(model, path_base)

	expr_dyn = generate_dynamics_expressions(model, nv = model.dim.q - 1)
	save_expressions(expr_dyn, path_dyn, overwrite=true)
	instantiate_dynamics!(model, path_dyn)

	# Lagrangian Model
	model_lag = deepcopy(ContactImplicitMPC.rigidbody)
	path_base_lag = joinpath(dir, "test_dynamics/base_lag.jld2")
	path_dyn_lag = joinpath(dir, "test_dynamics/dynamics_lag.jld2")

	expr_base_lag = generate_base_expressions(model_lag,
		M_analytical = false,
		mapping = ContactImplicitMPC.G_func,
		nv = model.dim.q - 1)

	save_expressions(expr_base, path_base_lag, overwrite=true)
	instantiate_base!(model_lag, path_base_lag)

	expr_dyn_lag = generate_dynamics_expressions(model_lag, nv = model.dim.q - 1)
	save_expressions(expr_dyn_lag, path_dyn_lag, overwrite=true)
	instantiate_dynamics!(model_lag, path_dyn_lag)

	# Data
	env = environment_3D_flat()
	nq = model.dim.q
	nv = model.dim.q - 1
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	h = 0.01

	# initial conditions
	r0 = [0.0; 0.0; 1.0]
	v0 = [7.5; 5.0; 0.0]

	quat0 = [1.0; 0.0; 0.0; 0.0]
	ω0 = [0.0; 0.0; 0.0]

	q0 = SVector{model.dim.q}([r0; quat0])
	q1 = SVector{model.dim.q}([r0 + v0 * h; 0.5 * h * L_multiply(quat0) * [sqrt((2.0 / h)^2.0 - ω0' * ω0); ω0]])
	q2 = SVector{model.dim.q}([r0 + v0 * 2h; 0.5 * 2h * L_multiply(quat0) * [sqrt((2.0 / (2h))^2.0 - ω0' * ω0); ω0]])
	q̇0 = zeros(SVector{nv})
	u1 = zeros(SVector{nu})
	w1 = zeros(SVector{nw})
	λ1 = zeros(SVector{nc * dim(env)})

	# Test equality
	@test norm(ContactImplicitMPC.M_fast(model, q0) -
		       ContactImplicitMPC.M_fast(model_lag, q0)) < 1e-10
	@test norm(ContactImplicitMPC.C_fast(model, q0, q̇0) -
			   ContactImplicitMPC.C_fast(model_lag, q0, q̇0)) < 1e-10
	# @test norm(ContactImplicitMPC.d_fast(model, [h], q0, q1, u1, w1, λ1, q2) -
	# 		   ContactImplicitMPC.d_fast(model_lag, [h], q0, q1, u1, w1, λ1, q2)) < 1e-10 # Needs to be fixed
end
