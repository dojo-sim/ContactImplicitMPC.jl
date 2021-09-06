mutable struct Simulation
	model::ContactModel
	env::Environment{<:World,<:FrictionCone}
	con::ContactMethods
	res::ResidualMethods
	rz::Any
	rθ::Any
end

function Simulation(model::ContactModel, env::Environment)
	Simulation(model, env, ContactMethods(), ResidualMethods(), zeros(0, 0), zeros(0, 0))
end

function get_simulation(model::String, env::String, sim_name::String;
		model_variable_name::String = model,
		dynamics_name::String = "dynamics",
		gen_base = true,
		gen_dyn = true,
		approx = false)

	#TODO: assert model exists

	dir_model = joinpath(@__DIR__, "..", "dynamics", model)
	dir_sim   = joinpath(@__DIR__, "..", "simulation", model, sim_name)

	dir_base = joinpath(dir_model, dynamics_name, "base.jld2")
	dir_dyn = joinpath(dir_model, dynamics_name, "dynamics.jld2")
	dir_res = joinpath(dir_sim, "residual.jld2")
	dir_jac = joinpath(dir_sim, "jacobians.jld2")

	model = deepcopy(eval(Symbol(model_variable_name)))
	env = deepcopy(eval(Symbol(env)))
	sim = Simulation(model, env)

	instantiate_base!(sim.model, dir_base)

	instantiate_dynamics!(sim.model, dir_dyn,
		derivs = approx)

	instantiate_residual!(sim,
		dir_res, dir_jac,
		jacobians = (approx ? :approx : :full))

	return sim
end

function z_initialize!(z, model::ContactModel, env::Environment{<:World,LinearizedCone}, q1)
	iq2 = index_q2(model, env, nquat = 0)
	z .= 1.0
	z[iq2] = deepcopy(q1)
end

function z_initialize!(z, model::ContactModel, env::Environment{<:World,NonlinearCone}, q1)
	iq2 = index_q2(model, env, nquat = 0)
	ib1 = index_b1(model, env, nquat = 0)
	iψ1 = index_ψ1(model, env, nquat = 0)
	iη1 = index_η1(model, env, nquat = 0)
	is2 = index_s2(model, env, nquat = 0)

    z .= 0.1
    z[iq2] = deepcopy(q1)

	# second-order cone initializations
	z[ib1] .= 0.1 # primal cones: vector part # TODO redundant
	z[iψ1] .= 1.0 # primal cones: scalar part
	z[iη1] .= 0.1 # dual cones: vector part # TODO redundant
	z[is2] .= 1.0 # dual cones: scalar part
	return z
end

function z_warmstart!(z, model::ContactModel, env::Environment{<:World,LinearizedCone}, q, a, idx_ineq)
	iq2 = index_q2(model, env, nquat = 0)
	z[iq2] = q
	# TODO SIMON sort out idx_ineq
	z[idx_ineq] .+= a * rand(length(idx_ineq))
	nothing
end

function z_warmstart!(z, model::ContactModel, env::Environment{<:World,NonlinearCone}, q, a, idx_ineq)
	@warn "warm start for second-order cone not implemented"
	z_initialize!(z, model, env, q)
	nothing
end

function θ_initialize!(θ, model::ContactModel, q0, q1, u1, w1, μ, h)
	nθ = num_data(model)

	iq0 = index_q0(model)
	iq1 = index_q1(model)
	iu1 = index_u1(model)
	iw1 = index_w1(model)
	iμ  = index_μ(model)
	ih  = index_h(model)

	θ = zeros(nθ)
	θ[iq0] .= deepcopy(q0)
	θ[iq1] .= deepcopy(q1)
	θ[iu1] .= deepcopy(u1)
	θ[iw1] .= deepcopy(w1)
	θ[iμ] .= deepcopy(μ)
	θ[ih] .= deepcopy(h)
	return θ
end

function E_func(model::ContactModel, env::Environment{<:World,LinearizedCone})
	nc = model.dim.c
	nb = nc * friction_dim(env)
	SMatrix{nc, nb}(kron(Diagonal(ones(nc)), ones(1, Int(nb / nc))))
end

function residual(model::ContactModel, env::Environment{<:World,LinearizedCone}, z, θ, κ)
	nc = model.dim.c
	nb = nc * friction_dim(env)
	nf = friction_dim(env)
	np = dim(env)

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	ϕ = ϕ_func(model, env, q2)

	k = kinematics(model, q2)
	λ1 = contact_forces(model, env, γ1, b1, q2, k)
	Λ1 = transpose(J_func(model, env, q2)) * λ1 #@@@@ maybe need to use J_fast
	vT_stack = velocity_stack(model, env, q1, q2, k, h)
	ψ_stack = transpose(E_func(model, env)) * ψ1

	[model.dyn.d(h, q0, q1, u1, w1, Λ1, q2);
	 s1 - ϕ;
	 η1 - vT_stack - ψ_stack;
	 s2 .- (μ[1] * γ1 .- E_func(model, env) * b1);
	 γ1 .* s1 .- κ;
	 b1 .* η1 .- κ;
	 ψ1 .* s2 .- κ]
end

function residual(model::ContactModel, env::Environment{<:World,NonlinearCone}, z, θ, κ)
	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	ϕ = ϕ_func(model, env, q2)
	k = kinematics(model, q2)
	λ1 = contact_forces(model, env, γ1, b1, q2, k)
	vT = velocity_stack(model, env, q1, q2, k, h)

	nc = model.dim.c
	nf = friction_dim(env)
	ne = dim(env)

	[dynamics(model, env, h, q0, q1, u1, w1, λ1, q2);
	 s1 - ϕ;
	 vcat([η1[(i - 1) * nf .+ (1:nf)] - vT[(i - 1) * nf .+ (1:nf)] for i = 1:nc]...);
	 s2 - μ[1] * γ1;
	 γ1 .* s1 .- κ;
	 vcat([
	 	second_order_cone_product( # [ψ1; b1] ∘ [s2; η1] - [κ; 0...0]
			[ψ1[i]; b1[(i-1) * nf .+ (1:nf)]],
			[s2[i]; η1[(i-1) * nf .+ (1:nf)]]
			) - [κ; zeros(nf)]
		for i = 1:nc]...)
	 ]
end

function ResidualMethods()
	function f()
		error("Not Implemented: use instantiate_residual!")
		return nothing
	end
	return ResidualMethods(fill(f, 3)...)
end
