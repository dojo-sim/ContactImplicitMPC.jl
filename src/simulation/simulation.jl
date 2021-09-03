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

function unpack_θ(model::ContactModel, θ)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c

	# Parameters
	off = 0
	q0 = θ[off .+ (1:nq)]
	off += nq
	q1 = θ[off .+ (1:nq)]
	off += nq
	u1 =  θ[off .+ (1:nu)]
	off += nu
	w1 =  θ[off .+ (1:nw)]
	off += nw
	μ =   θ[off .+ (1:1)]
	off += 1
	h =   θ[off .+ (1:1)]
	return q0, q1, u1, w1, μ, h
end

function pack_θ(model::ContactModel, q0, q1, u1, w1, μ, h)
	return [q0; q1; u1; w1; μ; h]
end

function unpack_z(model::ContactModel, env::Environment{<:World,LinearizedCone}, z)
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = nc * friction_dim(env)

	# system variables
	off = 0
	q2 =  z[off .+ (1:nq)]
	off += nq
	γ1 =  z[off .+ (1:nc)]
	off += nc
	b1 =  z[off .+ (1:nb)]
	off += nb
	ψ1 =  z[off .+ (1:nc)]
	off += nc
	η1 =  z[off .+ (1:nb)]
	off += nb
	s1 = z[off .+ (1:nc)]
	off += nc
	s2 = z[off .+ (1:nc)]
	off += nc
	return q2, γ1, b1, ψ1, η1, s1, s2
end

function unpack_z(model::ContactModel, env::Environment{<:World,NonlinearCone}, z)
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = nc * friction_dim(env)

	# system variables
	off = 0
	q2 =  z[off .+ (1:nq)]
	off += nq
	γ1 =  z[off .+ (1:nc)]
	off += nc
	b1 =  z[off .+ (1:nb)]
	off += nb
	η1 =  z[off .+ (1:(nb + nc))]
	off += (nb + nc)
	s1 = z[off .+ (1:nc)]
	off += nc
	s2 = z[off .+ (1:nc)]
	off += nc
	return q2, γ1, b1, η1, s1, s2
end

function pack_z(model::ContactModel, env::Environment{<:World,LinearizedCone}, q2, γ1, b1, ψ1, η1)
	s1 = ϕ_func(model, env, q2)
	s2 = model.μ_world * γ1 .- E_func(model, env) * b1
	return [q2; γ1; b1; ψ1; η1; s1; s2]
end

function pack_z(model::ContactModel, env::Environment{<:World,NonlinearCone}, q2, γ1, b1, η1)
	s1 = ϕ_func(model, q2)
	s2 = model.μ_world .* γ1
	return [q2; γ1; b1; η1; s1; s2]
end

function z_initialize!(z, model::ContactModel, env::Environment{<:World,LinearizedCone}, q1)
	nq = model.dim.q
    z .= 1.0
    z[1:nq] = q1
end

function z_initialize!(z, model::ContactModel, env::Environment{<:World,NonlinearCone}, q1)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)
	ne = dim(env)

    z .= 0.1
    z[1:nq] = q1

	# second-order cone initializations
	z[nq + nc + nb .+ (1:(nb + nc))] = vcat([[1.0; 0.1 * ones(ne - 1)] for i = 1:model.dim.c]...)
	z[nq + nc + nb + nb + nc + nc .+ (1:nc)] .= 1.0

	return z
end

function z_warmstart!(z, model::ContactModel, env::Environment{<:World,LinearizedCone}, q, a, idx_ineq)
	z[1:model.dim.q] = q
	z[idx_ineq] .+= a * rand(length(idx_ineq))
	nothing
end

function z_warmstart!(z, model::ContactModel, env::Environment{<:World,NonlinearCone}, q, a, idx_ineq)
	@warn "warm start for second-order cone not implemented"
	z_initialize!(z, model, env, q)
	nothing
end

function θ_initialize!(θ, model::ContactModel, q0, q1, u, w, μ, h)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	off = 0
    θ[1:nq] = q0
	off += nq
    θ[off .+ (1:nq)] = q1
	off += nq
    θ[off .+ (1:nu)] = u
	off += nu
    θ[off .+ (1:nw)] = w
	off += nw
	θ[off .+ (1:1)] .= μ
	off += 1
    θ[off .+ (1:1)] .= h
end

function num_var(model::ContactModel, env::Environment{<:World,LinearizedCone})
	dim = model.dim
	nb = dim.c * friction_dim(env)
	dim.q + dim.c + nb + dim.c + nb + 2 * dim.c
end

function num_var(model::ContactModel, env::Environment{<:World,NonlinearCone})
	dim = model.dim
	nb = dim.c * friction_dim(env)
	dim.q + dim.c + nb + (nb + dim.c) + 2 * dim.c
end

function num_bil(model::ContactModel, env::Environment{<:World,LinearizedCone})
	dim = model.dim
	nb = dim.c * friction_dim(env)
	nb + 2 * dim.c
end

function num_data(model::ContactModel)
	dim = model.dim
	dim.q + dim.q + dim.u + dim.w + 1 + 1
end

inequality_indices(model::ContactModel, env::Environment{<:World,LinearizedCone}) = collect(model.dim.q .+ (1:(num_var(model, env) - model.dim.q)))
function inequality_indices(model::ContactModel, env::Environment{<:World,NonlinearCone})
	nb = model.dim.c * friction_dim(env)
	collect([(model.dim.q .+ (1:model.dim.c))...,
	         (model.dim.q + model.dim.c + nb + nb + model.dim.c .+ (1:model.dim.c))...])
end

soc_indices(model::ContactModel, env::Environment{<:World,LinearizedCone}) = Vector{Int}[]

function soc_indices(model::ContactModel, env::Environment{<:World,NonlinearCone})
	nq = model.dim.q
	nc = model.dim.c
	ne = dim(env)
	nf = friction_dim(env)
	nb = nc * nf

	b_idx = nq + nc .+ (1:nb)
	η_idx = nq + nc + nb .+ (1:(nb + nc))
	s2_idx = nq + nc + nb + nb + nc + nc .+ (1:nc)

	pr_idx = [[s2_idx[i]; b_idx[(i - 1) * nf .+ (1:nf)]] for i = 1:nc]
	du_idx = [[η_idx[(i - 1) * ne .+ (1:ne)]...] for i = 1:nc]

	[pr_idx..., du_idx...]
end

function ort_indices(model::ContactModel, env::Environment{<:World,LinearizedCone})
	ix, iy1, iy2 = linearization_var_index(model, env)
	pr_idx = iy1
	du_idx = iy2
	return [pr_idx, du_idx]
end

function ort_indices(model::ContactModel, env::Environment{<:World,NonlinearCone})
	nq = model.dim.q
	nc = model.dim.c
	ne = dim(env)
	nf = friction_dim(env)
	nb = nc * nf

	b_idx = nq + nc .+ (1:nb)
	η_idx = nq + nc + nb .+ (1:(nb + nc))
	s2_idx = nq + nc + nb + nb + nc + nc .+ (1:nc)

	pr_idx = [[s2_idx[i]; b_idx[(i - 1) * nf .+ (1:nf)]] for i = 1:nc]
	du_idx = [[η_idx[(i - 1) * ne .+ (1:ne)]...] for i = 1:nc]

	return [pr_idx, du_idx]
end

function E_func(model::ContactModel, env::Environment{<:World,LinearizedCone})
	nc = model.dim.c
	nb = nc * friction_dim(env)
	SMatrix{nc, nb}(kron(Diagonal(ones(nc)), ones(1, Int(nb / nc))))
end

function residual(model::ContactModel, env::Environment{<:World,LinearizedCone}, z, θ, κ)
	nc = model.dim.c
	nb = nc * friction_dim(env)
	nf = Int(nb / nc)
	np = dim(env)

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, env, z)

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
	q2, γ1, b1, η1, s1, s2 = unpack_z(model, env, z)

	ϕ = ϕ_func(model, env, q2)
	k = kinematics(model, q2)
	λ1 = contact_forces(model, env, γ1, b1, q2, k)
	vT = velocity_stack(model, env, q1, q2, k, h)

	ne = dim(env)

	[dynamics(model, env, h, q0, q1, u1, w1, λ1, q2);
	 s1 - ϕ;
	 vcat([η1[(i - 1) * ne .+ (2:ne)] - vT[(i - 1) * (ne - 1) .+ (1:(ne - 1))] for i = 1:model.dim.c]...);
	 s2 - μ[1] * γ1;
	 γ1 .* s1 .- κ;
	 vcat([second_order_cone_product(η1[(i - 1) * ne .+ (1:ne)], [s2[i]; b1[(i-1) * (ne - 1) .+ (1:(ne - 1))]]) - [κ; zeros(ne - 1)] for i = 1:model.dim.c]...)]
end

function ResidualMethods()
	function f()
		error("Not Implemented: use instantiate_residual!")
		return nothing
	end
	return ResidualMethods(fill(f, 3)...)
end
