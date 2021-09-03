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

# function unpack_z(model::ContactModel, env::Environment{<:World,LinearizedCone}, z)
# 	nq = model.dim.q
# 	nu = model.dim.u
# 	nc = model.dim.c
# 	nb = nc * friction_dim(env)
#
# 	# system variables
# 	off = 0
# 	q2 =  z[off .+ (1:nq)]
# 	off += nq
# 	γ1 =  z[off .+ (1:nc)]
# 	off += nc
# 	b1 =  z[off .+ (1:nb)]
# 	off += nb
# 	ψ1 =  z[off .+ (1:nc)]
# 	off += nc
# 	η1 =  z[off .+ (1:nb)]
# 	off += nb
# 	s1 = z[off .+ (1:nc)]
# 	off += nc
# 	s2 = z[off .+ (1:nc)]
# 	off += nc
# 	return q2, γ1, b1, ψ1, η1, s1, s2
# end

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
	s1 =  z[off .+ (1:nc)]
	off += nc
	η1 =  z[off .+ (1:nb)]
	off += nb
	s2 =  z[off .+ (1:nc)]
	off += nc
	return q2, γ1, b1, ψ1, s1, η1, s2
end

function unpack_z(model::ContactModel, env::Environment{<:World,NonlinearCone}, z)
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = nc * friction_dim(env)

	# # system variables
	# off = 0
	# q2 =  z[off .+ (1:nq)]
	# off += nq
	# γ1 =  z[off .+ (1:nc)]
	# off += nc
	# b1 =  z[off .+ (1:nb)]
	# off += nb
	# s1 = z[off .+ (1:nc)]
	# off += nc
	# η1 =  z[off .+ (1:(nb + nc))]
	# off += (nb + nc)
	# s2 = z[off .+ (1:nc)]
	# off += nc
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
	s1 =  z[off .+ (1:nc)]
	off += nc
	η1 =  z[off .+ (1:nb)]
	off += nb
	s2 =  z[off .+ (1:nc)]
	off += nc
	return q2, γ1, b1, ψ1, s1, η1, s2
end

function pack_z(model::ContactModel, env::Environment{<:World,LinearizedCone}, q2, γ1, b1, ψ1, η1)
	s1 = ϕ_func(model, env, q2)
	s2 = model.μ_world * γ1 .- E_func(model, env) * b1
	return [q2; γ1; b1; ψ1; s1; η1; s2]
end

function pack_z(model::ContactModel, env::Environment{<:World,NonlinearCone}, q2, γ1, b1, ψ1, η1)
	s1 = ϕ_func(model, q2)
	s2 = model.μ_world .* γ1
	return [q2; γ1; b1; ψ1; s1; η1; s2]
end

function z_initialize!(z, model::ContactModel, env::Environment{<:World,LinearizedCone}, q1)
	iq2 = index_q2(model, env, nquat = 0)
	z .= 1.0
	z[iq2] = q1
end

function z_initialize!(z, model::ContactModel, env::Environment{<:World,NonlinearCone}, q1)
	# nq = model.dim.q
	# nc = model.dim.c
	# nb = nc * friction_dim(env)
	# ne = dim(env)

	iq2 = index_q2(model, env, nquat = 0)
	ib1 = index_b1(model, env, nquat = 0)
	iψ1 = index_ψ1(model, env, nquat = 0)
	iη1 = index_η1(model, env, nquat = 0)
	is2 = index_s2(model, env, nquat = 0)

    z .= 0.1
    z[iq2] = q1

	# second-order cone initializations
	z[ib1] .= 0.1 # primal cones: vector part # TODO redundant
	z[iψ1] .= 1.0 # primal cones: scalar part
	z[iη1] .= 0.1 # dual cones: vector part # TODO redundant
	z[is2] .= 1.0 # dual cones: scalar part

	# z[nq + nc + nb .+ (1:(nb + nc))] = vcat([[1.0; 0.1 * ones(ne - 1)] for i = 1:model.dim.c]...)
	# z[nq + nc + nb + nb + nc + nc .+ (1:nc)] .= 1.0
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
	dim.q + dim.c + nb + dim.c + dim.c + nb + dim.c
end

function num_var(model::ContactModel, env::Environment{<:World,NonlinearCone})
	dim = model.dim
	nb = dim.c * friction_dim(env)
	dim.q + dim.c + nb + dim.c + dim.c + nb + dim.c
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


# function ort_indices(model::ContactModel, env::Environment{<:World,LinearizedCone})
# 	ix, iy1, iy2 = linearization_var_index(model, env)
# 	pr_idx = iy1
# 	du_idx = iy2
# 	return [pr_idx, du_idx]
# end
#
# function ort_indices(model::ContactModel, env::Environment{<:World,NonlinearCone})
# 	nq = model.dim.q
# 	nc = model.dim.c
# 	ne = dim(env)
# 	nf = friction_dim(env)
# 	nb = nc * nf
#
# 	b_idx = nq + nc .+ (1:nb)
# 	η_idx = nq + nc + nb .+ (1:(nb + nc))
# 	s2_idx = nq + nc + nb + nb + nc + nc .+ (1:nc)
#
# 	pr_idx = [[s2_idx[i]; b_idx[(i - 1) * nf .+ (1:nf)]] for i = 1:nc]
# 	du_idx = [[η_idx[(i - 1) * ne .+ (1:ne)]...] for i = 1:nc]
#
# 	return [pr_idx, du_idx]
# end

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
	 # vcat([η1[(i - 1) * ne .+ (2:ne)] - vT[(i - 1) * (ne - 1) .+ (1:(ne - 1))] for i = 1:model.dim.c]...);
	 vcat([η1[(i - 1) * nf .+ (1:nf)] - vT[(i - 1) * nf .+ (1:nf)] for i = 1:nc]...);
	 γ1 .* s1 .- κ;
	 # vcat([second_order_cone_product(η1[(i - 1) * ne .+ (1:ne)], [s2[i]; b1[(i-1) * (ne - 1) .+ (1:(ne - 1))]]) - [κ; zeros(ne - 1)] for i = 1:model.dim.c]...)
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
