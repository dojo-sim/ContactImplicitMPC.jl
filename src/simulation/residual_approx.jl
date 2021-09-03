
function res_con(model::ContactModel, env::Environment{<:World,LinearizedCone}, z, θ, κ)
	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	[s1 .- ϕ_func(model, env, q2);
	 s2 .- (μ[1] * γ1 .- E_func(model, env) * b1);
	 γ1 .* s1 .- κ;
	 b1 .* η1 .- κ;
	 ψ1 .* s2 .- κ]
end

function rz_approx!(s, rz, z, θ)
	rz .= 0.0

	model = s.model
	env = s.env

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	k = kinematics(model, q2)
	λ1 = contact_forces(model, env, γ1, b1, q2, k)
	vT = velocity_stack(model, env, q1, q2, k, h)

	nb = model.dim.c * friction_dim(env)

	# Dynamics
	rz[1:model.dim.q, 1:model.dim.q] = model.dyn.∂q2(h, q0, q1, u1, w1, λ1, q2) + s.con.dJ(λ1, q2)
	rz[1:model.dim.q, 1:(model.dim.q + model.dim.c + nb)] += s.con.dλ1(q2) * s.con.dcf(γ1, b1, q2, k)

	# Maximum dissipation
	rz[model.dim.q + model.dim.c .+ (1:nb), model.dim.q + model.dim.c + nb .+ (1:model.dim.c + nb)] = s.con.mdψη(vT, ψ1, η1)
	rz[model.dim.q + model.dim.c .+ (1:nb), 1:model.dim.q] += s.con.mdvs(vT, ψ1, η1) * s.con.vsq2(q1, q2, k, h)

	# Other constraints
	s.con.rcz(view(rz, collect([(model.dim.q .+ (1:model.dim.c))..., ((model.dim.q + model.dim.c + nb + 1):num_var(model, env))...]), :), z, θ)
end

function rθ_approx!(s, rθ, z, θ)
	rθ .= 0.0

	model = s.model
	env = s.env

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	k = kinematics(model, q2)
	λ1 = contact_forces(model, env, γ1, b1, q2, k)
	vT = velocity_stack(model, env, q1, q2, k, h)

	nc = model.dim.c
	nb = nc * friction_dim(env)

	# Dynamics
	idx = collect([(1:2model.dim.q + model.dim.u + model.dim.w)..., num_data(model)])
	rθ[1:model.dim.q, idx] = model.dyn.dθ(h, q0, q1, u1, w1, λ1, q2)

	# Maximum dissipation
	idx = collect([(model.dim.q .+ (1:model.dim.q))..., num_data(model)])
	rθ[model.dim.q + model.dim.c .+ (1:nb), idx] = s.con.vsq1h(q1, q2, k, h)

	# Other constraints
	s.con.rcθ(view(rθ, collect([(model.dim.q .+ (1:model.dim.c))..., ((model.dim.q + model.dim.c + nb + 1):num_var(model, env))...]), :), z, θ)
end
